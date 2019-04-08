## Changes 02/11, updates for 2nd data freeze
## - M2 multi-sample calling
## - batch-specific PON
## - coverage stats per gene
## - new M2 orientation filter
## - force calling IDH/TERT

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## CallPON
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## This rule uses Mutect2 in tumor-only mode (artifact detection mode) to detect ALL
## variants in a given non-tumor sample
## No germline resource in included because this would exclude these variants from the
## PON
## See: 
## https://gatkforums.broadinstitute.org/gatk/discussion/11136/how-to-call-somatic-mutations-using-gatk4-mutect2
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule callpon:
    input:
        bam = ancient("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam"),
        intervallist = lambda wildcards: "{dir}/{interval}/scattered.interval_list".format(dir = config["mutect2"]["wgs_scatterdir"], interval = wildcards.interval)
    output:
        vcf = temp("results/mutect2/callpon/{aliquot_barcode}/{aliquot_barcode}.{interval}.pon.vcf"),
        idx = temp("results/mutect2/callpon/{aliquot_barcode}/{aliquot_barcode}.{interval}.pon.vcf.idx")
    params:
        mem = CLUSTER_META["callpon"]["mem"],
        readgroup_sample_tag = lambda wildcards: manifest.getRGSampleTagByAliquot(wildcards.aliquot_barcode)
    conda:
        "../envs/gatk4.yaml"
    threads:
        CLUSTER_META["callpon"]["ppn"]
    log:
        "logs/mutect2/callpon/{aliquot_barcode}.{interval}.log"
    benchmark:
        "benchmarks/mutect2/callpon/{aliquot_barcode}.{interval}.txt"
    message:
        "Calling Mutect2 in tumor-only mode to build a panel of normals\n"
        "Aliquot: {wildcards.aliquot_barcode}\n"
        "Interval: {wildcards.interval}"
    shell:
        "gatk --java-options -Xmx{params.mem}g Mutect2 \
            -R {config[reference_fasta]} \
            -I {input.bam} \
            -L {input.intervallist} \
            --tumor-sample {params.readgroup_sample_tag} \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            -O {output.vcf} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## MergePON
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## This is an intermediate "gather" step. Because `callpon` is run in scatter mode, using
## and interval list and spread across 50 jobs for each normal sample, this step uses
## MergeVCF to merge all intermediate VCF files into one VCF per normal sample
## See: 
## https://gatkforums.broadinstitute.org/gatk/discussion/11136/how-to-call-somatic-mutations-using-gatk4-mutect2
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule mergepon:
    input:
        lambda wildcards: expand("results/mutect2/callpon/{aliquot_barcode}/{aliquot_barcode}.{interval}.pon.vcf", aliquot_barcode = wildcards.aliquot_barcode, interval = WGS_SCATTERLIST)
    output:
        vcf = "results/mutect2/mergepon/{aliquot_barcode}.pon.vcf",
        idx = "results/mutect2/mergepon/{aliquot_barcode}.pon.vcf.idx"
    params:
        mem = CLUSTER_META["mergepon"]["mem"],
        input_files = lambda _, input: " ".join(["-I " + s for s in input])
    conda:
        "../envs/gatk4.yaml"
    threads:
        CLUSTER_META["mergepon"]["ppn"]
    log:
        "logs/mutect2/mergepon/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/mutect2/mergepon/{aliquot_barcode}.txt"
    message:
        "Merging VCF files (PON)\n"
        "Aliquot: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g MergeVcfs \
            {params.input_files} \
            -O {output.vcf} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## CreatePON
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Because PON calls on normal samples were done individually for each normal, this step
## merges calls from all samples to build a PON
## Protected output, this is the final PON file
## See: 
## https://gatkforums.broadinstitute.org/gatk/discussion/11136/how-to-call-somatic-mutations-using-gatk4-mutect2
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule createpon:
    input:
        lambda wildcards: expand("results/mutect2/mergepon/{aliquot_barcode}.pon.vcf", aliquot_barcode = manifest.getPONAliquotsByBatch(manifest.parseBatch(wildcards.aliquot_batch)))
    output:
        vcf = protected("results/mutect2/pon/{aliquot_batch}.vcf"),
        idx = protected("results/mutect2/pon/{aliquot_batch}.vcf.idx")
    params:
        mem = CLUSTER_META["createpon"]["mem"],
        vcfs = lambda _, input: " ".join(["--vcfs=" + s for s in input])
    conda:
        "../envs/gatk4.yaml"
    threads:
        CLUSTER_META["createpon"]["ppn"]
    log:
        "logs/mutect2/createpon/{aliquot_batch}.log"
    benchmark:
        "benchmarks/mutect2/createpon/{aliquot_batch}.txt"
    message:
        "Creating panel of normals from multiple Mutect2 VCFs\n"
        "Batch: {wildcards.aliquot_batch}"
    shell:
        "gatk --java-options -Xmx{params.mem}g CreateSomaticPanelOfNormals \
            {params.vcfs} \
            --duplicate-sample-strategy THROW_ERROR \
            --output {output.vcf} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Call SNV
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## This step uses Mutect2 to call variants on a tumor-normal pair, using a panel-of-normals
## and a germline resource as reference
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## GATK parameters taken from "GATK4_SomaticSNVindel_worksheet.pdf"
## Extract SM tag from BAM file
## Code snippet taken from https://github.com/IARCbioinfo/BAM-tricks
## https://software.broadinstitute.org/gatk/documentation/article.php?id=3891
## Calculate af-of-alleles: 1/(2*123136) = 0.00000406055 (6 sign. digits)
## 06/09/2018
## Removed -disable-read-filter MateOnSameContigOrNoMappedMateReadFilter
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule callsnv:
    input:
        tumor = lambda wildcards: ancient(expand("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam", aliquot_barcode = manifest.getTumorByCase(wildcards.case_barcode))),
        normal = lambda wildcards: ancient(expand("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam", aliquot_barcode = manifest.getNormalByCase(wildcards.case_barcode))),
        pon = lambda wildcards: ancient("results/mutect2/pon/{aliquot_batch}.vcf".format(aliquot_batch = manifest.getBatchByCase(wildcards.case_barcode))),
        ponidx = lambda wildcards: ancient("results/mutect2/pon/{aliquot_batch}.vcf.idx".format(aliquot_batch = manifest.getBatchByCase(wildcards.case_barcode))),
        intervallist = lambda wildcards: ancient("{dir}/{interval}/scattered.interval_list".format(dir = config["mutect2"]["wgs_scatterdir"], interval = wildcards.interval)),
        orientation_priors = lambda wildcards: ancient(expand("results/mutect2/filterorientation/{aliquot_barcode}.priors.tsv", aliquot_barcode = manifest.getTumorByCase(wildcards.case_barcode))),
    output:
        vcf = temp("results/mutect2/m2vcf-scatter/{case_barcode}.{interval}.vcf"),
        idx = temp("results/mutect2/m2vcf-scatter/{case_barcode}.{interval}.vcf.idx"),
        bam = temp("results/mutect2/m2bam-scatter/{case_barcode}.{interval}.bam")
    params:
        mem = CLUSTER_META["callsnv"]["mem"],
        sample_paths = lambda _, input: " ".join(["-I " + s for s in input["tumor"] + input["normal"]]),
        priors_paths = lambda _, input: " ".join(["--orientation-bias-artifact-priors " + s for s in input["orientation_priors"]]),
        normal_sample_tags = lambda wildcards: " ".join(["--normal-sample " + s for s in [manifest.getRGSampleTagByAliquot(normalsample) for normalsample in manifest.getNormalByCase(wildcards.case_barcode)]])
    conda:
        "../envs/gatk4.yaml"
    threads:
        CLUSTER_META["callsnv"]["ppn"]
    log:
        "logs/mutect2/callsnv/{case_barcode}.{interval}.log"
    benchmark:
        "benchmarks/mutect2/callsnv/{case_barcode}.{interval}.txt"
    message:
        "Calling SNVs (Mutect2)\n"
        "Case: {wildcards.case_barcode}\n"
        "Interval: {wildcards.interval}"
    shell:
        "gatk --java-options -Xmx{params.mem}g Mutect2 \
            -R {config[reference_fasta]} \
            {params.sample_paths} \
            {params.normal_sample_tags} \
            {params.priors_paths} \
            -L {input.intervallist} \
            --panel-of-normals {input.pon} \
            --germline-resource {config[mutect2][gnomad_vcf]} \
            --genotyping-mode GENOTYPE_GIVEN_ALLELES \
            --genotype-filtered-alleles true \
            --alleles {config[mutect2][given_alleles]} \
            --dont-use-soft-clipped-bases true \
            --standard-min-confidence-threshold-for-calling 20 \
            -O {output.vcf} \
            -bamout {output.bam} \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Single-sample SNV
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Calls SNVs in a single aliquot rather than a patient
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule sscallsnv:
    input:
        tumor = lambda wildcards: ancient("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode))),
        normal = lambda wildcards: ancient("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam".format(aliquot_barcode = manifest.getNormal(wildcards.pair_barcode))),
        pon = lambda wildcards: ancient("results/mutect2/pon/{aliquot_batch}.vcf".format(aliquot_batch = manifest.getBatch(manifest.getTumor(wildcards.pair_barcode)))),
        ponidx = lambda wildcards: ancient("results/mutect2/pon/{aliquot_batch}.vcf.idx".format(aliquot_batch = manifest.getBatch(manifest.getTumor(wildcards.pair_barcode)))),
        intervallist = lambda wildcards: ancient("{dir}/{interval}/scattered.interval_list".format(dir = config["mutect2"]["wgs_scatterdir"], interval = wildcards.interval)),
        orientation_priors = lambda wildcards: ancient("results/mutect2/filterorientation/{aliquot_barcode}.priors.tsv".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode)))
    output:
        vcf = temp("results/mutect2/ssm2vcf-scatter/{pair_barcode}.{interval}.vcf"),
        idx = temp("results/mutect2/ssm2vcf-scatter/{pair_barcode}.{interval}.vcf.idx")
    params:
        mem = CLUSTER_META["sscallsnv"]["mem"],
        normal_sample_tag = lambda wildcards: manifest.getRGSampleTagByAliquot(manifest.getNormal(wildcards.pair_barcode))
    conda:
        "../envs/gatk4.yaml"
    threads:
        CLUSTER_META["sscallsnv"]["ppn"]
    log:
        "logs/mutect2/sscallsnv/{pair_barcode}.{interval}.log"
    benchmark:
        "benchmarks/mutect2/sscallsnv/{pair_barcode}.{interval}.txt"
    message:
        "Single Sample Calling SNVs (Mutect2)\n"
        "Pair: {wildcards.pair_barcode}\n"
        "Interval: {wildcards.interval}"
    shell:
        "gatk --java-options -Xmx{params.mem}g Mutect2 \
            -R {config[reference_fasta]} \
            -I {input.tumor} \
            -I {input.normal} \
            --orientation-bias-artifact-priors {input.orientation_priors} \
            -L {input.intervallist} \
            --normal-sample {params.normal_sample_tag} \
            --panel-of-normals {input.pon} \
            --germline-resource {config[mutect2][gnomad_vcf]} \
            --genotyping-mode GENOTYPE_GIVEN_ALLELES \
            --genotype-filtered-alleles true \
            --alleles {config[mutect2][given_alleles]} \
            --dont-use-soft-clipped-bases true \
            --standard-min-confidence-threshold-for-calling 20 \
            -O {output.vcf} \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Merge SNV
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## This is an intermediate "gather" step. Like `callpon`, `callsnv` is run in scatter mode, 
## using an interval list and spread across 50 jobs for each normal sample, this step uses
## MergeVCF to merge all intermediate VCF files into one VCF per normal sample
## See: 
## https://gatkforums.broadinstitute.org/gatk/discussion/11136/how-to-call-somatic-mutations-using-gatk4-mutect2
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule mergesnv:
    input:
        lambda wildcards: expand("results/mutect2/m2vcf-scatter/{case_barcode}.{interval}.vcf", case_barcode = wildcards.case_barcode, interval = WGS_SCATTERLIST)
    output:
        vcf = protected("results/mutect2/m2vcf/{case_barcode}.vcf"),
        idx = protected("results/mutect2/m2vcf/{case_barcode}.vcf.idx")
    params:
        mem = CLUSTER_META["mergesnv"]["mem"],
        input_files = lambda _, input: " ".join(["-I " + s for s in input])
    conda:
        "../envs/gatk4.yaml"
    threads:
        CLUSTER_META["mergesnv"]["ppn"]
    log:
        "logs/mutect2/mergesnv/{case_barcode}.log"
    benchmark:
        "benchmarks/mutect2/mergesnv/{case_barcode}.txt"
    message:
        "Merging VCF files (M2)\n"
        "Case: {wildcards.case_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g MergeVcfs \
            {params.input_files} \
            -O {output.vcf} \
            > {log} 2>&1"


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Merge SNV (single sample)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Single sample merge
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule ssmergesnv:
    input:
        lambda wildcards: expand("results/mutect2/ssm2vcf-scatter/{pair_barcode}.{interval}.vcf", pair_barcode = wildcards.pair_barcode, interval = WGS_SCATTERLIST)
    output:
        protected("results/mutect2/ssm2vcf/{pair_barcode}.vcf")
    params:
        mem = CLUSTER_META["ssmergesnv"]["mem"],
        input_files = lambda _, input: " ".join(["-I " + s for s in input])
    conda:
        "../envs/gatk4.yaml"
    threads:
        CLUSTER_META["ssmergesnv"]["ppn"]
    log:
        "logs/mutect2/ssmergesnv/{pair_barcode}.log"
    benchmark:
        "benchmarks/mutect2/ssmergesnv/{pair_barcode}.txt"
    message:
        "Single Sample Merging VCF files (M2)\n"
        "Pair: {wildcards.pair_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g MergeVcfs \
            {params.input_files} \
            -O {output} \
            > {log} 2>&1"

rule mergem2bam:
    input:
        lambda wildcards: expand("results/mutect2/m2bam-scatter/{case_barcode}.{interval}.bam", case_barcode = wildcards.case_barcode, interval = WGS_SCATTERLIST)
    output:
        protected("results/mutect2/m2bam/{case_barcode}.bam")
    params:
        mem = CLUSTER_META["mergem2bam"]["mem"],
        input_files = lambda _, input: " ".join(["-I " + s for s in input])
    conda:
        "../envs/gatk4.yaml"
    threads:
        CLUSTER_META["mergem2bam"]["ppn"]
    log:
        "logs/mutect2/mergem2bam/{case_barcode}.log"
    benchmark:
        "benchmarks/mutect2/mergem2bam/{case_barcode}.txt"
    message:
        "Merging BAM files (M2)\n"
        "Case: {wildcards.case_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g MergeSamFiles \
            {params.input_files} \
            -O {output} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Summarize read support for known variant sites
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule pileupsummaries:
    input:
        ancient("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam")
    output:
        temp("results/mutect2/pileupsummaries/{aliquot_barcode}.pileupsummaries.txt")
    params:
        mem = CLUSTER_META["pileupsummaries"]["mem"]
    conda:
        "../envs/gatk4.yaml"
    threads:
        CLUSTER_META["pileupsummaries"]["ppn"]
    log:
        "logs/mutect2/pileupsummaries/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/mutect2/pileupsummaries/{aliquot_barcode}.txt"
    message:
        "Generating pileupsummaries\n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g GetPileupSummaries \
            -I {input} \
            -V {config[mutect2][tiny_vcf]} \
            -L {config[mutect2][tiny_vcf]} \
            -O {output} \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Estimate contamination
## Input: pileup summaries table
## Output: contamination table
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule calculatecontamination:
    input:
        tumortable = lambda wildcards: "results/mutect2/pileupsummaries/{aliquot_barcode}.pileupsummaries.txt".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode)),
        normaltable = lambda wildcards: "results/mutect2/pileupsummaries/{aliquot_barcode}.pileupsummaries.txt".format(aliquot_barcode = manifest.getNormal(wildcards.pair_barcode)),
    output:
        cont = "results/mutect2/contamination/{pair_barcode}.contamination.txt",
        segs = "results/mutect2/contamination/{pair_barcode}.segmentation.txt"
    params:
        mem = CLUSTER_META["calculatecontamination"]["mem"]
    threads:
        CLUSTER_META["calculatecontamination"]["ppn"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/mutect2/calculatecontamination/{pair_barcode}.log"
    benchmark:
        "benchmarks/mutect2/calculatecontamination/{pair_barcode}.txt"
    message:
        "Computing contamination\n"
        "Pair: {wildcards.pair_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g CalculateContamination \
            -I {input.tumortable} \
            --matched-normal {input.normaltable} \
            --output {output.cont} \
            --tumor-segmentation {output.segs} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## If variants have not been filtered, filter, else done
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule filtermutect:
    input:
        vcf = ancient("results/mutect2/m2vcf/{case_barcode}.vcf"),
        tab = lambda wildcards: expand("results/mutect2/contamination/{pair_barcode}.contamination.txt", pair_barcode = manifest.getPairsByCase(wildcards.case_barcode)),
        seg = lambda wildcards: expand("results/mutect2/contamination/{pair_barcode}.segmentation.txt", pair_barcode = manifest.getPairsByCase(wildcards.case_barcode))
    output:
        stats = protected("results/mutect2/m2filter/{case_barcode}.filterstats.tsv"),
        vcf = protected("results/mutect2/m2filter/{case_barcode}.filtered.vcf.gz"),
        tbi = protected("results/mutect2/m2filter/{case_barcode}.filtered.vcf.gz.tbi")
    params:
        mem = CLUSTER_META["filtermutect"]["mem"],
        tseg = lambda _, input: " ".join(["--tumor-segmentation " + s for s in input["seg"]]),
        ttab = lambda _, input: " ".join(["--contamination-table " + s for s in input["tab"]]),
    threads:
        CLUSTER_META["filtermutect"]["ppn"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/mutect2/filtermutect/{case_barcode}.log"
    benchmark:
        "benchmarks/mutect2/filtermutect/{case_barcode}.txt"
    message:
        "Filtering Mutect2 calls\n"
        "Case: {wildcards.case_barcode}"
    shell:    
        "gatk --java-options -Xmx{params.mem}g FilterMutectCalls \
            -V {input.vcf} \
            {params.tseg} \
            {params.ttab} \
            --stats {output.stats} \
            -O {output.vcf} \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            > {log} 2>&1"


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Single sample filter
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule ssfiltermutect:
    input:
        vcf = ancient("results/mutect2/ssm2vcf/{pair_barcode}.vcf"),
        tab = "results/mutect2/contamination/{pair_barcode}.contamination.txt",
        seg = "results/mutect2/contamination/{pair_barcode}.segmentation.txt"
    output:
        stats = protected("results/mutect2/ssm2filter/{pair_barcode}.filterstats.tsv"),
        vcf = protected("results/mutect2/ssm2filter/{pair_barcode}.filtered.vcf.gz"),
        tbi = protected("results/mutect2/ssm2filter/{pair_barcode}.filtered.vcf.gz.tbi")
    params:
        mem = CLUSTER_META["ssfiltermutect"]["mem"]
    threads:
        CLUSTER_META["ssfiltermutect"]["ppn"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/mutect2/ssfiltermutect/{pair_barcode}.log"
    benchmark:
        "benchmarks/mutect2/ssfiltermutect/{pair_barcode}.txt"
    message:
        "Single-Sample Filtering Mutect2 calls\n"
        "Pair: {wildcards.pair_barcode}"
    shell:    
        "gatk --java-options -Xmx{params.mem}g FilterMutectCalls \
            -V {input.vcf} \
            --tumor-segmentation {input.seg} \
            --contamination-table {input.tab} \
            --stats {output.stats} \
            -O {output.vcf} \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Collect metrics on sequencing context artifacts
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule collectartifacts:
    input:
        ancient("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam")
    output:
        tab = temp("results/mutect2/artifacts/{aliquot_barcode}.alt.tsv"),
        ref = temp("results/mutect2/artifacts/{aliquot_barcode}.ref.metrics"),
        alt = temp("results/mutect2/artifacts/{aliquot_barcode}.alt.metrics")
    params:
        prefix = "results/mutect2/artifacts/{aliquot_barcode}",
        mem = CLUSTER_META["collectartifacts"]["mem"]
    threads:
        CLUSTER_META["collectartifacts"]["ppn"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/mutect2/collectartifacts/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/mutect2/collectartifacts/{aliquot_barcode}.txt"
    message:
        "Collecting sequencing artifact metrics\n"
        "Aliquot: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g CollectF1R2Counts \
            -R {config[reference_fasta]} \
            -I {input} \
            -alt-table {output.tab} \
            -ref-hist {output.ref} \
            -alt-hist {output.alt} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Filter by orientation bias
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule filterorientation:
    input:
        tab = "results/mutect2/artifacts/{aliquot_barcode}.alt.tsv",
        ref = "results/mutect2/artifacts/{aliquot_barcode}.ref.metrics",
        alt = "results/mutect2/artifacts/{aliquot_barcode}.alt.metrics"
    output:
        "results/mutect2/filterorientation/{aliquot_barcode}.priors.tsv"
    params:
        mem = CLUSTER_META["filterorientation"]["mem"]
    threads:
        CLUSTER_META["filterorientation"]["ppn"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/mutect2/filterorientation/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/mutect2/filterorientation/{aliquot_barcode}.txt"
    message:
        "Calculating orientation bias\n"
        "Aliquot: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g LearnReadOrientationModel \
            -alt-table {input.tab} \
            -ref-hist {input.ref} \
            -alt-hist {input.alt} \
            -O {output} \
            > {log} 2>&1"

## END ##