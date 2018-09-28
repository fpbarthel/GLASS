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
        bam = "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam",
        intervallist = lambda wildcards: "{dir}/{interval}/scattered.interval_list".format(dir = config["wgs_scatterdir"], interval = wildcards.interval)
    output:
        temp("results/mutect2/callpon/{aliquot_barcode}/{aliquot_barcode}.{interval}.pon.vcf")
    params:
        mem = CLUSTER_META["callpon"]["mem"]
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
        "TUMOR_SM=`samtools view -H {input} | grep '^@RG' | sed \"s/.*SM:\\([^\\t]*\\).*/\\1/g\" | uniq`; \
        gatk --java-options -Xmx{params.mem}g Mutect2 \
            -R {config[reference_fasta]} \
            -I {input.bam} \
            -L {input.intervallist} \
            --tumor-sample $TUMOR_SM \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            -O {output} \
            > {log} 2>&1"

            #        TUMOR_SM=${TUMOR_SM:-\"TUMOR\"}; \

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
        temp("results/mutect2/mergepon/{aliquot_barcode}.pon.vcf")
    params:
        mem = CLUSTER_META["mergepon"]["mem"]
    threads:
        CLUSTER_META["mergepon"]["ppn"]
    log:
        "logs/mutect2/mergepon/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/mutect2/mergepon/{aliquot_barcode}.txt"
    message:
        "Merging VCF files (PON)\n"
        "Aliquot: {wildcards.aliquot_barcode}"
    run:
        input_cat = " ".join(["-I " + s for s in input])
        shell("gatk --java-options -Xmx{params.mem}g MergeVcfs \
            {input_cat} \
            -O {output} \
            > {log} 2>&1")

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
        lambda wildcards: expand("results/mutect2/mergepon/{aliquot_barcode}.pon.vcf", aliquot_barcode = manifest.getPONAliquots())
    output:
        protected("results/mutect2/pon/pon.vcf")
    params:
        mem = CLUSTER_META["createpon"]["mem"]
    threads:
        CLUSTER_META["createpon"]["ppn"]
    log:
        "logs/mutect2/createpon/createpon.log"
    benchmark:
        "benchmarks/mutect2/createpon/createpon.txt"
    message:
        "Creating panel of normals from multiple Mutect2 VCFs"
    run:
        vcfs = " ".join(["--vcfs=" + s for s in input])
        shell("gatk --java-options -Xmx{params.mem}g CreateSomaticPanelOfNormals \
            {vcfs} \
            --duplicate-sample-strategy THROW_ERROR \
            --output {output} \
            > {log} 2>&1")

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
        tumor = lambda wildcards: "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam".format(aliquot_barcode = manifest.getTumor(wildcards.pair_id)),
        normal = lambda wildcards: "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam".format(aliquot_barcode = manifest.getNormal(wildcards.pair_id)),
        pon = "results/mutect2/pon/pon.vcf",
        intervallist = lambda wildcards: "{dir}/{interval}/scattered.interval_list".format(dir = config["wgs_scatterdir"], interval = wildcards.interval)
    output:
        vcf = temp("results/mutect2/m2vcf-scatter/{pair_id}.{interval}.vcf"),
        bam = "results/mutect2/m2bam/{pair_id}.{interval}.bam"
    params:
        mem = CLUSTER_META["callsnv"]["mem"]
    conda:
        "../envs/gatk4.yaml"
    threads:
        CLUSTER_META["callsnv"]["ppn"]
    log:
        "logs/mutect2/callsnv/{pair_id}.{interval}.log"
    benchmark:
        "benchmarks/mutect2/callsnv/{pair_id}.{interval}.txt"
    message:
        "Calling SNVs (Mutect2)\n"
        "Pair: {wildcards.pair_id}\n"
        "Interval: {wildcards.interval}"
    shell:
        "TEST_NAM=`samtools view -H {input.tumor} | grep '^@RG' | sed \"s/.*SM:\\([^\\t]*\\).*/\\1/g\" | uniq`;"
        "CTRL_NAM=`samtools view -H {input.normal} | grep '^@RG' | sed \"s/.*SM:\\([^\\t]*\\).*/\\1/g\" | uniq`;"
        #"TEST_NAM=${TEST_NAM:-\"TUMOR\"};"
        #"CTRL_NAM=${CTRL_NAM:-\"NORMAL\"};"
        "gatk --java-options -Xmx{params.mem}g Mutect2 \
            -R {config[reference_fasta]} \
            -I {input.tumor} \
            -I {input.normal} \
            -L {input.intervallist} \
            --tumor-sample $TEST_NAM \
            --normal-sample $CTRL_NAM \
            --panel-of-normals {input.pon} \
            --germline-resource {config[gnomad_vcf]} \
            --af-of-alleles-not-in-resource {config[af_of_alleles_not_in_resource]} \
            --dont-use-soft-clipped-bases true \
            --standard-min-confidence-threshold-for-calling 20 \
            -O {output.vcf} \
            -bamout {output.bam} \
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
        lambda wildcards: expand("results/mutect2/m2vcf-scatter/{pair_id}.{interval}.vcf", pair_id = wildcards.pair_id, interval = WGS_SCATTERLIST)
    output:
        protected("results/mutect2/m2vcf/{pair_id}.vcf")
    params:
        mem = CLUSTER_META["mergesnv"]["mem"]
    threads:
        CLUSTER_META["mergesnv"]["ppn"]
    log:
        "logs/mutect2/mergesnv/{pair_id}.log"
    benchmark:
        "benchmarks/mutect2/mergesnv/{pair_id}.txt"
    message:
        "Merging VCF files (M2)\n"
        "Pair: {wildcards.pair_id}"
    run:
        input_cat = " ".join(["-I " + s for s in input])
        shell("gatk --java-options -Xmx{params.mem}g MergeVcfs \
            {input_cat} \
            -O {output} \
            > {log} 2>&1")

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Summarize read support for known variant sites
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# @sbamin What is small_exac_common_3_b37.vcf.gz? Haven't read further but search only gives me this link: https://crc.pitt.edu/variantcalling
# @barthf This is used in broad pipelines

rule pileupsummaries:
    input:
        "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam"
    output:
        "results/mutect2/pileupsummaries/{aliquot_barcode}.pileupsummaries.txt"
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
            -V {config[tiny_vcf]} \
            -O {output} \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Estimate contamination
## Input: pileup summaries table
## Output: contamination table
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# @sbamin We don't necessarily need contamination table based on variants across all chrs, and 2-3 large chrs are good enough.
# Idea is to query sufficient population level variants to estimate normal contamination of tumor sample.
# @barthf this is fast enough that I'm just gonna leave this for now

rule calculatecontamination:
    input:
        tumortable = lambda wildcards: "results/mutect2/pileupsummaries/{aliquot_barcode}.pileupsummaries.txt".format(aliquot_barcode = manifest.getTumor(wildcards.pair_id)),
        normaltable = lambda wildcards: "results/mutect2/pileupsummaries/{aliquot_barcode}.pileupsummaries.txt".format(aliquot_barcode = manifest.getNormal(wildcards.pair_id)),
    output:
        "results/mutect2/contamination/{pair_id}.contamination.txt"
    params:
        mem = CLUSTER_META["calculatecontamination"]["mem"]
    threads:
        CLUSTER_META["calculatecontamination"]["ppn"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/mutect2/calculatecontamination/{pair_id}.log"
    benchmark:
        "benchmarks/mutect2/calculatecontamination/{pair_id}.txt"
    message:
        "Computing contamination\n"
        "Pair: {wildcards.pair_id}"
    shell:
        "gatk --java-options -Xmx{params.mem}g CalculateContamination \
            -I {input.tumortable} \
            --matched-normal {input.normaltable} \
            -O {output} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## If variants have not been filtered, filter, else done
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule filtermutect:
    input:
        vcf = "results/mutect2/m2vcf/{pair_id}.vcf",
        tab = "results/mutect2/contamination/{pair_id}.contamination.txt"
    output:
        temp("results/mutect2/m2filter/{pair_id}.filtered.vcf")
    params:
        mem = CLUSTER_META["filtermutect"]["mem"]
    threads:
        CLUSTER_META["filtermutect"]["ppn"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/mutect2/filtermutect/{pair_id}.log"
    benchmark:
        "benchmarks/mutect2/filtermutect/{pair_id}.txt"
    message:
        "Filtering Mutect2 calls\n"
        "Pair: {wildcards.pair_id}"
    shell:    
        "gatk --java-options -Xmx{params.mem}g FilterMutectCalls \
            -V {input.vcf} \
            --contamination-table {input.tab} \
            -O {output} \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Collect metrics on sequencing context artifacts
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule collectartifacts:
    input:
        "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam"
    output:
        "results/mutect2/artifacts/{aliquot_barcode}.pre_adapter_detail_metrics.txt"
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
        "Sample: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g CollectSequencingArtifactMetrics \
            -I {input} \
            -O {params.prefix} \
            --FILE_EXTENSION \".txt\" \
            -R {config[reference_fasta]} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Filter by orientation bias
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule filterorientation:
    input:
        art = lambda wildcards: "results/mutect2/artifacts/{aliquot_barcode}.pre_adapter_detail_metrics.txt".format(aliquot_barcode = manifest.getTumor(wildcards.pair_id)),
        vcf = "results/mutect2/m2filter/{pair_id}.filtered.vcf"
    output:
        "results/mutect2/final/{pair_id}.final.vcf"
    params:
        mem = CLUSTER_META["filterorientation"]["mem"]
    threads:
        CLUSTER_META["filterorientation"]["ppn"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/mutect2/filterorientation/{pair_id}.log"
    benchmark:
        "benchmarks/mutect2/filterorientation/{pair_id}.txt"
    message:
        "Filtering Mutect2 calls by orientation bias\n"
        "Pair: {wildcards.pair_id}"
    shell:
        "gatk --java-options -Xmx{params.mem}g FilterByOrientationBias \
            -AM \"G/T\" \
            -AM \"C/T\" \
            -V {input.vcf} \
            -P {input.art} \
            -O {output} \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            > {log} 2>&1"

## END ##
