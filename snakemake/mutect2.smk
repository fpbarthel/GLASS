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
        bam = "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam",
        intervallist = lambda wildcards: "{dir}/{interval}/scattered.interval_list".format(dir=config["wgs_scatterdir"], interval=wildcards.interval)
    output:
        temp("results/callpon/{batch}/{aliquot_id}/{aliquot_id}.{interval}.pon.vcf")
    params:
        mem = CLUSTER_META["callpon"]["mem"]
    threads:
        CLUSTER_META["callpon"]["ppn"]
    log:
        "logs/callpon/{batch}.{aliquot_id}.{interval}.log"
    benchmark:
        "benchmarks/callpon/{batch}.{aliquot_id}.{interval}.txt"
    message:
        "Calling Mutect2 in tumor-only mode to build a panel of normals\n"
        "Batch: {wildcards.batch}\n"
        "Sample: {wildcards.aliquot_id}\n"
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
        lambda wildcards: expand("results/callpon/{batch}/{aliquot_id}/{aliquot_id}.{interval}.pon.vcf", batch=wildcards.batch, aliquot_id=wildcards.aliquot_id, interval=WGS_SCATTERLIST)
    output:
        temp("results/mergepon/{batch}/{aliquot_id}.pon.vcf")
    params:
        mem = CLUSTER_META["mergepon"]["mem"]
    threads:
        CLUSTER_META["mergepon"]["ppn"]
    log:
        "logs/mergepon/{aliquot_id}.log"
    benchmark:
        "benchmarks/mergepon/{aliquot_id}.txt"
    message:
        "Merging VCF files (PON)\n"
        "Sample: {wildcards.aliquot_id}"
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
## merges calls from all samples from a given batch to build a confident PON
## Protected output, this is the final batch-specifc PON file
## See: 
## https://gatkforums.broadinstitute.org/gatk/discussion/11136/how-to-call-somatic-mutations-using-gatk4-mutect2
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule createpon:
    input:
        lambda wildcards: expand("results/mergepon/{batch}/{aliquot_id}.pon.vcf", batch=wildcards.batch, aliquot_id=BATCH_TO_NORMAL[wildcards.batch])
    output:
        protected("results/pon/{batch}.pon.vcf")
    params:
        dir = "results/qc/multiqc",
        mem = CLUSTER_META["createpon"]["mem"]
    threads:
        CLUSTER_META["createpon"]["ppn"]
    log:
        "logs/createpon/{batch}.createpon.log"
    benchmark:
        "benchmarks/createpon/{batch}.createpon.txt"
    message:
        "Creating panel of normals from multiple Mutect2 VCFs\n"
        "Batch: {wildcards.batch}"
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
        tumor = lambda wildcards: "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["tumor_aliquot_id"]),
        normal = lambda wildcards: "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["normal_aliquot_id"]),
        pon = lambda wildcards: "results/pon/{batch}.pon.vcf".format(batch=PAIRS_DICT[wildcards.pair_id]["project_id"]),
        intervallist = lambda wildcards: "{dir}/{interval}/scattered.interval_list".format(dir=config["wgs_scatterdir"], interval=wildcards.interval)
    output:
        vcf = temp("results/m2vcf-scatter/{pair_id}.{interval}.vcf"),
        bam = "results/m2bam/{pair_id}.{interval}.bam"
    params:
        mem = CLUSTER_META["callsnv"]["mem"]
    threads:
        CLUSTER_META["callsnv"]["ppn"]
    log:
        "logs/callsnv/{pair_id}.log"
    benchmark:
        "benchmarks/callsnv/{pair_id}.txt"
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
        lambda wildcards: expand("results/m2vcf-scatter/{pair_id}.{interval}.vcf", pair_id=wildcards.pair_id, interval=WGS_SCATTERLIST)
    output:
        protected("results/m2vcf/{pair_id}.vcf")
    params:
        mem = CLUSTER_META["mergesnv"]["mem"]
    threads:
        CLUSTER_META["mergesnv"]["ppn"]
    log:
        "logs/mergesnv/{pair_id}.log"
    benchmark:
        "benchmarks/mergesnv/{pair_id}.txt"
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
        "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam"
    output:
        "results/pileupsummaries/{aliquot_id}.pileupsummaries.txt"
    params:
        mem = CLUSTER_META["pileupsummaries"]["mem"]
    threads:
        CLUSTER_META["pileupsummaries"]["ppn"]
    log:
        "logs/pileupsummaries/{aliquot_id}.log"
    benchmark:
        "benchmarks/pileupsummaries/{aliquot_id}.txt"
    message:
        "Generating pileupsummaries\n"
        "Sample: {wildcards.aliquot_id}"
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
        tumortable = lambda wildcards: "results/pileupsummaries/{aliquot_id}.pileupsummaries.txt".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["tumor_aliquot_id"]),
        normaltable = lambda wildcards: "results/pileupsummaries/{aliquot_id}.pileupsummaries.txt".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["normal_aliquot_id"]),
    output:
        "results/contamination/{pair_id}.contamination.txt"
    params:
        mem = CLUSTER_META["calculatecontamination"]["mem"]
    threads:
        CLUSTER_META["calculatecontamination"]["ppn"]
    log:
        "logs/calculatecontamination/{pair_id}.log"
    benchmark:
        "benchmarks/calculatecontamination/{pair_id}.txt"
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
        vcf = "results/m2vcf/{pair_id}.vcf",
        tab = "results/contamination/{pair_id}.contamination.txt"
    output:
        temp("results/m2filter/{pair_id}.filtered.vcf")
    params:
        mem = CLUSTER_META["filtermutect"]["mem"]
    threads:
        CLUSTER_META["filtermutect"]["ppn"]
    log:
        "logs/filtermutect/{pair_id}.log"
    benchmark:
        "benchmarks/filtermutect/{pair_id}.txt"
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
        "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam"
    output:
        "results/artifacts/{aliquot_id}.pre_adapter_detail_metrics.txt"
    params:
        prefix = "results/artifacts/{aliquot_id}",
        mem = CLUSTER_META["collectartifacts"]["mem"]
    threads:
        CLUSTER_META["collectartifacts"]["ppn"]
    log:
        "logs/collectartifacts/{aliquot_id}.log"
    benchmark:
        "benchmarks/collectartifacts/{aliquot_id}.txt"
    message:
        "Collecting sequencing artifact metrics\n"
        "Sample: {wildcards.aliquot_id}"
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
        art = lambda wildcards: "results/artifacts/{aliquot_id}.pre_adapter_detail_metrics.txt".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["tumor_aliquot_id"]),
        vcf = "results/m2filter/{pair_id}.filtered.vcf"
    output:
        "results/m2filter/{pair_id}.filtered2.vcf"
    params:
        mem = CLUSTER_META["filterorientation"]["mem"]
    threads:
        CLUSTER_META["filterorientation"]["ppn"]
    log:
        "logs/filterorientation/{pair_id}.log"
    benchmark:
        "benchmarks/filterorientation/{pair_id}.txt"
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

# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# ## VariantFiltration
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# ## Hard filter variants, removes "bad" variants from pre-VEP VCF
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# rule hardfilter:
#     input:
#         "results/m2filter/{pair_id}.filtered2.vcf"
#     output:
#         "results/m2filter/{pair_id}.filtered3.vcf"
#     params:
#         mem = CLUSTER_META["hardfilter"]["mem"]
#     threads:
#         CLUSTER_META["hardfilter"]["ppn"]
#     log:
#         "logs/hardfilter/{pair_id}.log"
#     benchmark:
#         "benchmarks/hardfilter/{pair_id}.txt"
#     message:
#         "Hard-filtering soft-filtered Mutect2 calls\n"
#         "Pair: {wildcards.pair_id}"
#     shell:
#         "gatk --java-options -Xmx{params.mem}g VariantFiltration \
#             -AM \"G/T\" \
#             -AM \"C/T\" \
#             -V {input.vcf} \
#             -P {input.art} \
#             -O {output} \
#             --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
#             > {log} 2>&1"

# @sbamin I suggest that we should put hold at VEP step until we finalize consensus calls. 
# Idea is too get how many filtered calls we have and with what level of confidence (based on consensus calls from other callers).
# To me, VEP with extended anntoations is of little value until vcfs are at freeze level.
# Also, we should have uniform processing past alignment step across both canine and human (adult + pediatric) life history projects.
# A few callers and filtering steps that we already have run on canine are at https://github.com/TheJacksonLaboratory/hourglass/issues/9.
# If we are in need of base annotations, like exonic vs promoter vs non-coding, that should be quick enough from VEP, 
# but I suggest to have a hold on extended annotations using vcf2maf as I guess we will end up doing those iteratively until we freeze vcfs.
