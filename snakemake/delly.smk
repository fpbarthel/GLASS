## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Delly
## At least one tumor sample and a matched control sample are required for SV discovery
## See: https://github.com/dellytools/delly
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule delly_call:
    input:
        tumor = lambda wildcards: "results/align/bqsr/{aliquot_id}.realn.mdup.bqsr.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["tumor_aliquot_id"]),
        normal = lambda wildcards: "results/align/bqsr/{aliquot_id}.realn.mdup.bqsr.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["normal_aliquot_id"])
    output:
        bcf = "results/delly/call/{pair_id}.bcf"
#        vcf = "results/delly/call/{pair_id}.vcf.gz"
    params:
#        vcftmp = "results/delly/call/{pair_id}.vcf",
        mem = CLUSTER_META["delly_call"]["mem"]
    threads:
        CLUSTER_META["delly_call"]["ppn"]
    conda:
        "../envs/delly.yaml"
    log:
        "logs/delly/call/{pair_id}.log"
    benchmark:
        "benchmarks/delly/call/{pair_id}.txt"
    message:
        "Calling DELLY on tumor/normal pair\n"
        "Pair: {wildcards.pair_id}"
    shell:
        "export OMP_NUM_THREADS=2; "
        "delly call \
            -n \
            -x {config[svmask_delly]} \
            -o {output.bcf} \
            -g {config[reference_fasta]} \
            {input.tumor} \
            {input.normal} \
            > {log} 2>&1; "
#        "bcftools view {output.bcf} > {params.vcftmp} && \
#            bgzip -i {params.vcftmp} && \
#            bcftools index -t {output.vcf}"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Pre-filter delly results
## Somatic pre-filtering requires a tab-delimited sample description file where the first 
## column is the sample id (as in the VCF/BCF file) and the second column is either tumor 
## or control.
## See: https://github.com/dellytools/delly
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule delly_prefilter:
    input:
        "results/delly/call/{pair_id}.bcf"
    output:
        bcf = "results/delly/filter/{pair_id}.prefilt.bcf"
    params:
        mem = CLUSTER_META["delly_filter"]["mem"]
    threads:
        CLUSTER_META["delly_filter"]["ppn"]
    conda:
        "../envs/delly.yaml"
    log:
        "logs/delly/filter/{pair_id}.log"
    benchmark:
        "benchmarks/delly/filter/{pair_id}.txt"
    message:
        "Pre-filtering somatic DELLY calls\n"
        "Pair: {wildcards.pair_id}"
    shell:
        "export OMP_NUM_THREADS=2; "
        "printf '{params.tumor_id}\ttumor\n{params.normal_id}\tcontrol\n' > {output.tsv}; "
        "delly filter \
            -f somatic \
            -o {output.bcf} \
            -s {output.tsv} \
            {input} \
            > {log} 2>&1; "

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Genotype found variants across a panel of controls
## Genotype pre-filtered somatic sites across a larger panel of control samples to efficiently 
## filter false postives and germline SVs. For performance reasons, this can be run in parallel 
## for each sample of the control panel and you may want to combine multiple pre-filtered 
## somatic site lists from multiple tumor samples.
## See: https://github.com/dellytools/delly
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# rule delly_genotype_controls:
#     input:
#         bcf = "results/delly/filter/{pair_id}.prefilt.bcf",
#         bam = "results/align/bqsr/{aliquot_id}.realn.mdup.bqsr.bam"
#     output:
#         bcf = "results/delly/gtcontrol/{pair_id}.prefilt.gtcontrol.bcf"
#     params:
#         mem = CLUSTER_META["delly_genotype_controls"]["mem"]
#     threads:
#         CLUSTER_META["delly_genotype_controls"]["ppn"]
#     conda:
#         "../envs/delly.yaml"
#     log:
#         "logs/delly/gtcontrol/{pair_id}.log"
#     benchmark:
#         "benchmarks/delly/gtcontrol/{pair_id}.txt"
#     message:
#         "Genotype found variants across a set of controls\n"
#         "Pair: {wildcards.pair_id}\n"
#         "Control sample: {wildcards.aliquot_id}"
#     shell:
#         "printf '{params.tumor_id}\ttumor\n{params.normal_id}\tcontrol\n' > {output.tsv}; "
#         "delly call \
#             -v {input.bcf} \
#             -g {config[reference_fasta]} \
#             -o {output.bcf} \
#             -s {output.tsv} \
#             {input} \
#             > {log} 2>&1; "


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Merge genotype BCF files
## See: https://github.com/dellytools/dellys
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# rule delly_merge_genotypes:
#     input:
#         bcf = "results/delly/filter/{pair_id}.prefilt.bcf",
#         bam = "results/align/bqsr/{aliquot_id}.realn.mdup.bqsr.bam"
#     output:
#         bcf = "results/delly/filter/{pair_id}.prefilt.bcf"
#     params:
#         mem = CLUSTER_META["delly_filter"]["mem"]
#     threads:
#         CLUSTER_META["delly_filter"]["ppn"]
#     conda:
#         "../envs/delly.yaml"
#     log:
#         "logs/delly/filter/{pair_id}.log"
#     benchmark:
#         "benchmarks/delly/filter/{pair_id}.txt"
#     message:
#         "Pre-filtering somatic DELLY calls\n"
#         "Pair: {wildcards.pair_id}"
#     shell:
#         "printf '{params.tumor_id}\ttumor\n{params.normal_id}\tcontrol\n' > {output.tsv}; "
#         "delly filter \
#             -f somatic \
#             -o {output.bcf} \
#             -s {output.tsv} \
#             {input} \
#             > {log} 2>&1; "


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Post-filter
## Post-filter for somatic SVs using all control samples.
## See: https://github.com/dellytools/delly
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 


## END ##