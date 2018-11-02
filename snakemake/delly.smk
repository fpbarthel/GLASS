## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Delly
## At least one tumor sample and a matched control sample are required for SV discovery
## See: https://github.com/dellytools/delly
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule delly_call:
    input:
        tumor = lambda wildcards: "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode)),
        normal = lambda wildcards: "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam".format(aliquot_barcode = manifest.getNormal(wildcards.pair_barcode))
    output:
        bcf = "results/delly/call/{pair_barcode}.bcf"
#        vcf = "results/delly/call/{pair_barcode}.vcf.gz"
    params:
#        vcftmp = "results/delly/call/{pair_barcode}.vcf",
        mem = CLUSTER_META["delly_call"]["mem"]
    threads:
        CLUSTER_META["delly_call"]["ppn"]
    conda:
        "../envs/delly.yaml"
    log:
        "logs/delly/call/{pair_barcode}.log"
    benchmark:
        "benchmarks/delly/call/{pair_barcode}.txt"
    message:
        "Calling DELLY on tumor/normal pair\n"
        "Pair: {wildcards.pair_barcode}"
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
        "results/delly/call/{pair_barcode}.bcf"
    output:
    	tsv = "results/delly/filter/{pair_barcode}.samples.tsv",
        bcf = "results/delly/filter/{pair_barcode}.prefilt.bcf"
    params:
        mem = CLUSTER_META["delly_prefilter"]["mem"],
        tumor_sm = lambda wildcards: manifest.getRGSampleTagByAliquot(manifest.getTumor(wildcards.pair_barcode)),
        normal_sm = lambda wildcards: manifest.getRGSampleTagByAliquot(manifest.getNormal(wildcards.pair_barcode))
    threads:
        CLUSTER_META["delly_prefilter"]["ppn"]
    conda:
        "../envs/delly.yaml"
    log:
        "logs/delly/filter/{pair_barcode}.log"
    benchmark:
        "benchmarks/delly/filter/{pair_barcode}.txt"
    message:
        "Pre-filtering somatic DELLY calls\n"
        "Pair: {wildcards.pair_barcode}"
    shell:
        "export OMP_NUM_THREADS=2; "
        "printf '{params.tumor_sm}\\ttumor\\n{params.normal_sm}\\tcontrol\\n' > {output.tsv}; "
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
#         bcf = "results/delly/filter/{pair_barcode}.prefilt.bcf",
#         bam = "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam"
#     output:
#         bcf = "results/delly/gtcontrol/{pair_barcode}.prefilt.gtcontrol.bcf"
#     params:
#         mem = CLUSTER_META["delly_genotype_controls"]["mem"]
#     threads:
#         CLUSTER_META["delly_genotype_controls"]["ppn"]
#     conda:
#         "../envs/delly.yaml"
#     log:
#         "logs/delly/gtcontrol/{pair_barcode}.log"
#     benchmark:
#         "benchmarks/delly/gtcontrol/{pair_barcode}.txt"
#     message:
#         "Genotype found variants across a set of controls\n"
#         "Pair: {wildcards.pair_barcode}\n"
#         "Control sample: {wildcards.aliquot_barcode}"
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
#         bcf = "results/delly/filter/{pair_barcode}.prefilt.bcf",
#         bam = "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam"
#     output:
#         bcf = "results/delly/filter/{pair_barcode}.prefilt.bcf"
#     params:
#         mem = CLUSTER_META["delly_filter"]["mem"]
#     threads:
#         CLUSTER_META["delly_filter"]["ppn"]
#     conda:
#         "../envs/delly.yaml"
#     log:
#         "logs/delly/filter/{pair_barcode}.log"
#     benchmark:
#         "benchmarks/delly/filter/{pair_barcode}.txt"
#     message:
#         "Pre-filtering somatic DELLY calls\n"
#         "Pair: {wildcards.pair_barcode}"
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