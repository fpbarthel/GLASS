## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Extracting read mappings from BAM file
## See: https://github.com/abyzovlab/CNVnator
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule cnvnator_tree:
    input:
        "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam"
    output:
        "results/cnvnator/tree/{aliquot_id}.root"
    params:
        mem = CLUSTER_META["cnvnator_tree"]["mem"]
    threads:
        CLUSTER_META["cnvnator_tree"]["ppn"]
    log:
        "logs/cnvnator/tree/{aliquot_id}.log"
    benchmark:
        "benchmarks/cnvnator/tree/{aliquot_id}.txt"
    message:
        "Extracting read mappings from BAM file\n"
        "Sample: {wildcards.aliquot_id}"
    shell:
        "cnvnator -root {output} \
            -tree {input} \
            -unique \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Generating a histogram
## See: https://github.com/abyzovlab/CNVnator
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule cnvnator_his:
    input:
        "results/cnvnator/tree/{aliquot_id}.root"
    output:
        ??
    params:
        binsize = 1000,
        mem = CLUSTER_META["cnvnator_his"]["mem"]
    threads:
        CLUSTER_META["cnvnator_his"]["ppn"]
    log:
        "logs/cnvnator/his/{aliquot_id}.log"
    benchmark:
        "benchmarks/cnvnator/his/{aliquot_id}.txt"
    message:
        "Generating a histogram\n"
        "Sample: {wildcards.aliquot_id}"
    shell:
        "cnvnator -root {input} \
            -his {params.binsize} \
            -unique \
            > {log} 2>&1"