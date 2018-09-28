## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Telomere content estimates from BAM file
## See: https://github.com/abyzovlab/CNVnator
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule telseq_run:
    input:
        "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam"
    output:
        protected("results/telseq/{aliquot_barcode}.telseq.txt")
    params:
        mem = CLUSTER_META["telseq_run"]["mem"]
    threads:
        CLUSTER_META["telseq_run"]["ppn"]
    conda:
        "../envs/telseq.yaml"
    log:
        "logs/telseq/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/telseq/{aliquot_barcode}.txt"
    message:
        "Telomere content estimates from BAM file\n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
        "telseq -o {output} \
            -r {config[telseq_r]} \
            {input} \
            > {log} 2>&1"

## END ##