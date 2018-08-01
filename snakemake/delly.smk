## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Delly
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule delly_call:
    input:
        tumor = lambda wildcards: "results/align/bqsr/{aliquot_id}.realn.mdup.bqsr.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["tumor_aliquot_id"]),
        normal = lambda wildcards: "results/align/bqsr/{aliquot_id}.realn.mdup.bqsr.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["normal_aliquot_id"])
    output:
        bcf = "results/delly/call/{pair_id}.bcf",
        vcf = "results/delly/call/{pair_id}.vcf.gz"
    params:
        vcftmp = "results/delly/call/{pair_id}.vcf",
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
        "delly call \
            -x {config[svmask_delly]} \
            -o {output.bcf} \
            -g {config[reference_fasta]} \
            {input.tumor} \
            {input.normal} \
            > {log} 2>&1; "
        "bcftools view {output.bcf} > {params.vcftmp} && \
            bgzip -i {params.vcftmp} && \
            bcftools index -t {output.vcf}"