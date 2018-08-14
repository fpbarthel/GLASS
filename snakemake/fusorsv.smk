## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## FusorSV
## Preparing for FusorSV by collecting VCF files from various callers
## See: https://github.com/timothyjamesbecker/FusorSV
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule fusorsv_prepare:
    input:
        delly = "results/delly/filter/{pair_id}.prefilt.bcf",
        lumpy = "results/lumpy/filter/{pair_id}.dict.svtyper.filtered.vcf",
        manta = "results/manta/{pair_id}/results/variants/somaticSV.vcf.gz"
    output:
        delly = "results/fusorsv/prepare/{pair_id}/{pair_id}.delly.vcf",
        lumpy = "results/fusorsv/prepare/{pair_id}/{pair_id}.lumpy.vcf",
        manta = "results/fusorsv/prepare/{pair_id}/{pair_id}.manta.vcf"
    params:
        mem = CLUSTER_META["fusorsv_prepare"]["mem"]
    threads:
        CLUSTER_META["fusorsv_prepare"]["ppn"]
    conda:
        "../envs/fusorsv.yaml"
    log:
        "logs/fusorsv/prepare/{pair_id}.log"
    benchmark:
        "benchmarks/fusorsv/prepare/{pair_id}.txt"
    message:
        "Preparing for FusorSV by collecting VCF files from various callers\n"
        "Pair: {wildcards.pair_id}"
    shell:
        "(bcftools view {input.delly} > {output.delly} && \
            bcftools view {input.lumpy} > {output.lumpy} && \
            bcftools view {input.manta} > {output.manta}) \
            > {log} 2>&1"

#        "bcftools view {output.bcf} > {params.vcftmp} && \
#            bgzip -i {params.vcftmp} && \
#            bcftools index -t {output.vcf}"

## END ##