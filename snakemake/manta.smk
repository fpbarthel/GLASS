## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Manta
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule manta_config:
    input:
        tumor = lambda wildcards: "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode)),
        normal = lambda wildcards: "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam".format(aliquot_barcode = manifest.getNormal(wildcards.pair_barcode))
    output:
        script = "results/manta/{pair_barcode}/runWorkflow.py"
    params:
        rundir = "results/manta/{pair_barcode}",
        mem = CLUSTER_META["manta_config"]["mem"]
    threads:
        CLUSTER_META["manta_config"]["ppn"]
    conda:
        "../envs/manta.yaml"
    log:
        "logs/manta/config/{pair_barcode}.log"
    benchmark:
        "benchmarks/manta/config/{pair_barcode}.txt"
    message:
        "Configuring Manta for tumor/normal pair\n"
        "Pair: {wildcards.pair_barcode}"
    shell:
        "configManta.py \
            --normalBam {input.normal} \
            --tumorBam {input.tumor} \
            --callRegions {config[svinclude_manta]} \
            --referenceFasta {config[reference_fasta]} \
            --runDir {params.rundir} \
            > {log} 2>&1; "

rule manta_execute:
    input:
        script = "results/manta/{pair_barcode}/runWorkflow.py",
        tumor = lambda wildcards: "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode)),
        normal = lambda wildcards: "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam".format(aliquot_barcode = manifest.getNormal(wildcards.pair_barcode))
    output:
        "results/manta/{pair_barcode}/results/variants/diploidSV.vcf.gz",
        "results/manta/{pair_barcode}/results/variants/somaticSV.vcf.gz",
        "results/manta/{pair_barcode}/results/variants/candidateSV.vcf.gz",
        "results/manta/{pair_barcode}/results/variants/candidateSmallIndels.vcf.gz"
    params:
        mem = CLUSTER_META["manta_execute"]["mem"]
    threads:
        CLUSTER_META["manta_execute"]["ppn"]
    conda:
        "../envs/manta.yaml"
    log:
        "logs/manta/execute/{pair_barcode}.log"
    benchmark:
        "benchmarks/manta/execute/{pair_barcode}.txt"
    message:
        "Running Manta for tumor/normal pair\n"
        "Pair: {wildcards.pair_barcode}"
    shell:
        "{input.script} \
            -m local \
            -j {threads} \
            -g {params.mem} \
            > {log} 2>&1; "

## END ##