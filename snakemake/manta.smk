## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Manta
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule manta_config:
    input:
        tumor = lambda wildcards: "results/align/bqsr/{aliquot_id}.realn.mdup.bqsr.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["tumor_aliquot_id"]),
        normal = lambda wildcards: "results/align/bqsr/{aliquot_id}.realn.mdup.bqsr.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["normal_aliquot_id"])
    output:
        script = "results/manta/{pair_id}/runWorkflow.py"
    params:
        rundir = "results/manta/{pair_id}",
        mem = CLUSTER_META["manta_config"]["mem"]
    threads:
        CLUSTER_META["manta_config"]["ppn"]
    conda:
        "envs/manta.yaml"
    log:
        "logs/manta/config/{pair_id}.log"
    benchmark:
        "benchmarks/manta/config/{pair_id}.txt"
    message:
        "Configuring Manta for tumor/normal pair\n"
        "Pair: {wildcards.pair_id}"
    shell:
        "configManta.py \
            --normalBam {input.normal} \
            --tumorBam {input.tumor} \
            --referenceFasta {config[reference_fasta]} \
            --runDir {params.rundir} \
            > {log} 2>&1; "

rule manta_execute:
    input:
        script = "results/manta/{pair_id}/runWorkflow.py",
        tumor = lambda wildcards: "results/align/bqsr/{aliquot_id}.realn.mdup.bqsr.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["tumor_aliquot_id"]),
        normal = lambda wildcards: "results/align/bqsr/{aliquot_id}.realn.mdup.bqsr.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["normal_aliquot_id"])
    output:
        "results/manta/{pair_id}/results/variants/diploidSV.vcf.gz",
        "results/manta/{pair_id}/results/variants/somaticSV.vcf.gz",
        "results/manta/{pair_id}/results/variants/candidateSV.vcf.gz",
        "results/manta/{pair_id}/results/variants/candidateSmallIndels.vcf.gz"
    params:
        mem = CLUSTER_META["manta_execute"]["mem"]
    threads:
        CLUSTER_META["manta_execute"]["ppn"]
    conda:
        "envs/manta.yaml"
    log:
        "logs/manta/execute/{pair_id}.log"
    benchmark:
        "benchmarks/manta/execute/{pair_id}.txt"
    message:
        "Running Manta for tumor/normal pair\n"
        "Pair: {wildcards.pair_id}"
    shell:
        "{input.script} \
            -m local \
            -j {threads} \
            -g {params.mem} \
            > {log} 2>&1; "