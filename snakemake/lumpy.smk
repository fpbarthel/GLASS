## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Extract and sort split reads
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule extractsplitter:
    input:
        "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam"
    output:
        unsorted = temp("results/lumpy/{aliquot_id}.realn.mdup.bqsr.splitters.unsorted.bam"),
        sorted = "results/lumpy/{aliquot_id}.realn.mdup.bqsr.splitters.sorted.bam"
    params:
        prefix = "results/lumpy/{aliquot_id}",
        mem = CLUSTER_META["extractsplitter"]["mem"]
    threads:
        CLUSTER_META["extractsplitter"]["ppn"]
    conda:
        "envs/lumpy-sv.yaml"
    log:
        "logs/extractsplitter/{aliquot_id}.log"
    benchmark:
        "benchmarks/extractsplitter/{aliquot_id}.txt"
    message:
        "Extracting and sorting split reads\n"
        "Sample: {wildcards.aliquot_id}"
    shell:
        "samtools view -h {input} | \
            extractSplitReads_BwaMem -i stdin | \
            samtools view -Sb - \
            -o {output.unsorted} \
            > {log} 2>&1;"
        "samtools sort \
            -o {output.sorted} \
            -O bam \
            -T {params.prefix} \
            {output.unsorted} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Extract and sort discordant reads
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule extractdiscordant:
    input:
        "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam"
    output:
        unsorted = temp("results/lumpy/{aliquot_id}.realn.mdup.bqsr.discordant.unsorted.bam"),
        sorted = "results/lumpy/{aliquot_id}.realn.mdup.bqsr.discordant.sorted.bam"
    params:
        prefix = "results/lumpy/{aliquot_id}",
        mem = CLUSTER_META["extractdiscordant"]["mem"]
    threads:
        CLUSTER_META["extractdiscordant"]["ppn"]
    conda:
        "envs/lumpy-sv.yaml"
    log:
        "logs/extractdiscordant/{aliquot_id}.log"
    benchmark:
        "benchmarks/extractdiscordant/{aliquot_id}.txt"
    message:
        "Extracting and sorting discordant reads\n"
        "Sample: {wildcards.aliquot_id}"
    shell:
        "samtools view -b -F 1294 {input} \
            -o {output.unsorted} \
            > {log} 2>&1;"
        "samtools sort \
            -o {output.sorted} \
            -O bam \
            -T {params.prefix} \
            {output.unsorted} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## LUMPY
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule lumpy_call:
    input:
        tumor = lambda wildcards: "results/align/bqsr/{aliquot_id}.realn.mdup.bqsr.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["tumor_aliquot_id"]),
        normal = lambda wildcards: "results/align/bqsr/{aliquot_id}.realn.mdup.bqsr.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["normal_aliquot_id"]),
        discordant_tumor = lambda wildcards: "results/lumpy/{aliquot_id}.realn.mdup.bqsr.discordant.sorted.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["tumor_aliquot_id"]),
        discordant_normal = lambda wildcards: "results/lumpy/{aliquot_id}.realn.mdup.bqsr.discordant.sorted.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["normal_aliquot_id"]),
        split_tumor = lambda wildcards: "results/lumpy/{aliquot_id}.realn.mdup.bqsr.splitters.sorted.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["tumor_aliquot_id"]),
        split_normal = lambda wildcards: "results/lumpy/{aliquot_id}.realn.mdup.bqsr.splitters.sorted.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["normal_aliquot_id"])
    output:
        vcf = "results/lumpy/call/{pair_id}.vcf.gz"
    params:
        vcftmp = "results/lumpy/call/{pair_id}.vcf",
        mem = CLUSTER_META["lumpy_call"]["mem"]
    threads:
        CLUSTER_META["lumpy_call"]["ppn"]
    conda:
        "envs/lumpy-sv.yaml"
    log:
        "logs/lumpy/call/{pair_id}.log"
    benchmark:
        "benchmarks/lumpy/call/{pair_id}.txt"
    message:
        "Calling LUMPY on tumor/normal pair\n"
        "Pair: {wildcards.pair_id}"
    shell:
        "lumpyexpress \
            -B {input.tumor},{input.normal} \
            -S {input.split_tumor},{input.split_normal} \
            -D {input.discordant_tumor},{input.discordant_normal} \
            -o {params.vcftmp} \
            > {log} 2>&1; "
        "bgzip -i {params.vcftmp} && \
            bftools index -t {output.vcf}"

## END ##