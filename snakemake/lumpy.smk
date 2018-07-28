## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Extract and sort split reads
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule extractsplitter:
    input:
        "results/align/bqsr/{aliquot_id}.realn.mdup.bqsr.bam"
    output:
        unsorted = temp("results/lumpy/split/{aliquot_id}.realn.mdup.bqsr.splitters.unsorted.bam"),
        sorted = "results/lumpy/split/{aliquot_id}.realn.mdup.bqsr.splitters.sorted.bam"
    params:
        prefix = "results/lumpy/split/{aliquot_id}",
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
        "results/align/bqsr/{aliquot_id}.realn.mdup.bqsr.bam"
    output:
        unsorted = temp("results/lumpy/discordant/{aliquot_id}.realn.mdup.bqsr.discordant.unsorted.bam"),
        sorted = "results/lumpy/discordant/{aliquot_id}.realn.mdup.bqsr.discordant.sorted.bam"
    params:
        prefix = "results/lumpy/discordant/{aliquot_id}",
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
## Added HEXDUMP defintion because it was undefined for some reason
## Added gatk UpdateVCFSequenceDictionary because bcftools index requires it
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule lumpy_call:
    input:
        tumor = lambda wildcards: "results/align/bqsr/{aliquot_id}.realn.mdup.bqsr.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["tumor_aliquot_id"]),
        normal = lambda wildcards: "results/align/bqsr/{aliquot_id}.realn.mdup.bqsr.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["normal_aliquot_id"]),
        discordant_tumor = lambda wildcards: "results/lumpy/discordant/{aliquot_id}.realn.mdup.bqsr.discordant.sorted.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["tumor_aliquot_id"]),
        discordant_normal = lambda wildcards: "results/lumpy/discordant/{aliquot_id}.realn.mdup.bqsr.discordant.sorted.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["normal_aliquot_id"]),
        split_tumor = lambda wildcards: "results/lumpy/split/{aliquot_id}.realn.mdup.bqsr.splitters.sorted.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["tumor_aliquot_id"]),
        split_normal = lambda wildcards: "results/lumpy/split/{aliquot_id}.realn.mdup.bqsr.splitters.sorted.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["normal_aliquot_id"])
    output:
        vcf = temp("results/lumpy/call/{pair_id}.dict.vcf.gz"),
        vcfsorted = protected("results/lumpy/call/{pair_id}.dict.sorted.vcf.gz")
    params:
        vcftmp = "results/lumpy/call/{pair_id}.vcf",
        vcftmpdict = "results/lumpy/call/{pair_id}.dict.vcf",
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
        "export HEXDUMP=`which hexdump || true`; "
        "lumpyexpress \
            -B {input.tumor},{input.normal} \
            -S {input.split_tumor},{input.split_normal} \
            -D {input.discordant_tumor},{input.discordant_normal} \
            -T {config[tempdir]}/{wildcards.pair_id} \
            -o {params.vcftmp} \
            > {log} 2>&1; "
        "gatk --java-options -Xmx{params.mem}g UpdateVCFSequenceDictionary \
            -V {params.vcftmp} \
            --source-dictionary {config[reference_dict]} \
            --replace true \
            -O {params.vcftmpdict} \
            >> {log} 2>&1; "
        "bgzip -i {params.vcftmpdict} && \
            bcftools sort -O z -o {output.vcfsorted} {output.vcf} && \
            bcftools index -t {output.vcfsorted} \
            >> {log} 2>&1"

## END ##