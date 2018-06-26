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

## END ##