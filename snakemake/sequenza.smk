## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run Sequenza
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

CHRS = [ str(x) for x in list(range(1,23)) + ['X', 'Y'] ]

rule bam2seqz:
    input:
        tumor = lambda wildcards: ancient("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode))),
        normal = lambda wildcards: ancient("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam".format(aliquot_barcode = manifest.getNormal(wildcards.pair_barcode)))
    output:
        temp(expand("results/sequenza/bam2seqz/{{pair_barcode}}/seqz_{chr}.gz", chr = CHRS))
    params:
        mem = CLUSTER_META["bam2seqz"]["mem"],
        prefix = "results/sequenza/bam2seqz/{pair_barcode}/seqz",
        chrs = lambda _: " ".join([s for s in CHRS])
    threads:
        CLUSTER_META["bam2seqz"]["ppn"]
    log:
        "logs/sequenza/bam2seqz/{pair_barcode}.log"
    conda:
        "../envs/sequenza.yaml"
    benchmark:
        "benchmarks/sequenza/bam2seqz/{pair_barcode}.txt"
    message:
        "Produce seqz files\n"
        "Pair: {wildcards.pair_barcode}\n"
    shell:"""
        sequenza-utils bam2seqz \
            -gc {config[sequenza][gc_ref]} \
            --fasta {config[reference_fasta]} \
            -n {input.normal} \
            -t {input.tumor} \
            --chromosome {params.chrs} \
            --qlimit {config[sequenza][min_mapq]} \
            --parallel {threads} \
            -o {params.prefix}.gz \
            > {log} 2>&1
        """

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Seqz binning
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule seqz_binning:
    input:
        "results/sequenza/bam2seqz/{pair_barcode}/seqz_{chr}.gz"
    output:
        temp("results/sequenza/seqz_binning/{pair_barcode}/seqz_{chr}.small.seqz.gz")
    params:
        mem = CLUSTER_META["seqz_binning"]["mem"]
    threads:
        CLUSTER_META["seqz_binning"]["ppn"]
    log:
        "logs/sequenza/seqz_binning/{pair_barcode}.{chr}.log"
    conda:
        "../envs/sequenza.yaml"
    benchmark:
        "benchmarks/sequenza/seqz_binning/{pair_barcode}.{chr}.txt"
    message:
        "Post-process by binning the original seqz files\n"
        "Pair: {wildcards.pair_barcode}\n"
        "Chr: {wildcards.chr}"
    shell:"""
        sequenza-utils seqz_binning \
            --seqz {input} \
            -w {config[sequenza][bin_size]} \
            -o {output} \
            > {log} 2>&1
        """

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Merge seqz files
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule mergeseqz:
    input:
        expand("results/sequenza/seqz_binning/{{pair_barcode}}/seqz_{chr}.small.seqz.gz", chr = CHRS)
    output:
        "results/sequenza/mergeseqz/{pair_barcode}.small.seqz.gz"
    params:
        mem = CLUSTER_META["mergeseqz"]["mem"],
        output = "results/sequenza/mergeseqz/{pair_barcode}.small.seqz"
    threads:
        CLUSTER_META["mergeseqz"]["ppn"]
    log:
        "logs/sequenza/mergeseqz/{pair_barcode}.log"
    conda:
        "../envs/gatk4.yaml"
    benchmark:
        "benchmarks/sequenza/mergeseqz/{pair_barcode}.txt"
    message:
        "Merge seqz files\n"
        "Pair: {wildcards.pair_barcode}"
    shell:"""
        zcat {input} | gawk '{{if (NR==1 || $1!="chromosome") {{print $0}}}}' 1> {params.output} 2> {log}
        gzip {params.output}
        """

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run Sequenza R part
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule seqzR:
    input:
        "results/sequenza/mergeseqz/{pair_barcode}.small.seqz.gz"
    output:
        "results/sequenza/seqzR/{pair_barcode}/{pair_barcode}_cellularity.ploidy.txt"
    params:
        prefix = "{pair_barcode}",
        outdir = "results/sequenza/seqzR/{pair_barcode}",
        kmin = config["sequenza"]["kmin"],
        break_method = config["sequenza"]["break_method"]
    threads:
        CLUSTER_META["seqzR"]["ppn"]
    log:
        "logs/sequenza/seqzR/{pair_barcode}.log"
    benchmark:
        "benchmarks/sequenza/seqzR/{pair_barcode}.txt"
    message:
        "Run sequenza R"
    script:
        "../R/snakemake/runSeqz.R"

## END ##
