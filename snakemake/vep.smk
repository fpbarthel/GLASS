## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## USE VEP
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule vep:
    input:
        tumor = lambda wildcards: "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["tumor_aliquot_id"]),
        normal = lambda wildcards: "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["normal_aliquot_id"]),
        vcf = "results/m2filter/{pair_id}.filtered2.vcf"
    output:
        "results/vep/{pair_id}.filtered2.anno.maf"
    params:
        mem = CLUSTER_META["vep"]["mem"]
    threads:
        CLUSTER_META["vep"]["ppn"]
    log:
        "logs/vep/{pair_id}.log"
    benchmark:
        "benchmarks/vep/{pair_id}.txt"
    message:
        "Running VEP (variant annotation) on filtered Mutect2 calls\n"
        "Pair: {wildcards.pair_id}"
    shell:
        "TEST_NAM=`samtools view -H {input.tumor} | grep '^@RG' | sed \"s/.*SM:\\([^\\t]*\\).*/\\1/g\" | uniq`;"
        "CTRL_NAM=`samtools view -H {input.normal} | grep '^@RG' | sed \"s/.*SM:\\([^\\t]*\\).*/\\1/g\" | uniq`;"
        "{config[vcf2maf]} \
            --input-vcf {input.vcf} \
            --output-maf {output} \
            --vep-path {config[veppath]} \
            --vep-data {config[vepdata]} \
            --vep-forks 2 \
            --ref-fasta {config[reference_fasta]} \
            --filter-vcf {config[gnomad_vcf]} \
            --tumor-id $TEST_NAM \
            --normal-id $CTRL_NAM \
            --species homo_sapiens \
            --ncbi-build GRCh37 \
            > {log} 2>&1"