
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run pileup
## See: http://dkoboldt.github.io/varscan/somatic-calling.html
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule pileup:
    input:
        "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam"
    output:
        temp("results/mpileup/{aliquot_id}.pileup")
    params:
        mem = CLUSTER_META["pileup"]["mem"]
    threads:
        CLUSTER_META["pileup"]["ppn"]
    log:
        "logs/pileup/{aliquot_id}.log"
    benchmark:
        "benchmarks/pileup/{aliquot_id}.txt"
    message:
        "Running samtools mpileup\n"
        "Aliquot ID: {wildcards.aliquot_id}"
    shell:
        "samtools mpileup \
            -q 1 \
            -f /fastscratch/verhaak-lab/GLASS-WG/human_g1k_v37_decoy.fasta \
            -o {output} \
            {input} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run VarScan2 to call SNVs on a tumor/normal pair
## See: http://dkoboldt.github.io/varscan/somatic-calling.html
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule varscan:
    input:
        tumor = lambda wildcards: "results/mpileup/{aliquot_id}.pileup".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["tumor_aliquot_id"]),
        normal = lambda wildcards: "results/mpileup/{aliquot_id}.pileup".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["normal_aliquot_id"])
    output:
        "results/varscan2/{pair_id}.snp.Somatic.hc",
        "results/varscan2/{pair_id}.snp.Somatic.lc",
        "results/varscan2/{pair_id}.snp.Germline",
        "results/varscan2/{pair_id}.snp.LOH"
    params:
        mem = CLUSTER_META["varscan"]["mem"],
        outputprefix = "results/varscan2/{pair_id}"
    threads:
        CLUSTER_META["varscan"]["ppn"]
    log:
        "logs/varscan/{pair_id}.log"
    benchmark:
        "benchmarks/varscan/{pair_id}.txt"
    message:
        "Calling SNVs (VarScan2)\n"
        "Pair: {wildcards.pair_id}"
    shell:
        "java -Xmx{params.mem}g -jar jar/VarScan.v2.4.3.jar somatic \
            {input.normal} \
            {input.tumor} \
            {params.outputprefix} \
            --min-coverage 8 \
            --min-coverage-normal 6 \
            --min-coverage-tumor 8 \
            --min-var-freq 0.10 \
            --min-freq-for-hom 0.75 \
            --tumor-purity 1.0 \
            --strand-filter 1 \
            --somatic-p-value 0.05 \
            --output-vcf 1 \
            > {log} 2>&1"

## END ##