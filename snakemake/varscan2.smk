
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
            -f {config[reference_fasta]} \
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
        "gatk --java-options -Xmx{params.mem}g Mutect2 \
            -R {config[reference_fasta]} \
            -I {input.tumor} \
            -I {input.normal} \
            -L {input.intervallist} \
            --tumor-sample $TEST_NAM \
            --normal-sample $CTRL_NAM \
            --panel-of-normals {input.pon} \
            --germline-resource {config[gnomad_vcf]} \
            --af-of-alleles-not-in-resource {config[af_of_alleles_not_in_resource]} \
            --dont-use-soft-clipped-bases true \
            --standard-min-confidence-threshold-for-calling 20 \
            -O {output.vcf} \
            -bamout {output.bam} \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            > {log} 2>&1"