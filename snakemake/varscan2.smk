
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run pileup
## See: http://dkoboldt.github.io/varscan/somatic-calling.html
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Depricated in favor of using a pipe to conserve disk space
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# rule pileup:
#     input:
#         "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam"
#     output:
#         temp("results/mpileup/{aliquot_id}.pileup")
#     params:
#         mem = CLUSTER_META["pileup"]["mem"]
#     threads:
#         CLUSTER_META["pileup"]["ppn"]
#     log:
#         "logs/pileup/{aliquot_id}.log"
#     benchmark:
#         "benchmarks/pileup/{aliquot_id}.txt"
#     message:
#         "Running samtools mpileup\n"
#         "Aliquot ID: {wildcards.aliquot_id}"
#     shell:
#         "samtools mpileup \
#             -q 1 \
#             -f /fastscratch/verhaak-lab/GLASS-WG/human_g1k_v37_decoy.fasta \
#             -o {output} \
#             {input} \
#             > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run VarScan2 to call SNVs on a tumor/normal pair
## See: http://dkoboldt.github.io/varscan/somatic-calling.html
## Using M2 interval_list (1-based) but converted to bed file (0-based)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

VARSCAN_SET = [ "snp.Somatic.hc", "snp.Somatic.lc", "snp.Germline", "snp.LOH" ]

rule varscan:
    input:
        tumor = lambda wildcards: "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["tumor_aliquot_id"]),
        normal = lambda wildcards: "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["normal_aliquot_id"]),
        intervalbed = lambda wildcards: "{dir}/{interval}/scattered.bed".format(dir=config["wgs_scatterdir"], interval=wildcards.interval)
    output:
        temp("results/varscan2/vs2-scatter/{pair_id}.{interval}.snp.vcf"),
        temp("results/varscan2/vs2-scatter/{pair_id}.{interval}.indel.vcf")
    params:
        mem = CLUSTER_META["varscan"]["mem"],
        outputprefix = "results/varscan2/{pair_id}.{interval}"
    threads:
        CLUSTER_META["varscan"]["ppn"]
    log:
        "logs/varscan/{pair_id}.{interval}.log"
    benchmark:
        "benchmarks/varscan/{pair_id}.{interval}.txt"
    message:
        "Calling SNVs (VarScan2)\n"
        "Pair: {wildcards.pair_id}\n"
        "Interval: {wildcards.interval}"
    shell:
        "TUMOR_MPILEUP=$(printf 'samtools mpileup -q 1 -f {config[reference_fasta]} -l {input.intervalbed} {input.tumor}');"
        "NORMAL_MPILEUP=$(printf 'samtools mpileup -q 1 -f {config[reference_fasta]} -l {input.intervalbed} {input.normal}');"
        "java -Xmx{params.mem}g -Djava.io.tmpdir={config[tempdir]} \
            -jar jar/VarScan.v2.4.3.jar somatic \
            <($NORMAL_MPILEUP) \
            <($TUMOR_MPILEUP) \
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

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Merge Varscan
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Copied and edited from M2-merge SNV rule
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule mergevarscan:
    input:
        snp = lambda wildcards: expand("results/varscan2/vs2-scatter/{pair_id}.{interval}.snp.vcf", pair_id=wildcards.pair_id, interval=WGS_SCATTERLIST),
        indel = lambda wildcards: expand("results/varscan2/vs2-scatter/{pair_id}.{interval}.indel.vcf", pair_id=wildcards.pair_id, interval=WGS_SCATTERLIST)
    output:
        snp = protected("results/varscan2/vcf/{pair_id}.snp.vcf"),
        indel = protected("results/varscan2/vcf/{pair_id}.indel.vcf")
    params:
        mem = CLUSTER_META["mergevarscan"]["mem"]
    threads:
        CLUSTER_META["mergevarscan"]["ppn"]
    log:
        "logs/mergevarscan/{pair_id}.log"
    benchmark:
        "benchmarks/mergevarscan/{pair_id}.txt"
    message:
        "Merging VCF files (M2)\n"
        "Pair: {wildcards.pair_id}"
    run:
        input_snps = " ".join(["-I " + s for s in input['snp']])
        input_indels = " ".join(["-I " + s for s in input['indel']])
        shell("gatk --java-options -Xmx{params.mem}g MergeVcfs \
            {input_snps} \
            -O {output.snp} \
            > {log} 2>&1")
        shell("gatk --java-options -Xmx{params.mem}g MergeVcfs \
            {input_indels} \
            -O {output.indel} \
            > {log} 2>&1")

## END ##