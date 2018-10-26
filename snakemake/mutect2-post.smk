## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Select variants
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule selectvariants:
    input:
        ancient("results/mutect2/filterorientation/{pair_barcode}.filterorientation.vcf")
    output:
        vcf = "results/mutect2/final/{pair_barcode}.final.vcf",
        geno = "results/mutect2/final/{pair_barcode}.final.dropgeno.vcf.gz",
        tbi = "results/mutect2/final/{pair_barcode}.final.dropgeno.vcf.gz.tbi"
    params:
        mem = CLUSTER_META["selectvariants"]["mem"]
    threads:
        CLUSTER_META["selectvariants"]["ppn"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/mutect2/selectvariants/{pair_barcode}.log"
    benchmark:
        "benchmarks/mutect2/selectvariants/{pair_barcode}.txt"
    message:
        "Select PASS calls\n"
        "Pair: {wildcards.pair_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g SelectVariants \
            -V {input} \
            -O {output.vcf} \
            --exclude-filtered true \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            > {log} 2>&1;"
        "bcftools view -Oz -G -o {output.geno} {output.vcf} >> {log} 2>&1;"
        "bcftools index -t {output.geno} >> {log} 2>&1;"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Merge consensus variant set
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule consensusvcf:
    input:
        ancient(expand("results/mutect2/final/{pair_barcode}.final.dropgeno.vcf.gz", pair_barcode = manifest.getSelectedPairs()))
    output:
        vcf = "results/mutect2/consensusvcf/consensus.vcf.gz",
        tbi = "results/mutect2/consensusvcf/consensus.vcf.gz.tbi",
        targets = "results/mutect2/consensusvcf/consensus.bed"
    params:
        mem = CLUSTER_META["consensusvcf"]["mem"]
    threads:
        CLUSTER_META["consensusvcf"]["ppn"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/mutect2/consensusvcf/consensusvcf.log"
    benchmark:
        "benchmarks/mutect2/consensusvcf/consensusvcf.txt"
    message:
        "Merge consensus variants"
    shell:
        "bcftools merge -Oz -o {output.vcf} -m none {input} > {log} 2>&1;"
        "bcftools index -t {output.vcf} >> {log} 2>&1;"
        "zcat {output.vcf} | \
            awk '{{OFS=\"\t\"; if (!/^#/){{print $1,$2-1,$2,$4\"/\"$5,\"+\"}}}}' \
            > {output.targets} 2>> {log}"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Convert VF to maf
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule vcf2maf:
    input:
        tumor = lambda wildcards: "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode)),
        normal = lambda wildcards: "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam".format(aliquot_barcode = manifest.getNormal(wildcards.pair_barcode)),
        vcf = "results/mutect2/final/{pair_barcode}.final.vcf"
    output:
        "results/mutect2/vcf2maf/{pair_barcode}.final.maf"
    params:
        mem = CLUSTER_META["vcf2maf"]["mem"],
        tumor_sample_tag = lambda wildcards: manifest.getRGSampleTagByAliquot(manifest.getTumor(wildcards.pair_barcode)),
        normal_sample_tag = lambda wildcards: manifest.getRGSampleTagByAliquot(manifest.getNormal(wildcards.pair_barcode))
    threads:
        CLUSTER_META["vcf2maf"]["ppn"]
    conda:
        "../envs/vcf2maf.yaml"
    log:
        "logs/mutect2/vcf2maf/{pair_barcode}.log"
    benchmark:
        "benchmarks/mutect2/vcf2maf/{pair_barcode}.txt"
    message:
        "Running VEP (variant annotation) on filtered Mutect2 calls and converting output to MAF\n"
        "Pair: {wildcards.pair_barcode}"
    shell:
        "vcf2maf.pl \
            --input-vcf {input.vcf} \
            --output-maf {output} \
            --vep-path {config[veppath]} \
            --vep-data {config[vepdata]} \
            --vep-forks 2 \
            --ref-fasta {config[reference_fasta]} \
            --filter-vcf {config[gnomad_vcf]} \
            --tumor-id {params.tumor_sample_tag} \
            --normal-id {params.normal_sample_tag} \
            --species homo_sapiens \
            --ncbi-build GRCh37 \
            > {log} 2>&1" ## {config[vcf2maf]}

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Annotate consensus variant set
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule annoconsensusvcf:
    input:
        "results/mutect2/consensusvcf/consensus.vcf.gz"
    output:
        "results/mutect2/annoconsensusvcf/consensus.tsv"
    params:
        mem = CLUSTER_META["annoconsensusvcf"]["mem"]
    threads:
        CLUSTER_META["annoconsensusvcf"]["ppn"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/mutect2/annoconsensusvcf/annoconsensusvcf.log"
    benchmark:
        "benchmarks/mutect2/annoconsensusvcf/annoconsensusvcf.txt"
    message:
        "Annotate consensus variants using VEP"
    shell:
        "/projects/barthf/opt/ensembl-vep/vep -i {input} -o {output} \
            --dir {config[vepdata]} \
            --cache \
            --assembly GRCh37 \
            --fork 4 \
            --variant_class\
            --symbol \
            --offline \
            --tab \
            > {log} 2>&1;"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Genotype sample using consensus variant list
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule genotypesample:
    input:
        bam = ancient("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam"),
        vcf = "results/mutect2/consensusvcf/consensus.vcf.gz",
        targets = "results/mutect2/consensusvcf/consensus.bed"
    output:
        vcf = "results/mutect2/freebayes/{aliquot_barcode}.gt.vcf"
        trigger = temp("results/mutect2/freebayes/{aliquot_barcode}.done")
    params:
        mem = CLUSTER_META["genotypesample"]["mem"]
    threads:
        CLUSTER_META["genotypesample"]["ppn"]
    conda:
        "../envs/freebayes.yaml"
    log:
        "logs/mutect2/genotypesample/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/mutect2/genotypesample/{aliquot_barcode}.txt"
    message:
        "Genotype consensus calls using freebayes\n"
        "Aliquot: {wildcards.aliquot_barcode}"
    shell:
        "freebayes \
            -f {config[reference_fasta]} \
            -t {input.targets} \
            -l \
            -@ {input.vcf} {input.bam} \
            > {output.vcf} 2> {log};"
        "touch {output.trigger}"
        
## END ##
