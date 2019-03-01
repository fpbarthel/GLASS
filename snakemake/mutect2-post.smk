## Changes 02/11, updates for 2nd data freeze
## - M2 multi-sample calling
## - batch-specific PON
## - coverage stats per gene


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Select variants
## 02/14 remove selectvariants and only drop GTs
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule selectvariants:
    input:
        ancient("results/mutect2/m2filter/{case_barcode}.filtered.vcf.gz")
    output:
        normalized = temp("results/mutect2/dropgt/{case_barcode}.filtered.normalized.vcf.gz"),
        sortd = temp("results/mutect2/dropgt/{case_barcode}.filtered.normalized.sorted.vcf.gz"),
        sortdind = temp("results/mutect2/dropgt/{case_barcode}.filtered.normalized.sorted.vcf.gz.tbi"),
        vcf = temp("results/mutect2/dropgt/{case_barcode}.filtered.normalized.sorted.dropgt.vcf.gz"),
        tbi = temp("results/mutect2/dropgt/{case_barcode}.filtered.normalized.sorted.dropgt.vcf.gz.tbi")
    params:
        mem = CLUSTER_META["selectvariants"]["mem"]
    threads:
        CLUSTER_META["selectvariants"]["ppn"]
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/mutect2/selectvariants/{case_barcode}.log"
    benchmark:
        "benchmarks/mutect2/selectvariants/{case_barcode}.txt"
    message:
        "Drop genotypes\n"
        "Case: {wildcards.case_barcode}"
    shell:
        "bcftools norm \
            -f {config[reference_fasta]} \
            --check-ref we \
            -m-both \
            {input} | \
         bcftools view \
            -Oz \
            -o {output.normalized} \
            2>> {log};"
       
        "bcftools sort \
            -Oz \
            -o {output.sortd} \
            {output.normalized} \
            2>> {log};"
        
        "bcftools index \
            -t {output.sortd} \
            2>> {log};"

        "bcftools view -Oz -G -o {output.vcf} {output.sortd} >> {log} 2>&1;"
        "bcftools index -t {output.vcf} >> {log} 2>&1;"


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Funcotator
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# rule vcf2maf:
#     input:
#         ancient("results/mutect2/m2filter/{case_barcode}.filtered.vcf.gz")
#     output:
#         vcf = "results/mutect2/final/{case_barcode}.final.vcf.gz",
#         geno = "results/mutect2/final/{case_barcode}.final.dropgeno.vcf.gz",
#         tbi = "results/mutect2/final/{case_barcode}.final.dropgeno.vcf.gz.tbi"
#     params:
#         mem = CLUSTER_META["selectvariants"]["mem"]
#     threads:
#         CLUSTER_META["selectvariants"]["ppn"]
#     conda:
#         "../envs/gatk4.yaml"
#     log:
#         "logs/mutect2/selectvariants/{case_barcode}.log"
#     benchmark:
#         "benchmarks/mutect2/selectvariants/{case_barcode}.txt"
#     message:
#         "Select PASS calls\n"
#         "Case: {wildcards.case_barcode}"
#     shell:
#         "gatk --java-options -Xmx{params.mem}g SelectVariants \
#             -V {input} \
#             -O {output.vcf} \
#             --exclude-filtered true \
#             --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
#             > {log} 2>&1;"
#         "bcftools view -Oz -G -o {output.geno} {output.vcf} >> {log} 2>&1;"
#         "bcftools index -t {output.geno} >> {log} 2>&1;"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## M2-post
## Use bcftools to post-process mutect final VCF
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# rule mutect2postprocess:
#     input:
#         vcf = "results/mutect2/final/{case_barcode}.final.vcf"
#     output:
#         normalized = temp("results/mutect2/m2post/{case_barcode}.normalized.vcf.gz"),
#         final = protected("results/mutect2/m2post/{case_barcode}.normalized.sorted.vcf.gz"),
#         tbi = protected("results/mutect2/m2post/{case_barcode}.normalized.sorted.vcf.gz.tbi")
#     params:
#         mem = CLUSTER_META["mutect2postprocess"]["mem"]
#     threads:
#         CLUSTER_META["mutect2postprocess"]["ppn"]
#     conda:
#         "../envs/freebayes.yaml"
#     log:
#         "logs/mutect2/mutect2postprocess/{pair_barcode}.log"
#     benchmark:
#         "benchmarks/mutect2/mutect2postprocess/{pair_barcode}.txt"
#     message:
#         "Post-process Mutect2 calls using bcftools\n"
#         "Pair barcode: {wildcards.pair_barcode}"
#     shell:
#         "bcftools norm \
#             -f {config[reference_fasta]} \
#             --check-ref ws \
#             -m-both \
#             {input.vcf} | \
#          vt decompose_blocksub - | \
#          bcftools norm -d none | \
#          bcftools view \
#             -Oz \
#             -o {output.normalized} \
#             2>> {log};"
       
#         "bcftools sort \
#             -Oz \
#             -o {output.final} \
#             {output.normalized} \
#             2>> {log};"
        
#         "bcftools index \
#             -t {output.final} \
#             2>> {log};"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Merge consensus variant set
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule consensusvcf:
    input:
        vcfs = expand("results/mutect2/dropgt/{case_barcode}.filtered.normalized.sorted.dropgt.vcf.gz", case_barcode = manifest.getSelectedCases()),
        tbis = expand("results/mutect2/dropgt/{case_barcode}.filtered.normalized.sorted.dropgt.vcf.gz.tbi", case_barcode = manifest.getSelectedCases())
    output:
        merged = temp("results/mutect2/consensusvcf/consensus.vcf.gz"),
        normalized = temp("results/mutect2/consensusvcf/consensus.normalized.vcf.gz"),
        final = temp("results/mutect2/consensusvcf/consensus.normalized.sorted.vcf.gz"),
        tbi = temp("results/mutect2/consensusvcf/consensus.normalized.sorted.vcf.gz.tbi"),
        
    params:
        mem = CLUSTER_META["consensusvcf"]["mem"]
    threads:
        CLUSTER_META["consensusvcf"]["ppn"]
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/mutect2/consensusvcf/consensusvcf.log"
    benchmark:
        "benchmarks/mutect2/consensusvcf/consensusvcf.txt"
    message:
        "Merge consensus variants"
    shell:
        "bcftools merge \
            -m none {input.vcfs} \
            -Oz \
            -o {output.merged} \
            2>> {log};"
        
        "bcftools norm \
            -f {config[reference_fasta]} \
            --check-ref we \
            -m-both \
            {output.merged} | \
         bcftools view \
            -Oz \
            -o {output.normalized} \
            2>> {log};"
       
        "bcftools sort \
            -Oz \
            -o {output.final} \
            {output.normalized} \
            2>> {log};"
        
        "bcftools index \
            -t {output.final} \
            2>> {log};"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Annotate consensus variant set
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule annoconsensusvcf:
    input:
        vcf = "results/mutect2/consensusvcf/consensus.normalized.sorted.vcf.gz",
        tbi = "results/mutect2/consensusvcf/consensus.normalized.sorted.vcf.gz.tbi"
    output:
        vcf = temp("results/mutect2/consensusvcf/consensus.normalized.sorted.funcotated.vcf"),
        idx = temp("results/mutect2/consensusvcf/consensus.normalized.sorted.funcotated.vcf.idx"),
        gz = protected("results/mutect2/consensusvcf/consensus.normalized.sorted.funcotated.vcf.gz"),
        tbi = protected("results/mutect2/consensusvcf/consensus.normalized.sorted.funcotated.vcf.gz.tbi")
    params:
        mem = CLUSTER_META["annoconsensusvcf"]["mem"]
    threads:
        CLUSTER_META["annoconsensusvcf"]["ppn"]
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/mutect2/annoconsensusvcf/annoconsensusvcf.log"
    benchmark:
        "benchmarks/mutect2/annoconsensusvcf/annoconsensusvcf.txt"
    message:
        "Annotate consensus variants using Funcotator"
    shell:
        "gatk Funcotator \
            --variant {input.vcf} \
            --reference {config[reference_fasta]} \
            --ref-version hg19 \
            --data-sources-path {config[funcotator_dir]} \
            --output {output.vcf} \
            --output-file-format VCF \
            --remove-filtered-variants false \
            2>> {log};"
        
        "bcftools view \
            -Oz \
            -o {output.gz} \
            {output.vcf} \
            2>> {log};"

        "bcftools index \
            -t {output.gz} \
            2>> {log};"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Annotate consensus variant set
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule maf2db:
    input:
        vcf = "results/mutect2/consensusvcf/consensus.normalized.sorted.funcotated.vcf.gz"
    output:
        tsv = "results/mutect2/consensusvcf/consensus.normalized.sorted.funcotated.tsv"
    params:
        mem = CLUSTER_META["maf2db"]["mem"]
    threads:
        CLUSTER_META["maf2db"]["ppn"]
    #conda:
    #    "../envs/r.yaml"
    log:
        "logs/mutect2/maf2db/maf2db.log"
    benchmark:
        "benchmarks/mutect2/maf2db/maf2db.txt"
    message:
        "Copy variants to TSV for uploading to database"
    script:
        "../R/snakemake/snv2db.R"


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Upload genotype to database
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule geno2db:
    input:
        vcf = "results/mutect2/dropgt/{case_barcode}.filtered.normalized.sorted.vcf.gz"
    output:
        geno = "results/mutect2/geno2db/{case_barcode}.geno.tsv",
        info = "results/mutect2/geno2db/{case_barcode}.info.tsv"
    params:
        mem = CLUSTER_META["geno2db"]["mem"]
    threads:
        CLUSTER_META["geno2db"]["ppn"]
    #conda:
    #    "../envs/r.yaml"
    log:
        "logs/mutect2/geno2db/{case_barcode}.log"
    benchmark:
        "benchmarks/mutect2/geno2db/{case_barcode}.txt"
    message:
        "Copy M2 calls to TSV for uploading to database\n"
        "Aliquot: {wildcards.case_barcode}"
    script:
        "../R/snakemake/geno2db.R"

# ## END ##
