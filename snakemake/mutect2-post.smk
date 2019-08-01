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
        sortd = "results/mutect2/dropgt/{case_barcode}.filtered.normalized.sorted.vcf.gz",
        sortdind = "results/mutect2/dropgt/{case_barcode}.filtered.normalized.sorted.vcf.gz.tbi",
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
## Single sample Select variants
## 02/14 remove selectvariants and only drop GTs
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule ssselectvariants:
    input:
        ancient("results/mutect2/ssm2filter/{pair_barcode}.filtered.vcf.gz")
    output:
        normalized = temp("results/mutect2/ssdropgt/{pair_barcode}.filtered.normalized.vcf.gz"),
        sortd = temp("results/mutect2/ssdropgt/{pair_barcode}.filtered.normalized.sorted.vcf.gz"),
        sortdind = temp("results/mutect2/ssdropgt/{pair_barcode}.filtered.normalized.sorted.vcf.gz.tbi"),
        vcf = temp("results/mutect2/ssdropgt/{pair_barcode}.filtered.normalized.sorted.dropgt.vcf.gz"),
        tbi = temp("results/mutect2/ssdropgt/{pair_barcode}.filtered.normalized.sorted.dropgt.vcf.gz.tbi")
    params:
        mem = CLUSTER_META["ssselectvariants"]["mem"]
    threads:
        CLUSTER_META["ssselectvariants"]["ppn"]
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/mutect2/ssselectvariants/{pair_barcode}.log"
    benchmark:
        "benchmarks/mutect2/ssselectvariants/{pair_barcode}.txt"
    message:
        "Single Sample Drop genotypes\n"
        "Case: {wildcards.pair_barcode}"
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
## Merge consensus variant set
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule consensusvcf:
    input:
        vcfs = expand("results/mutect2/dropgt/{case_barcode}.filtered.normalized.sorted.dropgt.vcf.gz", case_barcode = manifest.getSelectedCases()),
        tbis = expand("results/mutect2/dropgt/{case_barcode}.filtered.normalized.sorted.dropgt.vcf.gz.tbi", case_barcode = manifest.getSelectedCases())
    output:
        merged = temp("results/mutect2/consensusvcf/consensus.vcf.gz"),
        normalized = temp("results/mutect2/consensusvcf/consensus.normalized.vcf.gz"),
        final = "results/mutect2/consensusvcf/consensus.normalized.sorted.vcf.gz",
        tbi = "results/mutect2/consensusvcf/consensus.normalized.sorted.vcf.gz.tbi",
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
## Annotate consensus variant set (with VEP via vcf2maf)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule annotate_vep:
    input:
        vcf = "results/mutect2/consensusvcf/consensus.normalized.sorted.vcf.gz"
    output:
        vcf_uncompressed = temp("results/mutect2/annoconsensusvcf/consensus.normalized.sorted.uncompressed.vcf"),
        vcf_reformatted = temp("results/mutect2/annoconsensusvcf/consensus.normalized.sorted.vcf"),
        maf = protected("results/mutect2/annoconsensusvcf/consensus.normalized.sorted.vep.maf")
    params:
        mem = CLUSTER_META["annoconsensusvcf"]["mem"]
    threads:
        CLUSTER_META["annoconsensusvcf"]["ppn"]
    conda:
        "../envs/vcf2maf.yaml"
    log:
        "logs/mutect2/annoconsensusvcf/annoconsensusvcf.log"
    benchmark:
        "benchmarks/mutect2/annoconsensusvcf/annoconsensusvcf.txt"
    message:
        "Annotate consensus variants using VEP"
    shell:
        "bcftools view \
            -Ov \
            -o {output.vcf_uncompressed} {input.vcf} \
            > {log} 2>&1;"

        "vcf2vcf.pl \
            --input-vcf {output.vcf_uncompressed} \
            --output-vcf {output.vcf_reformatted} \
            --ref-fasta {config[reference_fasta]} \
            >> {log} 2>&1;"
        
        "vcf2maf.pl \
            --input-vcf {output.vcf_reformatted} \
            --output-maf {output.maf} \
            --vep-path {config[veppath]} \
            --vep-data {config[vepdata]} \
            --vep-forks 8 \
            --ref-fasta {config[reference_fasta]} \
            --filter-vcf {config[gnomad_vcf]} \
            --tumor-id TUMOR \
            --normal-id NORMAL \
            --species homo_sapiens \
            --ncbi-build GRCh37 \
            >> {log} 2>&1"
            
#Manually performing this step using vep_upload R script
#
#rule vep2db:
#    input:
#        maf = "results/mutect2/annoconsensusvcf/consensus.normalized.sorted.vep.maf",
#        vcf = "results/mutect2/consensusvcf/consensus.normalized.sorted.vcf.gz"
#    output:
#        tsv = "results/mutect2/maf2db/consensus.normalized.sorted.vep.tsv"
#    params:
#        mem = CLUSTER_META["vep2db"]["mem"]
#    threads:
#        CLUSTER_META["vep2db"]["ppn"]
#    conda:
#        "../envs/vep2db.yaml"
#    log:
#        "logs/mutect2/vep2db/vep2db.log"
#    benchmark:
#        "benchmarks/mutect2/vep2db/vep2db.txt"
#    message:
#        "Copy variants to remote"
#    script:
#        "../R/snakemake/vep_upload.R"
            
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Upload genotype to database
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 


rule geno2db:
    input:
        vcf = "results/mutect2/dropgt/{case_barcode}.filtered.normalized.sorted.vcf.gz",
        ssvcf = lambda wildcards: expand("results/mutect2/ssdropgt/{pair_barcode}.filtered.normalized.sorted.vcf.gz", pair_barcode = manifest.getPairsByCase(wildcards.case_barcode))
    output:
        geno = "results/mutect2/geno2db/{case_barcode}.geno.tsv",
        info = "results/mutect2/geno2db/{case_barcode}.info.tsv",
        mfcov = "results/mutect2/geno2db/{case_barcode}.mfcov.tsv",
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
