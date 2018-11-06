## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Select variants
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# rule selectvariants:
#     input:
#         ancient("results/mutect2/filterorientation/{pair_barcode}.filterorientation.vcf")
#     output:
#         vcf = "results/mutect2/final/{pair_barcode}.final.vcf",
#         geno = "results/mutect2/final/{pair_barcode}.final.dropgeno.vcf.gz",
#         tbi = "results/mutect2/final/{pair_barcode}.final.dropgeno.vcf.gz.tbi"
#     params:
#         mem = CLUSTER_META["selectvariants"]["mem"]
#     threads:
#         CLUSTER_META["selectvariants"]["ppn"]
#     conda:
#         "../envs/gatk4.yaml"
#     log:
#         "logs/mutect2/selectvariants/{pair_barcode}.log"
#     benchmark:
#         "benchmarks/mutect2/selectvariants/{pair_barcode}.txt"
#     message:
#         "Select PASS calls\n"
#         "Pair: {wildcards.pair_barcode}"
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
## Convert VF to maf
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# rule vcf2maf:
#     input:
#         tumor = lambda wildcards: "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode)),
#         normal = lambda wildcards: "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam".format(aliquot_barcode = manifest.getNormal(wildcards.pair_barcode)),
#         vcf = "results/mutect2/final/{pair_barcode}.final.vcf"
#     output:
#         "results/mutect2/vcf2maf/{pair_barcode}.final.maf"
#     params:
#         mem = CLUSTER_META["vcf2maf"]["mem"],
#         tumor_sample_tag = lambda wildcards: manifest.getRGSampleTagByAliquot(manifest.getTumor(wildcards.pair_barcode)),
#         normal_sample_tag = lambda wildcards: manifest.getRGSampleTagByAliquot(manifest.getNormal(wildcards.pair_barcode))
#     threads:
#         CLUSTER_META["vcf2maf"]["ppn"]
#     conda:
#         "../envs/vcf2maf.yaml"
#     log:
#         "logs/mutect2/vcf2maf/{pair_barcode}.log"
#     benchmark:
#         "benchmarks/mutect2/vcf2maf/{pair_barcode}.txt"
#     message:
#         "Running VEP (variant annotation) on filtered Mutect2 calls and converting output to MAF\n"
#         "Pair: {wildcards.pair_barcode}"
#     shell:
#         "vcf2maf.pl \
#             --input-vcf {input.vcf} \
#             --output-maf {output} \
#             --vep-path {config[veppath]} \
#             --vep-data {config[vepdata]} \
#             --vep-forks 2 \
#             --ref-fasta {config[reference_fasta]} \
#             --filter-vcf {config[gnomad_vcf]} \
#             --tumor-id {params.tumor_sample_tag} \
#             --normal-id {params.normal_sample_tag} \
#             --species homo_sapiens \
#             --ncbi-build GRCh37 \
#             > {log} 2>&1" ## {config[vcf2maf]}

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## M2-post
## Use bcftools to post-process mutect final VCF
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule mutect2postprocess:
    input:
        vcf = "results/mutect2/final/{pair_barcode}.final.vcf"
    output:
        normalized = temp("results/mutect2/m2post/{pair_barcode}.normalized.vcf.gz"),
        final = protected("results/mutect2/m2post/{pair_barcode}.normalized.sorted.vcf.gz"),
        tbi = protected("results/mutect2/m2post/{pair_barcode}.normalized.sorted.vcf.gz.tbi")
    params:
        mem = CLUSTER_META["mutect2postprocess"]["mem"]
    threads:
        CLUSTER_META["mutect2postprocess"]["ppn"]
    conda:
        "../envs/freebayes.yaml"
    log:
        "logs/mutect2/mutect2postprocess/{pair_barcode}.log"
    benchmark:
        "benchmarks/mutect2/mutect2postprocess/{pair_barcode}.txt"
    message:
        "Post-process Mutect2 calls using bcftools\n"
        "Pair barcode: {wildcards.pair_barcode}"
    shell:
        "bcftools norm \
            -f {config[reference_fasta]} \
            --check-ref ws \
            -m-both \
            {input.vcf} | \
         vt decompose_blocksub - | \
         bcftools norm -d none | \
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
## Merge consensus variant set
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# rule testme:
#     conda:
#         "../envs/freebayes.yaml"
#     shell:
#         "echo done"

rule consensusvcf:
    input:
        vcfs = expand("results/mutect2/final/{pair_barcode}.final.dropgeno.vcf.gz", pair_barcode = manifest.getSelectedPairs()),
        tbis = expand("results/mutect2/final/{pair_barcode}.final.dropgeno.vcf.gz.tbi", pair_barcode = manifest.getSelectedPairs())
    output:
        merged = "results/mutect2/consensusvcf/consensus.vcf.gz",
        normalized = "results/mutect2/consensusvcf/consensus.normalized.vcf.gz",
        final = "results/mutect2/consensusvcf/consensus.normalized.sorted.vcf.gz",
        tbi = "results/mutect2/consensusvcf/consensus.normalized.sorted.vcf.gz.tbi",
        bed = "results/mutect2/consensusvcf/consensus.normalized.sorted.bed"
    params:
        mem = CLUSTER_META["consensusvcf"]["mem"]
    threads:
        CLUSTER_META["consensusvcf"]["ppn"]
    conda:
        "../envs/freebayes.yaml"
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
            --check-ref ws \
            -m-both \
            {output.merged} | \
         vt decompose_blocksub - | \
         bcftools norm -d none | \
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

        "zcat {output.final} | awk '{{OFS=\"\t\"; \
            if (!/^#/ && (length($4) > 1 || length($5) > 1))\
            {{ print $1,$2-sqrt((length($4)-length($5))^2)-1,$2+sqrt((length($4)-length($5))^2)+1,$4\"/\"$5,\"+\" }} \
            else if (!/^#/) \
            {{ print $1,$2-1,$2,$4\"/\"$5,\"+\" }} \
            }}' > {output.bed} 2>> {log};"


        #     "bcftools merge \
        #     -Oz -o {output.vcf} \
        #     -m none {input} \
        #     > {log} 2>&1;"
        # "bcftools index \
        #     -t {output.vcf} \
        #     >> {log} 2>&1;"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Annotate consensus variant set
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule annoconsensusvcf:
    input:
        vcf = "results/mutect2/consensusvcf/consensus.normalized.sorted.vcf.gz"
    output:
        vcf_uncompressed = temp("results/mutect2/annoconsensusvcf/consensus.normalized.sorted.uncompressed.vcf"),
        vcf_reformatted = temp("results/mutect2/annoconsensusvcf/consensus.normalized.sorted.vcf"),
        maf = "results/mutect2/annoconsensusvcf/consensus.normalized.sorted.maf"
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

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Annotate consensus variant set
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule maf2db:
    input:
        maf = "results/mutect2/annoconsensusvcf/consensus.normalized.sorted.maf",
        vcf = "results/mutect2/consensusvcf/consensus.normalized.sorted.vcf.gz"
    output:
        tsv = "results/mutect2/maf2db/consensus.normalized.sorted.tsv"
    params:
        mem = CLUSTER_META["maf2db"]["mem"]
    threads:
        CLUSTER_META["maf2db"]["ppn"]
    conda:
        "../envs/r.yaml"
    log:
        "logs/mutect2/maf2db/maf2db.log"
    benchmark:
        "benchmarks/mutect2/maf2db/maf2db.txt"
    message:
        "Copy variants to remote"
    script:
        "maf2db.R"


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Genotype sample using consensus variant list
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule genotypesample:
    input:
        bam = ancient("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam"),
        vcf = "results/mutect2/consensusvcf/consensus.normalized.sorted.vcf.gz"
    output:
        vcf = temp("results/mutect2/genotypes/{aliquot_barcode}.consensus.uncompressed.vcf"),
        vcfgz = temp("results/mutect2/genotypes/{aliquot_barcode}.vcf.gz"),
        normalized = temp("results/mutect2/genotypes/{aliquot_barcode}.normalized.vcf.gz"),
        final = protected("results/mutect2/genotypes/{aliquot_barcode}.normalized.sorted.vcf.gz"),
        tbi = protected("results/mutect2/genotypes/{aliquot_barcode}.normalized.sorted.vcf.gz.tbi")
    params:
        mem = CLUSTER_META["genotypesample"]["mem"],
        vcf = "results/mutect2/genotypes/{aliquot_barcode}.vcf",
        readgroup_sample_tag = lambda wildcards: manifest.getRGSampleTagByAliquot(wildcards.aliquot_barcode)
    threads:
        CLUSTER_META["genotypesample"]["ppn"]
    conda:
        "../envs/vcf2maf.yaml"
    log:
        "logs/mutect2/genotypesample/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/mutect2/genotypesample/{aliquot_barcode}.txt"
    message:
        "Genotype consensus calls\n"
        "Aliquot: {wildcards.aliquot_barcode}"
    shell:
        "bcftools view \
            -Ov \
            -o {output.vcf} {input.vcf} \
            > {log} 2>&1;"

        "vcf2vcf.pl \
            --input-vcf {output.vcf} \
            --output-vcf {params.vcf} \
            --tumor-bam {input.bam} \
            --vcf-tumor-id {params.readgroup_sample_tag} \
            --ref-fasta {config[reference_fasta]} \
            >> {log} 2>&1;"

        "bgzip {params.vcf} 2>> {log};"

        "bcftools norm \
            -f {config[reference_fasta]} \
            --check-ref ws \
            -m-both \
            {output.vcfgz} | \
         vt decompose_blocksub - | \
         bcftools norm -d none | \
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
## Genotype sample using freebayes
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule freebayes:
    input:
        bam = ancient("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam"),
        vcf = "results/mutect2/consensusvcf/consensus.normalized.sorted.vcf.gz",
        targets = "results/mutect2/consensusvcf/consensus.normalized.sorted.bed"
    output:
        vcfgz = protected("results/mutect2/freebayes/{aliquot_barcode}.vcf.gz"),
        normalized = temp("results/mutect2/freebayes/{aliquot_barcode}.normalized.vcf.gz"),
        final = protected("results/mutect2/freebayes/{aliquot_barcode}.normalized.sorted.vcf.gz"),
        tbi = protected("results/mutect2/freebayes/{aliquot_barcode}.normalized.sorted.vcf.gz.tbi")
    params:
        vcf = "results/mutect2/freebayes/{aliquot_barcode}.vcf",
        mem = CLUSTER_META["freebayes"]["mem"]
    threads:
        CLUSTER_META["freebayes"]["ppn"]
    conda:
        "../envs/freebayes.yaml"
    log:
        "logs/mutect2/freebayes/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/mutect2/freebayes/{aliquot_barcode}.txt"
    message:
        "Genotype consensus calls using freebayes\n"
        "Aliquot: {wildcards.aliquot_barcode}"
    shell:
        "freebayes \
            -f {config[reference_fasta]} \
            -t {input.targets} \
            -l \
            -@ {input.vcf} \
            {input.bam} \
            > {params.vcf} 2> {log};"

        "bgzip {params.vcf} 2>> {log};"

        "bcftools norm \
            -f {config[reference_fasta]} \
            --check-ref ws \
            -m-both \
            {output.vcfgz} | \
         vt decompose_blocksub - | \
         bcftools norm -d none | \
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
## Genotype sample using freebayes
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# rule freebayes_batch:
#     input:
#         bam = ancient("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam"),
#         vcf = "results/mutect2/consensusvcf/consensus.normalized.sorted.{batch}.vcf.gz",
#         targets = "results/mutect2/consensusvcf/consensus.normalized.sorted.{batch}.bed"
#     output:
#         vcfgz = protected("results/mutect2/freebayes/batch{batch}/{aliquot_barcode}.vcf.gz"),
#         normalized = temp("results/mutect2/freebayes/batch{batch}/{aliquot_barcode}.normalized.vcf.gz"),
#         final = protected("results/mutect2/freebayes/batch{batch}/{aliquot_barcode}.normalized.sorted.vcf.gz"),
#         tbi = protected("results/mutect2/freebayes/batch{batch}/{aliquot_barcode}.normalized.sorted.vcf.gz.tbi")
#     params:
#         vcf = "results/mutect2/freebayes/{aliquot_barcode}.vcf",
#         mem = CLUSTER_META["freebayes"]["mem"]
#     threads:
#         CLUSTER_META["freebayes"]["ppn"]
#     conda:
#         "../envs/freebayes.yaml"
#     log:
#         "logs/mutect2/freebayes/{aliquot_barcode}.log"
#     benchmark:
#         "benchmarks/mutect2/freebayes/{aliquot_barcode}.txt"
#     message:
#         "Genotype consensus calls using freebayes\n"
#         "Aliquot: {wildcards.aliquot_barcode}"
#     shell:
#         "freebayes \
#             -f {config[reference_fasta]} \
#             -t {input.targets} \
#             -l \
#             -@ {input.vcf} \
#             {input.bam} \
#             > {params.vcf} 2> {log};"

#         "bgzip {params.vcf} 2>> {log};"

#         "bcftools norm \
#             -f {config[reference_fasta]} \
#             --check-ref ws \
#             -m-both \
#             {output.vcfgz} | \
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
## Upload genotype to database
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule geno2db:
    input:
        freebayes = "results/mutect2/freebayes/{aliquot_barcode}.normalized.sorted.vcf.gz",
        consensus = "results/mutect2/consensusvcf/consensus.normalized.sorted.vcf.gz"
    output:
        tsv = "results/mutect2/geno2db/{aliquot_barcode}.normalized.sorted.tsv"
    params:
        mutect2 = lambda wildcards: "results/mutect2/m2post/{}.normalized.sorted.vcf.gz".format(manifest.getFirstPair(wildcards.aliquot_barcode)) if manifest.getFirstPair(wildcards.aliquot_barcode) is not None else "",
        mem = CLUSTER_META["geno2db"]["mem"]
    threads:
        CLUSTER_META["geno2db"]["ppn"]
    conda:
        "../envs/r.yaml"
    log:
        "logs/mutect2/geno2db/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/mutect2/geno2db/{aliquot_barcode}.txt"
    message:
        "Add freebayes calls to the database\n"
        "Aliquot: {wildcards.aliquot_barcode}"
    script:
        "geno2db.R"

# ## END ##
