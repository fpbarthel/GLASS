## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Optitype pipeline
## Author: Frederick Varn
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Extract Chr 6 and unmapped reads to speed pipeline
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule samtoolsview:
    input:
        ancient("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam")
    output:
        chr6 = temp("results/optitype/samtoolsview/{aliquot_barcode}/{aliquot_barcode}.chr6.bam"),
        unmapped = temp("results/optitype/samtoolsview/{aliquot_barcode}/{aliquot_barcode}.unmapped.bam")
    conda:
        "../envs/optitype.yaml"
    log:
        "logs/samtoolsview/{aliquot_barcode}.log"
    message:
        "Extracting potential HLA reads \n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
        "samtools view -b -f3 {input} \"6\" > {output.chr6} 2>> {log};"
        "samtools view -b -f13 {input} > {output.unmapped} 2>> {log}"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Merge Chr 6 and unmapped reads into a single .bam file
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule samtoolsmerge:
    input:
        chr6 = "results/optitype/samtoolsview/{aliquot_barcode}/{aliquot_barcode}.chr6.bam",
        unmapped = "results/optitype/samtoolsview/{aliquot_barcode}/{aliquot_barcode}.unmapped.bam"
    output:
        temp("results/optitype/samtoolsview/{aliquot_barcode}/{aliquot_barcode}.hlareads.bam")
    conda:
        "../envs/optitype.yaml"
    log:
        "logs/samtoolsview/{aliquot_barcode}.log"
    message:
        "Merging potential HLA reads into one file \n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
        "(samtools merge -f {output} {input.chr6} {input.unmapped}) 2>{log}"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Convert merged .bam file to fastq
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule samtofastq:
    input:
        "results/optitype/samtoolsview/{aliquot_barcode}/{aliquot_barcode}.hlareads.bam"
    output:
        F = temp("results/optitype/samtofastq/{aliquot_barcode}/{aliquot_barcode}_r1.fastq"),
        F2 = temp("results/optitype/samtofastq/{aliquot_barcode}/{aliquot_barcode}_r2.fastq")
    params:
        mem = 12
    conda:
        "../envs/optitype.yaml"
    log:
        "logs/samtofastq/{aliquot_barcode}.log"
    message:
        "Converting potential HLA reads into fastq \n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
        "(gatk --java-options \"-Xmx{params.mem}g\" SamToFastq \
        -I {input} \
        -F {output.F} \
        -F2 {output.F2}) 2>{log}"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Use Razers3 to align fastq files to HLA reference
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule razers3:
    input:
        F = "results/optitype/samtofastq/{aliquot_barcode}/{aliquot_barcode}_r1.fastq",
        F2 = "results/optitype/samtofastq/{aliquot_barcode}/{aliquot_barcode}_r2.fastq",
    output:	
        out1 = temp("results/optitype/HLA_reads/{aliquot_barcode}.fished_1.bam"),
        out2 = temp("results/optitype/HLA_reads/{aliquot_barcode}.fished_2.bam")
    conda:
        "../envs/optitype.yaml"
    log:
        "logs/razers3/{aliquot_barcode}.log"
    message:
        "Aligning potential HLA reads to HLA reference \n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
        "razers3 -i 95 -m 1 -dr 0 -tc 1 -o \
        {output.out1} \
        {config[hla_ref]} \
        {input.F} 2>> {log};"
        "razers3 -i 95 -m 1 -dr 0 -tc 1 -o \
        {output.out2} \
        {config[hla_ref]} \
        {input.F2} 2>> {log}"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Convert HLA-aligned .bam file to fastq file
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule bam2fq:
    input:
        F = "results/optitype/HLA_reads/{aliquot_barcode}.fished_1.bam",
        F2 = "results/optitype/HLA_reads/{aliquot_barcode}.fished_2.bam"
    output:
        F = "results/optitype/HLA_reads/{aliquot_barcode}.fished_1.fastq",
        F2 = "results/optitype/HLA_reads/{aliquot_barcode}.fished_2.fastq"
    conda:
        "../envs/optitype.yaml"
    log:
        "logs/bam2fq/{aliquot_barcode}.log"
    message:
        "Converting aligned HLA file to fastq \n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
        "samtools bam2fq {input.F} > {output.F} 2>> {log};"
        "samtools bam2fq {input.F2} > {output.F2} 2>> {log}"
	
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Call HLAtype with OptiType
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule optitype:
    input:
        F = "results/optitype/HLA_reads/{aliquot_barcode}.fished_1.fastq",
        F2 = "results/optitype/HLA_reads/{aliquot_barcode}.fished_2.fastq"
    output:
        pdf = "results/optitype/HLA_calls/{aliquot_barcode}/{aliquot_barcode}_coverage_plot.pdf",
        res = "results/optitype/HLA_calls/{aliquot_barcode}/{aliquot_barcode}_result.tsv"
    params:
        output_dir = "results/optitype/HLA_calls/{aliquot_barcode}/",
        prefix = "{aliquot_barcode}"
    conda:
        "../envs/optitype.yaml"
    log:
        "logs/optitype/{aliquot_barcode}.log"
    message:
        "Calling HLA type \n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
        "(OptiTypePipeline.py \
        -i {input.F} \
        {input.F2} \
        --config conf/optitype_config.ini \
        --prefix {params.prefix} \
        --dna \
        -v \
        -o {params.output_dir}) 2>> {log}"
