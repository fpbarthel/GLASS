## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Snakefile for GLASS-WG pipeline
## Floris Barthel 2018
## Development branch
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

configfile: "config.yaml"

## Although this statement goes against all coding conventions, we want it here because we want to run
## everything on a temporary storage while we keep this script safe on a permanent drive
workdir: config["workdir"]

## GDC token file for authentication
KEYFILE 	= config["gdc_token"]

## Metadata
SAMPLES_META 	= json.load(open(config["sample_json"]))
CLUSTER_META 	= json.load(open(config["cluster_json"]))

## Variables
FILENAMES 	= [item['file_name'] for item in SAMPLES_META]
UUIDS 		= [item['id'] for item in SAMPLES_META]
SAMPLES 	= [item['sample_id'] for item in SAMPLES_META]

## Targets
DOWNLOAD_TARGETS 	= ["download/{}/{}".format(uuid, filename) for uuid, filename in zip(UUIDS, FILENAMES)] 
REVERTSAM_TARGETS 	= ["revertsam/{}.revertsam.bam".format(sample_id) for sample_id in zip(SAMPLES)] 

## Set targets to final set of targets (as to skip intermediate steps if already complete)
TARGETS = REVERTSAM_TARGETS

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Master rule
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
rule all:
	input: TARGETS
  
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Download BAM file from GDC
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule download:
    output:
        "download/{uuid}/{filename}.bam"
    threads: CLUSTER_META["download"]["ppn"]
    message:
    	"Downloading UUID {wildcards.uuid} (file {wildcards.filename}) from GDC"
    log:
    	"logs/download/{uuid}.log"
    shell:
    	"gdc-client download -d download -n {threads} -t ~/gdc_token.key {wildcards.uuid} 2> {log}"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Revert GDC-aligned legacy BAM to unaligned SAM file
## Clear BAM attributes, such as re-calibrated base quality scores, alignment information
## Moreover, BAM file is split into per-readgroup BAM files
## See: https://gatkforums.broadinstitute.org/gatk/discussion/6484
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule revertsam:
 	input:
 		"download/{uuid}/{filename}.bam"
 	params:
 		dir = "revertsam/{sample_id}"
 	log: 
 		"logs/revertsam/{sample_id}.log"
 	threads: CLUSTER_META["revertsam"]["ppn"]
 	shell:
 		"java {config[java_opt]} \
 			-jar {config[gatk_jar]} RevertSam \
 			--INPUT={input} \
 			--OUTPUT=${params.dir} \
 			--OUTPUT_BY_READGROUP=true \
 			--RESTORE_ORIGINAL_QUALITIES=true \
 			--VALIDATION_STRINGENCY=SILENT \
 			--ATTRIBUTE_TO_CLEAR=AS \
 			--ATTRIBUTE_TO_CLEAR=FT \
 			--ATTRIBUTE_TO_CLEAR=CO \
 			--ATTRIBUTE_TO_CLEAR=XT \
 			--ATTRIBUTE_TO_CLEAR=XN \
 			--ATTRIBUTE_TO_CLEAR=OC \
 			--ATTRIBUTE_TO_CLEAR=OP \
 			--SANITIZE=true \
 			--SORT_ORDER=queryname \
 			--TMP_DIR={config[tempdir]} \
 			2> {log}"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Mark Illumina Adapters
## Add XT tag to read records to mark the 5' start position of adapter sequences
## Adapter sequences are then removed by subsequent steps
## See: https://gatkforums.broadinstitute.org/gatk/discussion/6483/how-to-map-and-clean-up-short-read-sequence-data-efficiently
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule markadapters:
	input:
 		"revertsam/{sample_id}/{sample_id}.{rg}.revertsam.bam"
 	output:
 		bam 	= "markadapters/{sample_id}/{sample_id}.{rg}.revertsam.markadapters.bam",
 		metric 	= "markadapters/{sample_id}/{sample_id}.{rg}.markadapters.metrics.txt"
 	log: 
 		"logs/markadapters/{sample_id}.{rg}.log"
 	threads: CLUSTER_META["markadapters"]["ppn"]
 	shell:
 		"java {config[java_opt]} \
 			-jar {config[gatk_jar]} MarkIlluminaAdapters \
 			I={input} \
 			O={output.bam} \
 			M={output.metric} \
 			TMP_DIR={config[tempdir]} \
 			2> {log}"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## BAM to FASTQ
## Converts cleaned-up uBAM to FASTQ format, one FASTQ pair per readgroup
## See: https://gatkforums.broadinstitute.org/gatk/discussion/6483/how-to-map-and-clean-up-short-read-sequence-data-efficiently
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule samtofastq:
	input:
		"markadapters/{sample_id}/{sample_id}.{rg}.revertsam.markadapters.bam"
	output:
		r1 = "fq/{sample_id}/{sample_id}.{rg}_R1.fq.gz",
		r2 = "fq/{sample_id}/{sample_id}.{rg}_R2.fq.gz"
		un = "fq/{sample_id}/{sample_id}.{rg}_unpaired.fq.gz"
	shell:
		"java {config[java_opt]} \
 			-jar {config[gatk_jar]} SamToFastq \
			INPUT={input} \
			VALIDATION_STRINGENCY=SILENT \
			FASTQ={output.r1} \
			SECOND_END_FASTQ={output.r2} \
			UNPAIRED_FASTQ={output.un} \
			TMP_DIR={config[tempdir]} \
			2> {log}"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Align reads using BWA-MEM
## This is optimzed for >70 bp PE reads and this pipeline will need update if we intend
## to use it with shorter reads
## See: https://gatkforums.broadinstitute.org/gatk/discussion/6483/how-to-map-and-clean-up-short-read-sequence-data-efficiently
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule aln:
	input:
		r1 = "fq/{sample_id}/{sample_id}.{rg}_R1.fq.gz",
		r2 = "fq/{sample_id}/{sample_id}.{rg}_R2.fq.gz"
	output:
		"bwa/{sample_id}.bam"
	shell:
		"bwa mem -M -t {threads} -R {config[fasta]} {input.r1} {input.r2} > {output} 2> {log}"

# rule samtools_sort:
#     input:
#         "mapped_reads/{sample}.bam"
#     output:
#         "sorted_reads/{sample}.bam"
#     shell:
#         "samtools sort -T sorted_reads/{wildcards.sample} "
#         "-O bam {input} > {output}"

# rule samtools_index:
#     input:
#         "sorted_reads/{sample}.bam"
#     output:
#         "sorted_reads/{sample}.bam.bai"
#     shell:
#         "samtools index {input}"

# rule bcftools_call:
#     input:
#         fa="data/genome.fa",
#         bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
#         bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
#     output:
#         "calls/all.vcf"
#     log:
#     	"logs/bcftools/all_calls.log"
#     shell:
#         "samtools mpileup -g -f {input.fa} {input.bam} | "
#         "bcftools call -mv - > {output} 2> {log}"

# rule report:
#     input:
#         "calls/all.vcf"
#     output:
#         "report.html"
#     run:
#         from snakemake.utils import report
#         with open(input[0]) as vcf:
#             n_calls = sum(1 for l in vcf if not l.startswith("#"))

#         report("""
#         An example variant calling workflow
#         ===================================

#         Reads were mapped to the Yeast
#         reference genome and variants were called jointly with
#         SAMtools/BCFtools.

#         This resulted in {n_calls} variants (see Table T1_).
#         """, output[0], T1=input[0])