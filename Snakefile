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
    threads:
    	CLUSTER_META["download"]["ppn"]
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
 	threads:
 		CLUSTER_META["revertsam"]["ppn"]
 	shell:
 		"gatk --java-options {config[java_opt]} RevertSam \
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
 	threads:
 		CLUSTER_META["markadapters"]["ppn"]
 	shell:
 		"gatk --java-options {config[java_opt]} MarkIlluminaAdapters \
 			--INPUT={input} \
 			--OUTPUT={output.bam} \
 			--METRICS={output.metric} \
 			--TMP_DIR={config[tempdir]} \
 			2> {log}"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## (1) BAM to FASTQ
## Converts cleaned-up uBAM to FASTQ format, one FASTQ pair per readgroup
##
## (2) Align reads using BWA-MEM
## This is optimzed for >70 bp PE reads and this pipeline will need update if we intend
## to use it with shorter reads
##
## (3) Merge BAM Alignment
## Restore altered data and apply and adjust meta information lost during alignment
##
## See: https://gatkforums.broadinstitute.org/gatk/discussion/6483/how-to-map-and-clean-up-short-read-sequence-data-efficiently
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule samtofastq_bwa_mergebamalignment:
	input:
		"markadapters/{sample_id}/{sample_id}.{rg}.revertsam.markadapters.bam"
	output:
		"bwa/{sample_id}/{sample_id}.{rg}.bam"
	log: 
 		"logs/samtofastq_bwa_mergebamalignment/{sample_id}.{rg}.log"
 	threads:
 		CLUSTER_META["samtofastq_bwa_mergebamalignment"]["ppn"]
	shell:
		"gatk --java-options {config[samtofastq_java_opt]} SamToFastq \
			INPUT={input} \
			FASTQ=/dev/stdout \
			CLIPPING_ATTRIBUTE=XT \
			CLIPPING_ACTION=2 \
			INTERLEAVE=true \
			NON_PF=true \
			VALIDATION_STRINGENCY=SILENT \
			TMP_DIR={config[tempdir]} | \
		bwa mem -M -t {threads} -p {config[fasta]} /dev/stdin | \
		gatk --java-options {config[mergebamalignment_java_opt]} MergeBamAlignment \
			--REFERENCE_SEQUENCE={config[fasta]} \
			--UNMAPPED_BAM={input} \
			--ALIGNED_BAM=/dev/stdin \
			--OUTPUT={output} \
			--CREATE_INDEX=true \
			--ADD_MATE_CIGAR=true \
			--CLIP_ADAPTERS=false \
			--CLIP_OVERLAPPING_READS=true \
			--INCLUDE_SECONDARY_ALIGNMENTS=true \
			--MAX_INSERTIONS_OR_DELETIONS=-1 \
			--PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
			--ATTRIBUTES_TO_RETAIN=XS \
			--TMP_DIR={config[tempdir]} \
			2> {log}"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Mark Duplicates & merge readgroups
## This step marks duplicate reads
## See: https://gatkforums.broadinstitute.org/gatk/discussion/2799
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule markduplicates:
	input:
		expand("bwa/{sample_id}/{sample_id}.{rg}.bam", rg=READGROUPS[sample_id])
	output:
		bam = "markduplicates/{sample_id}.dedup.bam",
		metrics = "markduplicates/{sample_id}.metrics.txt"
	log:
		"logs/markduplicates/{sample_id}.log"
	threads:
		CLUSTER_META["markduplicates"]["ppn"]
	shell:
		"gatk --java-options {config[samtofastq_java_opt]} MarkDuplicates \
			--INPUT={input} \
			--OUTPUT={output.bam} \
			--METRICS_FILE={output.metrics} \
			--CREATE_INDEX=true \
			2> {log}"

## END ##
