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
KEYFILE     = config["gdc_token"]

## Metadata
SAMPLES_META    = json.load(open(config["sample_json"]))
CLUSTER_META    = json.load(open(config["cluster_json"]))

## Variables
FILENAMES   = [item['file_name'] for item in SAMPLES_META]
UUIDS       = [item['id'] for item in SAMPLES_META]
SAMPLES     = [item['sample_id'] for item in SAMPLES_META]

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Master rule
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule all:
    input: expand("bqsr/{ID}.realn.dedup.bqsr.bam", ID=SAMPLES)
  
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
    benchmark:
        "benchmarks/download/{uuid}.txt"
    shell:
        "gdc-client download \
            -d download \
            -n {threads} \
            -t ~/gdc_token.key \
            {wildcards.uuid} \
            2> {log}"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Function that resolves UUID and FILENAME from SAMPLE_ID
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

def get_gdc_bam_filename(wildcards):
    i = SAMPLES.index(wildcards.sample_id)
    return "download/{}/{}".format(UUIDS[i], FILENAMES[i])

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Revert GDC-aligned legacy BAM to unaligned SAM file
## Clear BAM attributes, such as re-calibrated base quality scores, alignment information
## Moreover, BAM file is split into per-readgroup BAM files
## See: https://gatkforums.broadinstitute.org/gatk/discussion/6484
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule revertsam:
    input:
        get_gdc_bam_filename
    output:
        temp(dynamic("revertsam/{sample_id}/{readgroup}.bam"))
    params:
        dir = "revertsam/{sample_id}"
    log: 
        "logs/revertsam/{sample_id}.log"
    threads:
        CLUSTER_META["revertsam"]["ppn"]
    benchmark:
        "benchmarks/revertsam/{sample_id}.txt"
    message:
        "Reverting {wildcards.sample_id} back to unaligned BAM file, stripping any previous "
        "pre-processing and restoring original base quality scores. Output files are split "
        "by readgroup."
    shell:
        "gatk --java-options {config[standard_java_opt]} RevertSam \
            --INPUT={input} \
            --OUTPUT={params.dir} \
            --OUTPUT_BY_READGROUP=true \
            --OUTPUT_BY_READGROUP_FILE_FORMAT=bam \
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
        "revertsam/{sample_id}/{readgroup}.bam"
    output:
        bam = temp("markadapters/{sample_id}/{sample_id}.{readgroup}.revertsam.markadapters.bam"),
        metric = "markadapters/{sample_id}/{sample_id}.{readgroup}.markadapters.metrics.txt"
    threads:
        CLUSTER_META["markadapters"]["ppn"]
    log: 
        dynamic("logs/markadapters/{sample_id}.{readgroup}.log")
    benchmark:
        "benchmarks/markadapters/{sample_id}.{readgroup}.txt"
    message:
        "Adding XT tags to {wildcards.sample_id} RGID: {wildcards.readgroup}. This marks Illumina "
        "Adapters and allows them to be removed in later steps."
    shell:
        "gatk --java-options {config[standard_java_opt]} MarkIlluminaAdapters \
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
## ie. restore all the original readgroup tags
##
## See: https://gatkforums.broadinstitute.org/gatk/discussion/6483/how-to-map-and-clean-up-short-read-sequence-data-efficiently
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule samtofastq_bwa_mergebamalignment:
    input:
        "markadapters/{sample_id}/{sample_id}.{readgroup}.revertsam.markadapters.bam"
    output:
        temp("bwa/{sample_id}/{sample_id}.{readgroup}.realn.bam")
    threads:
        CLUSTER_META["samtofastq_bwa_mergebamalignment"]["ppn"]
    log: 
        "logs/samtofastq_bwa_mergebamalignment/{sample_id}.{readgroup}.log"
    benchmark:
        "benchmarks/revertsam/{sample_id}.{readgroup}.txt"
    message:
        "BAM to FASTQ --> BWA-MEM --> Merge BAM Alignment.\n"
        "Sample: {wildcards.sample_id}\n"
        "Readgroup: {wildcards.readgroup}\n"
        "The first step converts the reverted BAM to an interleaved FASTQ, removing Illumina "
        "adapters. The output is then piped to BWA-MEM and aligned. Aligned reads are merged "
        "with the original pre-aligned BAM to preserve original metadata, including read groups"
    shell:
        "gatk --java-options {config[samtofastq_java_opt]} SamToFastq \
            --INPUT={input} \
            --FASTQ=/dev/stdout \
            --CLIPPING_ATTRIBUTE=XT \
            --CLIPPING_ACTION=2 \
            --INTERLEAVE=true \
            --NON_PF=true \
            --VALIDATION_STRINGENCY=SILENT \
            --TMP_DIR={config[tempdir]} | \
        bwa mem -M -t {threads} -p {config[reference_fasta]} /dev/stdin | \
        gatk --java-options {config[mergebamalignment_java_opt]} MergeBamAlignment \
            --REFERENCE_SEQUENCE={config[reference_fasta]} \
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
        dynamic("bwa/{sample_id}/{sample_id}.{readgroup}.realn.bam") 
    output:
        bam = "markduplicates/{sample_id}.realn.dedup.bam",
        metrics = "markduplicates/{sample_id}.metrics.txt"
    threads:
        CLUSTER_META["markduplicates"]["ppn"]
    log:
        "logs/markduplicates/{sample_id}.log"
    benchmark:
        "benchmarks/markduplicates/{sample_id}.txt"
    message:
        "Readgroup-specific BAM files are combined into a single BAM for {wildcards.sample_id}. "
        "Potential PCR duplicates are marked."
    run:
        multi_input = " ".join(["--INPUT=" + s for s in input])
        shell("gatk --java-options {config[standard_java_opt]} MarkDuplicates \
            {multi_input} \
            --OUTPUT={output.bam} \
            --METRICS_FILE={output.metrics} \
            --CREATE_INDEX=true \
            2> {log}")

rule baserecalibrator:
    input:
        "markduplicates/{sample_id}.realn.dedup.bam"
    output:
        "bqsr/{sample_id}.bqsr.txt"
    threads:
        CLUSTER_META["baserecalibrator"]["ppn"]
    log:
        "logs/bqsr/{sample_id}.recal.log"
    benchmark:
        "benchmarks/bqsr/{sample_id}.recal.txt"
    message:
        "Calculating base recalibration scores for {wildcards.sample_id}."
    shell:
        "gatk --java-options {config[standard_java_opt]} BaseRecalibrator \
            -R {config[reference_fasta]} \
            -I {input} \
            -O {output} \
            --known-sites {config[gnomad_vcf]}"

rule applybqsr:
    input:
        bam = "markduplicates/{sample_id}.realn.dedup.bam",
        bqsr = "bqsr/{sample_id}.bqsr.txt"
    output:
        protected("bqsr/{sample_id}.realn.dedup.bqsr.bam")
    threads:
        CLUSTER_META["applybqsr"]["ppn"]
    log:
        "logs/bqsr/{sample_id}.apply.log"
    benchmark:
        "benchmarks/bqsr/{sample_id}.apply.txt"
    message:
        "Applying base recalibration scores to {wildcards.sample_id} and generating final "
        "aligned BAM file"
    shell:
        "gatk --java-options {config[standard_java_opt]} ApplyBQSR \
            -R {config[reference_fasta]} \
            -I {input.bam} \
            -OQ true \
            -O {output} \
            -bqsr {input.bqsr}"

# rule coverage:
#     input:
#         "bqsr/{sample_id}.realn.dedup.bqsr.bam"
#     threads:
#         CLUSTER_META["coverage"]["ppn"]
#     log:
#         "logs/coverage/{sample_id}.log"
#     shell:

## END ##
