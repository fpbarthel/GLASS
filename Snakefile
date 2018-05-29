## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Snakefile for GLASS-WG pipeline
## Floris Barthel 2018
## Development branch
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

import pandas as pd
import itertools
import os

## Touch function taken from stackoverflow
## Link: https://stackoverflow.com/questions/1158076/implement-touch-using-python
def touch(fname, mode=0o666, dir_fd=None, **kwargs):
    flags = os.O_CREAT | os.O_APPEND
    with os.fdopen(os.open(fname, flags=flags, mode=mode, dir_fd=dir_fd)) as f:
        os.utime(f.fileno() if os.utime in os.supports_fd else fname,
            dir_fd=None if os.supports_fd else dir_fd, **kwargs)


configfile: "config.yaml"

## Although this statement goes against all coding conventions, we want it here because we want to run
## everything on a temporary storage while we keep this script safe on a permanent drive
workdir: config["workdir"]

## GDC token file for authentication
KEYFILE     = config["gdc_token"]

## Metadata
SAMPLES_META    = json.load(open(config["sample_json"]))
CLUSTER_META    = json.load(open(config["cluster_json"]))

## JSON processing

BAM_FILES = {}
BAM_FILES_UUIDS = {}
BAM_READGROUPS = {}
ALL_READGROUPS = {}
READGROUP_SAMPLE = {}
FQ_FILES = {}
SAMPLES = []

RGPL = {}
RGPU = {}
RGLB = {}
RGDT = {}
RGSM = {}
RGCN = {}

for case in SAMPLES_META:
    for sample in case["samples"]:
        SAMPLES.append(sample["sample_id"])
        FQ_FILES[sample["sample_id"]] = {}
        RGPL[sample["sample_id"]] = {}
        RGPU[sample["sample_id"]] = {}
        RGLB[sample["sample_id"]] = {}
        RGDT[sample["sample_id"]] = {}
        RGSM[sample["sample_id"]] = {}
        RGCN[sample["sample_id"]] = {}
        for file in sample["files"]:
            if file["file_format"] == "BAM":
                BAM_FILES[sample["sample_id"]] = file["file_name"]
                BAM_FILES_UUIDS[sample["sample_id"]] = file["file_uuid"]
                BAM_READGROUPS[sample["sample_id"]] = [readgroup["rg_ID"] for readgroup in file["readgroups"]]
                ALL_READGROUPS[sample["sample_id"]] = [readgroup["rg_ID"] for readgroup in file["readgroups"]]
                for readgroup in file["readgroups"]:
                    READGROUP_SAMPLE[readgroup["rg_ID"]] = sample["sample_id"]
            if file["file_format"] == "FQ":
                FQ_FILES[sample["sample_id"]][file["readgroups"][0]["rg_ID"]] = file["file_name"].split(",")
                RGPL[sample["sample_id"]][file["readgroups"][0]["rg_ID"]] = file["readgroups"][0]["rg_PL"]
                RGPU[sample["sample_id"]][file["readgroups"][0]["rg_ID"]] = file["readgroups"][0]["rg_PU"]
                RGLB[sample["sample_id"]][file["readgroups"][0]["rg_ID"]] = file["readgroups"][0]["rg_LB"]
                RGDT[sample["sample_id"]][file["readgroups"][0]["rg_ID"]] = file["readgroups"][0]["rg_DT"]
                RGSM[sample["sample_id"]][file["readgroups"][0]["rg_ID"]] = file["readgroups"][0]["rg_SM"]
                RGCN[sample["sample_id"]][file["readgroups"][0]["rg_ID"]] = file["readgroups"][0]["rg_CN"]
                if sample["sample_id"] in ALL_READGROUPS:
                    ALL_READGROUPS[sample["sample_id"]].append(file["readgroups"][0]["rg_ID"])
                else:
                    ALL_READGROUPS[sample["sample_id"]] = [ file["readgroups"][0]["rg_ID"] ]

FQ_FILES = dict((sample,readgroup) for sample,readgroup in FQ_FILES.items() if len(readgroup)>0)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Master rule
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule all:
    input: "results/qc/multiqc/multiqc_report.html"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Download only rule
## Run snakemake with 'snakemake download_only' to activate
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule download_only:
    input: expand("download/{uuid}/{file}", zip, uuid=BAM_FILES_UUIDS.values(), file=BAM_FILES.values())
 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Download BAM file from GDC
## GDC key needs to be re-downloaded and updated from time to time
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
            -t {config[gdc_token]} \
            {wildcards.uuid} \
            2> {log}"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## RevertSAM and FASTQ-2-uBAM both output uBAM files to the same directory
## Snakemake needs to be informed which rule takes presidence
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

#ruleorder: revertsam > fq2ubam

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Revert GDC-aligned legacy BAM to unaligned SAM file
## Clear BAM attributes, such as re-calibrated base quality scores, alignment information
## Moreover, BAM file is split into per-readgroup BAM files
## See: https://gatkforums.broadinstitute.org/gatk/discussion/6484
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## NOTE THAT
## --MAX_DISCARD_FRACTION=0.05
## SHOULD BE USED FOR TESTING PURPOSES ONLY!!
## This has now been added to "config[revertsam_extra_args]"
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## 05/29/18
## Because of a snakemake limitation with dynamic output (as this rule outputs an variable)
## number of files depending on the number of readgroups, we have added python code to
## create empty files ("touch") for those readgroups not relevant for this analysis
## Issue: https://bitbucket.org/snakemake/snakemake/issues/865/pre-determined-dynamic-output
## Unlikely this will get fixed any time soon
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule revertsam:
    input:
        lambda wildcards: "download/{uuid}/{file}".format(uuid=BAM_FILES_UUIDS[wildcards.sample_id], file=BAM_FILES[wildcards.sample_id])
    output:
        map = "results/ubam/{sample_id}/{sample_id}.output_map.txt",
        bams = temp(expand("results/ubam/{{sample_id}}/{{sample_id}}.{rg}.bam", rg=list(itertools.chain.from_iterable(BAM_READGROUPS.values()))))
    params:
        dir = "results/ubam/{sample_id}"
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
    run:
        ## Create a readgroup name / filename mapping file
        rgmap = pd.DataFrame(
            {
                "READ_GROUP_ID": BAM_READGROUPS[wildcards["sample_id"]],
                "OUTPUT": ["results/ubam/{sample}/{sample}.{rg}.bam".format(sample=wildcards["sample_id"], rg=rg) for rg in BAM_READGROUPS[wildcards["sample_id"]]]
            },
            columns = ["READ_GROUP_ID", "OUTPUT"]
        )
        rgmap.to_csv(output["map"], sep="\t", index=False)

        ## Create empty files ("touch") for readgroups not in this BAM file
        ## Workaround for issue documented here: https://bitbucket.org/snakemake/snakemake/issues/865/pre-determined-dynamic-output
        other_rg_f = ["results/ubam/{sample}/{sample}.{rg}.bam".format(sample=wildcards["sample_id"],rg=rg) for sample, rgs in BAM_READGROUPS.items() for rg in rgs if sample not in wildcards["sample_id"]]
        for f in other_rg_f:
            touch(f)

        shell("gatk --java-options {config[standard_java_opt]} RevertSam \
            --INPUT={input} \
            --OUTPUT_BY_READGROUP=true \
            --OUTPUT_BY_READGROUP_FILE_FORMAT=bam \
            --OUTPUT_MAP={output.map} \
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
            --SORT_ORDER=queryname {config[revertsam_extra_args]}\
            --TMP_DIR={config[tempdir]} \
            2> {log}")

# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# ## Convert from FASTQ pair to uBAM
# ## This step eases follow up steps
# ## See: https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule fq2ubam:
    input:
        R1 = lambda wildcards: "fastq/{sample}/{file}".format(sample=wildcards.sample_id, file=FQ_FILES[wildcards.sample_id][wildcards.readgroup][0]),
        R2 = lambda wildcards: "fastq/{sample}/{file}".format(sample=wildcards.sample_id, file=FQ_FILES[wildcards.sample_id][wildcards.readgroup][1])
    output:
        temp("results/ubam/{sample_id}/{sample_id}.{readgroup}.bam")
    params:
        RGID = lambda wildcards: wildcards.readgroup,
        RGPL = lambda wildcards: RGPL[wildcards.sample_id][wildcards.readgroup],
        RGPU = lambda wildcards: RGPU[wildcards.sample_id][wildcards.readgroup],
        RGLB = lambda wildcards: RGLB[wildcards.sample_id][wildcards.readgroup],
        RGDT = lambda wildcards: RGDT[wildcards.sample_id][wildcards.readgroup],
        RGSM = lambda wildcards: RGSM[wildcards.sample_id][wildcards.readgroup],
        RGCN = lambda wildcards: RGCN[wildcards.sample_id][wildcards.readgroup]
    log:
        "logs/fq2ubam/{sample_id}.{readgroup}.log"
    threads:
        CLUSTER_META["fq2ubam"]["ppn"]
    benchmark:
        "benchmarks/fq2ubam/{sample_id}.{readgroup}.txt"
    message:
        "Converting FASTQ file {wildcards.sample_id}, {wildcards.readgroup} "
        "to uBAM format."
    run:
        ## ISODATE=`date +%Y-%m-%dT%H:%M:%S%z`; \
        shell("gatk --java-options {config[standard_java_opt]} FastqToSam \
            --FASTQ={input.R1} \
            --FASTQ2={input.R2} \
            --OUTPUT={output} \
            --READ_GROUP_NAME=\"{params.RGID}\" \
            --PLATFORM_UNIT=\"{params.RGPU}\" \
            --SAMPLE_NAME=\"{params.RGSM}\" \
            --PLATFORM=\"{params.RGPL}\" \
            --LIBRARY_NAME=\"{params.RGLB}\" \
            --SEQUENCING_CENTER=\"{params.RGCN}\" \
            --SORT_ORDER=queryname \
            2> {log}")
        #            --RUN_DATE=\"{params.RGDT}\" \

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run FASTQC on uBAM
## URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule fastqc:
    input:
        "results/ubam/{sample_id}/{sample_id}.{readgroup}.bam"
    output:
        "results/qc/{sample_id}/{sample_id}.{readgroup}_fastqc.html"
    params:
        dir = "results/qc/{sample_id}"
    threads:
        CLUSTER_META["fastqc"]["ppn"]
    log:
        "logs/fastqc/{sample_id}.{readgroup}.log"
    benchmark:
        "benchmarks/fastqc/{sample_id}.{readgroup}.txt"
    message:
        "Running FASTQC for {wildcards.sample_id} {wildcards.readgroup}"
    shell:
        "fastqc \
            -o {params.dir} \
            -f bam \
            {input} \
            2> {log}"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Mark Illumina Adapters
## Add XT tag to read records to mark the 5' start position of adapter sequences
## Adapter sequences are then removed by subsequent steps
## See: https://gatkforums.broadinstitute.org/gatk/discussion/6483/how-to-map-and-clean-up-short-read-sequence-data-efficiently
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule markadapters:
    input:
        "results/ubam/{sample_id}/{sample_id}.{readgroup}.bam"
    output:
        bam = temp("results/markadapters/{sample_id}/{sample_id}.{readgroup}.revertsam.markadapters.bam"),
        metric = "results/markadapters/{sample_id}/{sample_id}.{readgroup}.markadapters.metrics.txt"
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
        "results/markadapters/{sample_id}/{sample_id}.{readgroup}.revertsam.markadapters.bam"
    output:
        temp("results/bwa/{sample_id}/{sample_id}.{readgroup}.realn.bam")
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
        lambda wildcards: expand("results/bwa/{sample}/{sample}.{rg}.realn.bam", sample=wildcards.sample_id, rg=ALL_READGROUPS[wildcards.sample_id]) #get_readgroup_BAMs_sample
    output:
        bam = temp("results/markduplicates/{sample_id}.realn.dedup.bam"),
        metrics = "results/markduplicates/{sample_id}.metrics.txt"
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

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Recalibrate base quality scores
## This steps computes a bare recalibration table
## See: https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
rule baserecalibrator:
    input:
        "results/markduplicates/{sample_id}.realn.dedup.bam"
    output:
        "results/bqsr/{sample_id}.bqsr.txt"
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
            --known-sites {config[gnomad_vcf]} \
            2> {log}"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Apply BQSR
## This step applies a base recalibration table to an input BAM file
## Formerly "PrintReads"
## See: https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
rule applybqsr:
    input:
        bam = "results/markduplicates/{sample_id}.realn.dedup.bam",
        bqsr = "results/bqsr/{sample_id}.bqsr.txt"
    output:
        protected("results/bqsr/{sample_id}.realn.dedup.bqsr.bam")
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
            -bqsr {input.bqsr} \
            --create-output-bam-md5 true \
            2> {log}"

rule wgsmetrics:
    input:
        "results/bqsr/{sample_id}.realn.dedup.bqsr.bam"
    output:
        "results/qc/{sample_id}.WgsMetrics.txt"
    threads:
        CLUSTER_META["wgsmetrics"]["ppn"]
    log:
        "logs/wgsmetrics/{sample_id}.WgsMetrics.log"
    benchmark:
        "benchmarks/wgsmetrics/{sample_id}.WgsMetrics.txt"
    message:
        "Computing WGS Metrics for {wildcards.sample_id}"
    shell:
        "gatk --java-options {config[standard_java_opt]} CollectWgsMetrics \
            -R {config[reference_fasta]} \
            -I {input} \
            -O {output} \
            --USE_FAST_ALGORITHM true \
            2> {log}"

rule validatebam:
    input:
        "results/bqsr/{sample_id}.realn.dedup.bqsr.bam"
    output:
        "results/qc/{sample_id}.ValidateSamFile.txt"
    threads:
        CLUSTER_META["validatebam"]["ppn"]
    log:
        "logs/validatebam/{sample_id}.ValidateSamFile.log"
    benchmark:
        "benchmarks/validatebam/{sample_id}.ValidateSamFile.txt"
    message:
        "Validating BAM file {wildcards.sample_id}"
    shell:
        "gatk --java-options {config[standard_java_opt]} ValidateSamFile \
            -I {input} \
            -O {output} \
            -M SUMMARY \
            2> {log}"

rule multiqc:
    input:
        expand("results/qc/{sample}.ValidateSamFile.txt", sample=ALL_READGROUPS.keys()),
        expand("results/qc/{sample}.WgsMetrics.txt", sample=ALL_READGROUPS.keys()),
        lambda wildcards: ["results/qc/{sample}/{sample}.{rg}_fastqc.html".format(sample=sample, rg=readgroup)
          for sample, readgroups in ALL_READGROUPS.items()
          for readgroup in readgroups] 
    output:
        "results/qc/multiqc/multiqc_report.html"
    params:
        dir = "results/qc/multiqc"
    threads:
        CLUSTER_META["multiqc"]["ppn"]
    log:
        "logs/multiqc/multiqc.log"
    benchmark:
        "benchmarks/multiqc/multiqc.txt"
    message:
        "Running MultiQC"
    shell:
        "multiqc -o {params.dir} {config[workdir]}; cp -R {params.dir}/* {config[html_dir]}"

## END ##
