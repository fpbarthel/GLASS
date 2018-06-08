## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Snakefile for GLASS-WG pipeline
## Authors: Floris Barthel, Samir Amin
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

## Turn an unnamed list of dicts into a nammed list of dicts
## Taken from stackoverflow
## https://stackoverflow.com/questions/4391697/find-the-index-of-a-dict-within-a-list-by-matching-the-dicts-value
def build_dict(seq, key):
    return dict((d[key], dict(d, index=index)) for (index, d) in enumerate(seq))

## Although this statement goes against all coding conventions, we want it here because we want to run
## everything on a temporary storage while we keep this script safe on a permanent drive
## TEMPORARY
configfile: "config.yaml"
workdir: config["workdir"]

## GDC token file for authentication
KEYFILE     = config["gdc_token"]

## Metadata
SAMPLES_META    = json.load(open(config["sample_json"]))
CLUSTER_META    = json.load(open(config["cluster_json"]))

CASES 		= json.load(open(config["cases_json"]))
SAMPLES 	= json.load(open(config["samples_json"]))
ALIQUOTS 	= json.load(open(config["aliquots_json"]))
FILES 		= json.load(open(config["files_json"]))
READGROUPS 	= json.load(open(config["readgroups_json"]))
PAIRS 		= json.load(open(config["pairs_json"]))

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## JSON processing
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# @sbamin Unless already implemented, we should be explicitly checking json input for 
# 1. non-empty variables, and 
# 2. unique RGID and RGPU but an identical RGSM tags, e.g., https://github.com/TheJacksonLaboratory/glass_wgs_alignment/blob/d72fb20659bd20fddf952d331533b9ffd88d446e/runner/preprocess_fqs.R#L25 
# We can either check it upfront while making json or more preferable to check just before snakemake submits a workflow per case or sample.
# That way, snakemake should STOP with error or emit WARN for non-compliant RG format. 
# This is more practical if input is FQ and not BAM unless we already have RG info for BAM file.

### NOTE NEED TO SPERATE BAM AND FASTQ READGROUPS

## Validate CASES JSON
## CASES should be unique
## CHECK THAT ALL CASE_ID VALUES IN CASES ARE UNIQUE
## TO-DO


## Validate FILES JSON
## FILES -> FILE_UUID should be unique
## CHECK THAT ALL FILE_UUID VALUES IN FILES ARE UNIQUE
## Check that input files exist
## TO-DO


## Validate PAIRS JSON
## PAIR -> PAIR_ID should be unique
## CHECK THAT ALL PAIR_ID VALUES IN PAIR ARE UNIQUE
## IN PROGRESS


## CASES -> DICT
CASES_DICT = build_dict(CASES, "case_id")


## FILES -> DICT
FILES_DICT = build_dict(FILES, "file_uuid")


## Pair IDs are unique, PAIRS -> DICT
PAIRS_DICT = build_dict(PAIRS, "pair_id")


## Aliquot IDs and BAM files map 1:1

ALIQUOT_TO_BAM_PATH = {}
for file in FILES:
    if file["file_format"] == "BAM":
        ALIQUOT_TO_BAM_PATH[ file["aliquot_id"] ] = file["file_path"]


## Aliquots and RGIDs map 1:many

ALIQUOT_TO_RGID = {}        
for readgroup in READGROUPS:
    if readgroup["aliquot_id"] not in ALIQUOT_TO_RGID:
        ALIQUOT_TO_RGID[ readgroup["aliquot_id"] ] = [ readgroup["RGID"] ]
    else:
        ALIQUOT_TO_RGID[ readgroup["aliquot_id"] ].append(readgroup["RGID"])


## Batches and normal aliquot IDs map 1:many
## Normal aliquot IDs are repeated across multiple pairs from same case
## Each pair has one normal and one tumor

BATCH_TO_NORMAL = {}
for pair in PAIRS:
    pair["project_id"] = CASES_DICT[ pair["case_id"] ]["project_id"]
    PAIRS_DICT[ pair["pair_id"] ]["project_id"] = pair["project_id"]
    if pair["project_id"] not in BATCH_TO_NORMAL:
        BATCH_TO_NORMAL[ pair["project_id"] ] = [ pair["normal_aliquot_id"] ]
    elif pair["normal_aliquot_id"] not in BATCH_TO_NORMAL[ pair["project_id"] ]:
        BATCH_TO_NORMAL[ pair["project_id"] ].append(pair["normal_aliquot_id"])
        

## Readgroup information and 
## Aliquots and RGIDs map 1:many
## RGIDs are unique within an aliquot
## Aliquot IDs and fastQ files map 1:many

ALIQUOT_TO_READGROUP = {} 
for readgroup in READGROUPS:
    if readgroup["aliquot_id"] not in ALIQUOT_TO_READGROUP:
        ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ] = { readgroup["RGID"] : readgroup }
    else:
        ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ][ readgroup["RGID"] ] = readgroup
    ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ][ readgroup["RGID"] ]["file_path"] = FILES_DICT[ ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ][ readgroup["RGID"] ]["file_uuid"] ]["file_path"]
    ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ][ readgroup["RGID"] ]["file_format"] = FILES_DICT[ ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ][ readgroup["RGID"] ]["file_uuid"] ]["file_format"]

## List of scatterlist items to iterate over
## Each Mutect2 run spawns 50 jobs based on this scatterlist

WGS_SCATTERLIST = ["temp_{num}_of_50".format(num=str(j+1).zfill(4)) for j in range(50)]

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
   input: expand("{file}", file=ALIQUOT_TO_BAM_PATH.values())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## SNV rule
## Run snakemake with target 'snv'
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule snv:
    input: expand("results/vep/{pair_id}.filtered2.anno.maf", pair_id=PAIRS_DICT.keys())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Download BAM file from GDC
## GDC key needs to be re-downloaded and updated from time to time
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

# rule download:
#     output:
#         "download/{uuid}/{filename}.bam"
#     threads:
#         CLUSTER_META["download"]["ppn"]
#     message:
#         "Downloading from GDC\n"
#         "UUID {wildcards.uuid}\n"
#         "File {wildcards.filename}"
#     log:
#         "logs/download/{uuid}.log"
#     benchmark:
#         "benchmarks/download/{uuid}.txt"
#     shell:
#         "gdc-client download \
#             -d download \
#             -n {threads} \
#             -t {config[gdc_token]} \
#             {wildcards.uuid} \
#             > {log} 2>&1"

# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# ## RevertSAM and FASTQ-2-uBAM both output uBAM files to the same directory
# ## Snakemake needs to be informed which rule takes presidence
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# #ruleorder: revertsam > fq2ubam

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

# rule revertsam: ## 5/8 this rule is causing trouble, comment for now
#     input:
#         lambda wildcards: ALIQUOT_TO_BAM_PATH[wildcards.aliquot_id]
#     output:
#         map = "results/ubam/{aliquot_id}/{aliquot_id}.output_map.txt",
#         bams = temp(expand("results/ubam/{{aliquot_id}}/{{aliquot_id}}.{rg}.bam", rg=list(itertools.chain.from_iterable(ALIQUOT_TO_RGID.values()))))
#     params:
#         dir = "results/ubam/{aliquot_id}",
#         mem = CLUSTER_META["revertsam"]["mem"]
#     log: 
#         "logs/revertsam/{aliquot_id}.log"
#     threads:
#         CLUSTER_META["revertsam"]["ppn"]
#     benchmark:
#         "benchmarks/revertsam/{aliquot_id}.txt"
#     message:
#         "Reverting sample back to unaligned BAM file, stripping any previous "
#         "pre-processing and restoring original base quality scores. Output files are split "
#         "by readgroup.\n"
#         "Sample: {wildcards.aliquot_id}"
#     run:
#         ## Create a readgroup name / filename mapping file
#         rgmap = pd.DataFrame(
#             {
#                 "READ_GROUP_ID": BAM_READGROUPS[wildcards["aliquot_id"]],
#                 "OUTPUT": ["results/ubam/{aliquot_id}/{aliquot_id}.{rg}.bam".format(aliquot_id=wildcards["aliquot_id"], rg=rg) for rg in BAM_READGROUPS[wildcards["aliquot_id"]]]
#             },
#             columns = ["READ_GROUP_ID", "OUTPUT"]
#         )
#         rgmap.to_csv(output["map"], sep="\t", index=False)

#         ## Create empty files ("touch") for readgroups not in this BAM file
#         ## Workaround for issue documented here: https://bitbucket.org/snakemake/snakemake/issues/865/pre-determined-dynamic-output
#         other_rg_f = ["results/ubam/{aliquot_id}/{aliquot_id}.{rg}.bam".format(aliquot_id=wildcards["aliquot_id"],rg=rg) for sample, rgs in BAM_READGROUPS.items() for rg in rgs if sample not in wildcards["aliquot_id"]]
#         for f in other_rg_f:
#             touch(f)

#         shell("gatk --java-options -Xmx{params.mem}g RevertSam \
#             --INPUT={input} \
#             --OUTPUT_BY_READGROUP=true \
#             --OUTPUT_BY_READGROUP_FILE_FORMAT=bam \
#             --OUTPUT_MAP={output.map} \
#             --RESTORE_ORIGINAL_QUALITIES=true \
#             --VALIDATION_STRINGENCY=SILENT \
#             --ATTRIBUTE_TO_CLEAR=AS \
#             --ATTRIBUTE_TO_CLEAR=FT \
#             --ATTRIBUTE_TO_CLEAR=CO \
#             --ATTRIBUTE_TO_CLEAR=XT \
#             --ATTRIBUTE_TO_CLEAR=XN \
#             --ATTRIBUTE_TO_CLEAR=OC \
#             --ATTRIBUTE_TO_CLEAR=OP \
#             --SANITIZE=true \
#             --SORT_ORDER=queryname {config[revertsam_extra_args]}\
#             --TMP_DIR={config[tempdir]} \
#             > {log} 2>&1")

# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# ## Convert from FASTQ pair to uBAM
# ## This step eases follow up steps
# ## See: https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# rule fq2ubam:
#     input:
#         R1 = lambda wildcards: "fastq/{sample}/{file}".format(sample=wildcards.aliquot_id, file=FQ_FILES[wildcards.aliquot_id][wildcards.readgroup][0]),
#         R2 = lambda wildcards: "fastq/{sample}/{file}".format(sample=wildcards.aliquot_id, file=FQ_FILES[wildcards.aliquot_id][wildcards.readgroup][1])
#     output:
#         temp("results/ubam/{aliquot_id}/{aliquot_id}.{readgroup}.bam")
#     params:
#         RGID = lambda wildcards: wildcards.readgroup,
#         RGPL = lambda wildcards: RGPL[wildcards.aliquot_id][wildcards.readgroup],
#         RGPU = lambda wildcards: RGPU[wildcards.aliquot_id][wildcards.readgroup],
#         RGLB = lambda wildcards: RGLB[wildcards.aliquot_id][wildcards.readgroup],
#         RGDT = lambda wildcards: RGDT[wildcards.aliquot_id][wildcards.readgroup],
#         RGSM = lambda wildcards: RGSM[wildcards.aliquot_id][wildcards.readgroup],
#         RGCN = lambda wildcards: RGCN[wildcards.aliquot_id][wildcards.readgroup],
#         mem = CLUSTER_META["fq2ubam"]["mem"]
#     log:
#         "logs/fq2ubam/{aliquot_id}.{readgroup}.log"
#     threads:
#         CLUSTER_META["fq2ubam"]["ppn"]
#     benchmark:
#         "benchmarks/fq2ubam/{aliquot_id}.{readgroup}.txt"
#     message:
#         "Converting FASTQ file to uBAM format\n"
#         "Sample: {wildcards.aliquot_id}\n"
#         "Readgroup: {wildcards.readgroup}"
#     run:
#         ## ISODATE=`date +%Y-%m-%dT%H:%M:%S%z`; \
#         shell("gatk --java-options -Xmx{params.mem}g FastqToSam \
#             --FASTQ={input.R1} \
#             --FASTQ2={input.R2} \
#             --OUTPUT={output} \
#             --READ_GROUP_NAME=\"{params.RGID}\" \
#             --PLATFORM_UNIT=\"{params.RGPU}\" \
#             --SAMPLE_NAME=\"{params.RGSM}\" \
#             --PLATFORM=\"{params.RGPL}\" \
#             --LIBRARY_NAME=\"{params.RGLB}\" \
#             --SEQUENCING_CENTER=\"{params.RGCN}\" \
#             --SORT_ORDER=queryname \
#             > {log} 2>&1")
#         #            --RUN_DATE=\"{params.RGDT}\" \

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run FASTQC on uBAM
## URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# @sbamin Besides html, fastqc should provide a tab-delimited output, 
# e.g., https://gist.github.com/chapmanb/3953983, 
# https://gitlab.com/gmapps/railab_chipseq/tree/master/scripts 
# A step further, that can be programmatically added to emit WARN or STOP if converted FQ fails to PASS base filters, e.g., per base or tile seq quality.

rule fastqc:
    input:
        "results/ubam/{aliquot_id}/{aliquot_id}.{readgroup}.bam"
    output:
        "results/qc/{aliquot_id}/{aliquot_id}.{readgroup}_fastqc.html"
    params:
        dir = "results/qc/{aliquot_id}",
        mem = CLUSTER_META["fastqc"]["mem"]
    threads:
        CLUSTER_META["fastqc"]["ppn"]
    log:
        "logs/fastqc/{aliquot_id}.{readgroup}.log"
    benchmark:
        "benchmarks/fastqc/{aliquot_id}.{readgroup}.txt"
    message:
        "Running FASTQC\n"
        "Sample: {wildcards.aliquot_id}\n"
        "Readgroup: {wildcards.readgroup}"
    shell:
        "fastqc \
            -o {params.dir} \
            -f bam \
            {input} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Mark Illumina Adapters
## Add XT tag to read records to mark the 5' start position of adapter sequences
## Adapter sequences are then removed by subsequent steps
## See: https://gatkforums.broadinstitute.org/gatk/discussion/6483/how-to-map-and-clean-up-short-read-sequence-data-efficiently
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule markadapters:
    input:
        "results/ubam/{aliquot_id}/{aliquot_id}.{readgroup}.bam"
    output:
        bam = temp("results/markadapters/{aliquot_id}/{aliquot_id}.{readgroup}.revertsam.markadapters.bam"),
        metric = "results/markadapters/{aliquot_id}/{aliquot_id}.{readgroup}.markadapters.metrics.txt"
    params:
        mem = CLUSTER_META["markadapters"]["mem"]
    threads:
        CLUSTER_META["markadapters"]["ppn"]
    params:
        mem = CLUSTER_META["markadapters"]["mem"]
    log: 
        dynamic("logs/markadapters/{aliquot_id}.{readgroup}.log")
    benchmark:
        "benchmarks/markadapters/{aliquot_id}.{readgroup}.txt"
    message:
        "Adding XT tags. This marks Illumina Adapters and allows them to be removed in later steps\n"
        "Sample: {wildcards.aliquot_id}\n"
        "Readgroup: {wildcards.readgroup}"
    shell:
        "gatk --java-options -Xmx{params.mem}g MarkIlluminaAdapters \
            --INPUT={input} \
            --OUTPUT={output.bam} \
            --METRICS={output.metric} \
            --TMP_DIR={config[tempdir]} \
            > {log} 2>&1"

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
##
## Update 06/01: Added markadapter metrics as input even though not required. Metrics file
## is saved all the way at the end of markadapters step, and adding it makes sure that the
## input data is good
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# @sbamin I am new to this step. Are read groups, including RGID per read incorporated in final merge bam
# I am more familiar with bwa mem -R '@RG\tID:foo\tSM:bar' format. 
# https://github.com/TheJacksonLaboratory/glass_wgs_alignment/blob/d72fb20659bd20fddf952d331533b9ffd88d446e/runner/preprocess_fqs.R#L103
# @sbamin If original bam file (and so coverted fastq) are using ref genome with chr prefix for chromosomes, 
# does workflow take care of removing chr prefix (or vice versa) while aligning to ref genome from GATK b37 legacy bundle?

rule samtofastq_bwa_mergebamalignment:
    input:
        bam = "results/markadapters/{aliquot_id}/{aliquot_id}.{readgroup}.revertsam.markadapters.bam",
        metric = "results/markadapters/{aliquot_id}/{aliquot_id}.{readgroup}.markadapters.metrics.txt"
    output:
        bam = temp("results/bwa/{aliquot_id}/{aliquot_id}.{readgroup}.realn.bam"),
        bai = temp("results/bwa/{aliquot_id}/{aliquot_id}.{readgroup}.realn.bai")
    threads:
        CLUSTER_META["samtofastq_bwa_mergebamalignment"]["ppn"]
    params:
        mem = CLUSTER_META["samtofastq_bwa_mergebamalignment"]["mem"]
    log: 
        "logs/samtofastq_bwa_mergebamalignment/{aliquot_id}.{readgroup}.log"
    benchmark:
        "benchmarks/revertsam/{aliquot_id}.{readgroup}.txt"
    message:
        "BAM to FASTQ --> BWA-MEM --> Merge BAM Alignment.\n"
        "The first step converts the reverted BAM to an interleaved FASTQ, removing Illumina "
        "adapters. The output is then piped to BWA-MEM and aligned. Aligned reads are merged "
        "with the original pre-aligned BAM to preserve original metadata, including read groups.\n"
        "Sample: {wildcards.aliquot_id}\n"
        "Readgroup: {wildcards.readgroup}"
    shell:
        "gatk --java-options {config[samtofastq_java_opt]} SamToFastq \
            --INPUT={input.bam} \
            --FASTQ=/dev/stdout \
            --CLIPPING_ATTRIBUTE=XT \
            --CLIPPING_ACTION=2 \
            --INTERLEAVE=true \
            --NON_PF=true \
            --TMP_DIR={config[tempdir]} | \
         bwa mem -M -t {threads} -p {config[reference_fasta]} /dev/stdin | \
         gatk --java-options {config[mergebamalignment_java_opt]} MergeBamAlignment \
            --ALIGNED_BAM=/dev/stdin \
            --UNMAPPED_BAM={input.bam} \
            --OUTPUT={output.bam} \
            --REFERENCE_SEQUENCE={config[reference_fasta]} \
            --CREATE_INDEX=true \
            --ADD_MATE_CIGAR=true \
            --CLIP_ADAPTERS=false \
            --CLIP_OVERLAPPING_READS=true \
            --INCLUDE_SECONDARY_ALIGNMENTS=true \
            --MAX_INSERTIONS_OR_DELETIONS=-1 \
            --PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
            --ATTRIBUTES_TO_RETAIN=XS \
            --TMP_DIR={config[tempdir]} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Mark Duplicates & merge readgroups
## This step marks duplicate reads
## See: https://gatkforums.broadinstitute.org/gatk/discussion/2799
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule markduplicates:
    input:
        lambda wildcards: expand("results/bwa/{sample}/{sample}.{rg}.realn.bam", sample=wildcards.aliquot_id, rg=ALIQUOT_TO_RGID[wildcards.aliquot_id])
    output:
        bam = temp("results/markduplicates/{aliquot_id}.realn.mdup.bam"),
        bai = temp("results/markduplicates/{aliquot_id}.realn.mdup.bai"),
        metrics = "results/markduplicates/{aliquot_id}.metrics.txt"
    params:
        mem = CLUSTER_META["markduplicates"]["mem"]
    threads:
        CLUSTER_META["markduplicates"]["ppn"]
    log:
        "logs/markduplicates/{aliquot_id}.log"
    benchmark:
        "benchmarks/markduplicates/{aliquot_id}.txt"
    message:
        "Readgroup-specific BAM files are combined into a single BAM. "
        "Potential PCR duplicates are marked.\n"
        "Sample: {wildcards.aliquot_id}"
    run:
        multi_input = " ".join(["--INPUT=" + s for s in input])
        shell("gatk --java-options -Xmx{params.mem}g MarkDuplicates \
            {multi_input} \
            --OUTPUT={output.bam} \
            --METRICS_FILE={output.metrics} \
            --CREATE_INDEX=true \
            > {log} 2>&1")

# @sbamin A few notes on IndelRealignment step at annotated link: https://hyp.is/8_20bK-aEeerk1MduBFv6w/gatkforums.broadinstitute.org/gatk/discussion/7847 
# IndelRealinger adds OC:Z tag to realigned reads (GATK Doc # 7156), shifts MAPQ by +10. 
# For high-quality data (>75bp read length, 30x plus, PCR-free lib), and MAPQ > 20, there seems to have minimal impact of +MAPQ10. 
# Perhaps, we should get clarification from GATK team based on GATK Doc # 7847 on pros and cons of keeping IndelRealigner step in the workflow. 
# If there is minimal impact, we should add IndelRealigner. I have not yet seen but OC:Z tag may be required by other indel callers.

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Recalibrate base quality scores
## This steps computes a bare recalibration table
## See: https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule baserecalibrator:
    input:
        "results/markduplicates/{aliquot_id}.realn.mdup.bam"
    output:
        "results/bqsr/{aliquot_id}.bqsr.txt"
    params:
        mem = CLUSTER_META["baserecalibrator"]["mem"]
    threads:
        CLUSTER_META["baserecalibrator"]["ppn"]
    log:
        "logs/bqsr/{aliquot_id}.recal.log"
    benchmark:
        "benchmarks/bqsr/{aliquot_id}.recal.txt"
    message:
        "Calculating base recalibration scores.\n"
        "Sample: {wildcards.aliquot_id}"
    shell:
        "gatk --java-options -Xmx{params.mem}g BaseRecalibrator \
            -R {config[reference_fasta]} \
            -I {input} \
            -O {output} \
            --known-sites {config[gnomad_vcf]} \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Apply BQSR
## This step applies a base recalibration table to an input BAM file
## Formerly "PrintReads"
## See: https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule applybqsr:
    input:
        bam = "results/markduplicates/{aliquot_id}.realn.mdup.bam",
        bqsr = "results/bqsr/{aliquot_id}.bqsr.txt"
    output:
        protected("results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam")
    params:
        mem = CLUSTER_META["applybqsr"]["mem"]
    threads:
        CLUSTER_META["applybqsr"]["ppn"]
    log:
        "logs/bqsr/{aliquot_id}.apply.log"
    benchmark:
        "benchmarks/bqsr/{aliquot_id}.apply.txt"
    message:
        "Applying base recalibration scores and generating final BAM file\n"
        "Sample: {wildcards.aliquot_id}"
    shell:
        "gatk --java-options -Xmx{params.mem}g ApplyBQSR \
            -R {config[reference_fasta]} \
            -I {input.bam} \
            -OQ true \
            -O {output} \
            -bqsr {input.bqsr} \
            --create-output-bam-md5 true \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## WGS Metrics
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Calculate coverage metrics for the WGS BAM file. These metrics are shown by MultiQC
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule wgsmetrics:
    input:
        "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam"
    output:
        "results/qc/{aliquot_id}.WgsMetrics.txt"
    params:
        mem = CLUSTER_META["wgsmetrics"]["mem"]
    threads:
        CLUSTER_META["wgsmetrics"]["ppn"]
    log:
        "logs/wgsmetrics/{aliquot_id}.WgsMetrics.log"
    benchmark:
        "benchmarks/wgsmetrics/{aliquot_id}.WgsMetrics.txt"
    message:
        "Computing WGS Metrics\n"
        "Sample: {wildcards.aliquot_id}"
    shell:
        "gatk --java-options -Xmx{params.mem}g CollectWgsMetrics \
            -R {config[reference_fasta]} \
            -I {input} \
            -O {output} \
            --USE_FAST_ALGORITHM true \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Validate BAM file
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Final check to ensure no errors in final analysis-ready BAM file
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule validatebam:
    input:
        "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam"
    output:
        "results/qc/{aliquot_id}.ValidateSamFile.txt"
    params:
        mem = CLUSTER_META["validatebam"]["mem"]
    threads:
        CLUSTER_META["validatebam"]["ppn"]
    log:
        "logs/validatebam/{aliquot_id}.ValidateSamFile.log"
    benchmark:
        "benchmarks/validatebam/{aliquot_id}.ValidateSamFile.txt"
    message:
        "Validating BAM file\n"
        "Sample: {wildcards.aliquot_id}"
    shell:
        "gatk --java-options -Xmx{params.mem}g ValidateSamFile \
            -I {input} \
            -O {output} \
            -M SUMMARY \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## MultiQC
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## MultiQC (and its Click dependency) can have problwms with Python locale settings
## Fix these by specifying the following options in .bash_profile
## export LC_ALL=en_US.UTF-8
## export LANG=en_US.UTF-8
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule multiqc:
    input:
        expand("results/qc/{sample}.ValidateSamFile.txt", sample=ALIQUOT_TO_RGID.keys()),
        expand("results/qc/{sample}.WgsMetrics.txt", sample=ALIQUOT_TO_RGID.keys()),
        lambda wildcards: ["results/qc/{sample}/{sample}.{rg}_fastqc.html".format(sample=sample, rg=readgroup)
          for sample, readgroups in ALIQUOT_TO_RGID.items()
          for readgroup in readgroups] 
    output:
        "results/qc/multiqc/multiqc_report.html"
    params:
        dir = "results/qc/multiqc",
        mem = CLUSTER_META["samtofastq_bwa_mergebamalignment"]["mem"]
    threads:
        CLUSTER_META["multiqc"]["ppn"]
    log:
        "logs/multiqc/multiqc.log"
    benchmark:
        "benchmarks/multiqc/multiqc.txt"
    message:
        "Running MultiQC"
    shell:
        "multiqc -o {params.dir} {config[workdir]}/results \
            > {log} 2>&1; \
            cp -R {params.dir}/* {config[html_dir]}"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## CallPON
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## This rule uses Mutect2 in tumor-only mode (artifact detection mode) to detect ALL
## variants in a given non-tumor sample
## No germline resource in included because this would exclude these variants from the
## PON
## See: 
## https://gatkforums.broadinstitute.org/gatk/discussion/11136/how-to-call-somatic-mutations-using-gatk4-mutect2
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule callpon:
    input:
        bam = "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam",
        intervallist = lambda wildcards: "{dir}/{interval}/scattered.interval_list".format(dir=config["wgs_scatterdir"], interval=wildcards.interval)
    output:
        temp("results/callpon/{batch}/{aliquot_id}/{aliquot_id}.{interval}.pon.vcf")
    params:
        mem = CLUSTER_META["callpon"]["mem"]
    threads:
        CLUSTER_META["callpon"]["ppn"]
    log:
        "logs/callpon/{batch}.{aliquot_id}.{interval}.log"
    benchmark:
        "benchmarks/callpon/{batch}.{aliquot_id}.{interval}.txt"
    message:
        "Calling Mutect2 in tumor-only mode to build a panel of normals\n"
        "Batch: {wildcards.batch}\n"
        "Sample: {wildcards.aliquot_id}\n"
        "Interval: {wildcards.interval}"
    shell:
        "TUMOR_SM=`samtools view -H {input} | grep '^@RG' | sed \"s/.*SM:\\([^\\t]*\\).*/\\1/g\" | uniq`; \
        gatk --java-options -Xmx{params.mem}g Mutect2 \
            -R {config[reference_fasta]} \
            -I {input.bam} \
            -L {input.intervallist} \
            --tumor-sample $TUMOR_SM \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            -O {output} \
            > {log} 2>&1"

            #        TUMOR_SM=${TUMOR_SM:-\"TUMOR\"}; \

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## MergePON
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## This is an intermediate "gather" step. Because `callpon` is run in scatter mode, using
## and interval list and spread across 50 jobs for each normal sample, this step uses
## MergeVCF to merge all intermediate VCF files into one VCF per normal sample
## See: 
## https://gatkforums.broadinstitute.org/gatk/discussion/11136/how-to-call-somatic-mutations-using-gatk4-mutect2
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule mergepon:
    input:
        lambda wildcards: expand("results/callpon/{batch}/{aliquot_id}/{aliquot_id}.{interval}.pon.vcf", batch=wildcards.batch, aliquot_id=wildcards.aliquot_id, interval=WGS_SCATTERLIST)
    output:
        temp("results/mergepon/{batch}/{aliquot_id}.pon.vcf")
    params:
        mem = CLUSTER_META["mergepon"]["mem"]
    threads:
        CLUSTER_META["mergepon"]["ppn"]
    log:
        "logs/mergepon/{aliquot_id}.log"
    benchmark:
        "benchmarks/mergepon/{aliquot_id}.txt"
    message:
        "Merging VCF files (PON)\n"
        "Sample: {wildcards.aliquot_id}"
    run:
        input_cat = " ".join(["-I " + s for s in input])
        shell("gatk --java-options -Xmx{params.mem}g MergeVcfs \
            {input_cat} \
            -O {output} \
            > {log} 2>&1")

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## CreatePON
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Because PON calls on normal samples were done individually for each normal, this step
## merges calls from all samples from a given batch to build a confident PON
## Protected output, this is the final batch-specifc PON file
## See: 
## https://gatkforums.broadinstitute.org/gatk/discussion/11136/how-to-call-somatic-mutations-using-gatk4-mutect2
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule createpon:
    input:
        lambda wildcards: expand("results/mergepon/{batch}/{aliquot_id}.pon.vcf", batch=wildcards.batch, aliquot_id=BATCH_TO_NORMAL[wildcards.batch])
    output:
        protected("results/pon/{batch}.pon.vcf")
    params:
        dir = "results/qc/multiqc",
        mem = CLUSTER_META["createpon"]["mem"]
    threads:
        CLUSTER_META["createpon"]["ppn"]
    log:
        "logs/createpon/{batch}.createpon.log"
    benchmark:
        "benchmarks/createpon/{batch}.createpon.txt"
    message:
        "Creating panel of normals from multiple Mutect2 VCFs\n"
        "Batch: {wildcards.batch}"
    run:
        vcfs = " ".join(["--vcfs=" + s for s in input])
        shell("gatk --java-options -Xmx{params.mem}g CreateSomaticPanelOfNormals \
            {vcfs} \
            --duplicate-sample-strategy THROW_ERROR \
            --output {output} \
            > {log} 2>&1")

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Call SNV
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## This step uses Mutect2 to call variants on a tumor-normal pair, using a panel-of-normals
## and a germline resource as reference
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## GATK parameters taken from "GATK4_SomaticSNVindel_worksheet.pdf"
## Extract SM tag from BAM file
## Code snippet taken from https://github.com/IARCbioinfo/BAM-tricks
## https://software.broadinstitute.org/gatk/documentation/article.php?id=3891
## Calculate af-of-alleles: 1/(2*123136) = 0.00000406055 (6 sign. digits)
## 06/09/2018
## Removed -disable-read-filter MateOnSameContigOrNoMappedMateReadFilter
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule callsnv:
    input:
        tumor = lambda wildcards: "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["tumor_aliquot_id"]),
        normal = lambda wildcards: "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["normal_aliquot_id"]),
        pon = lambda wildcards: "results/pon/{batch}.pon.vcf".format(batch=PAIRS_DICT[wildcards.pair_id]["project_id"]),
        intervallist = lambda wildcards: "{dir}/{interval}/scattered.interval_list".format(dir=config["wgs_scatterdir"], interval=wildcards.interval)
    output:
        vcf = temp("results/m2vcf-scatter/{pair_id}.{interval}.vcf"),
        bam = "results/m2bam/{pair_id}.{interval}.bam"
    params:
        mem = CLUSTER_META["callsnv"]["mem"]
    threads:
        CLUSTER_META["callsnv"]["ppn"]
    log:
        "logs/callsnv/{pair_id}.log"
    benchmark:
        "benchmarks/callsnv/{pair_id}.txt"
    message:
        "Calling SNVs (Mutect2)\n"
        "Pair: {wildcards.pair_id}\n"
        "Interval: {wildcards.interval}"
    shell:
        "TEST_NAM=`samtools view -H {input.tumor} | grep '^@RG' | sed \"s/.*SM:\\([^\\t]*\\).*/\\1/g\" | uniq`;"
        "CTRL_NAM=`samtools view -H {input.normal} | grep '^@RG' | sed \"s/.*SM:\\([^\\t]*\\).*/\\1/g\" | uniq`;"
        #"TEST_NAM=${TEST_NAM:-\"TUMOR\"};"
        #"CTRL_NAM=${CTRL_NAM:-\"NORMAL\"};"
        "gatk --java-options -Xmx{params.mem}g Mutect2 \
            -R {config[reference_fasta]} \
            -I {input.tumor} \
            -I {input.normal} \
            -L {input.intervallist} \
            --tumor-sample $TEST_NAM \
            --normal-sample $CTRL_NAM \
            --panel-of-normals {input.pon} \
            --germline-resource {config[gnomad_vcf]} \
            --af-of-alleles-not-in-resource {config[af_of_alleles_not_in_resource]} \
            --dont-use-soft-clipped-bases true \
            --standard-min-confidence-threshold-for-calling 20 \
            -O {output.vcf} \
            -bamout {output.bam} \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Merge SNV
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## This is an intermediate "gather" step. Like `callpon`, `callsnv` is run in scatter mode, 
## using an interval list and spread across 50 jobs for each normal sample, this step uses
## MergeVCF to merge all intermediate VCF files into one VCF per normal sample
## See: 
## https://gatkforums.broadinstitute.org/gatk/discussion/11136/how-to-call-somatic-mutations-using-gatk4-mutect2
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule mergesnv:
    input:
        lambda wildcards: expand("results/m2vcf-scatter/{pair_id}.{interval}.vcf", pair_id=wildcards.pair_id, interval=WGS_SCATTERLIST)
    output:
        protected("results/m2vcf/{pair_id}.vcf")
    params:
        mem = CLUSTER_META["mergesnv"]["mem"]
    threads:
        CLUSTER_META["mergesnv"]["ppn"]
    log:
        "logs/mergesnv/{pair_id}.log"
    benchmark:
        "benchmarks/mergesnv/{pair_id}.txt"
    message:
        "Merging VCF files (M2)\n"
        "Pair: {wildcards.pair_id}"
    run:
        input_cat = " ".join(["-I " + s for s in input])
        shell("gatk --java-options -Xmx{params.mem}g MergeVcfs \
            {input_cat} \
            -O {output} \
            > {log} 2>&1")

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Summarize read support for known variant sites
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# @sbamin What is small_exac_common_3_b37.vcf.gz? Haven't read further but search only gives me this link: https://crc.pitt.edu/variantcalling
# @barthf This is used in broad pipelines

rule pileupsummaries:
    input:
        "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam"
    output:
        "results/pileupsummaries/{aliquot_id}.pileupsummaries.txt"
    params:
        mem = CLUSTER_META["pileupsummaries"]["mem"]
    threads:
        CLUSTER_META["pileupsummaries"]["ppn"]
    log:
        "logs/pileupsummaries/{aliquot_id}.log"
    benchmark:
        "benchmarks/pileupsummaries/{aliquot_id}.txt"
    message:
        "Generating pileupsummaries\n"
        "Sample: {wildcards.aliquot_id}"
    shell:
        "gatk --java-options -Xmx{params.mem}g GetPileupSummaries \
            -I {input} \
            -V {config[tiny_vcf]} \
            -O {output} \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Estimate contamination
## Input: pileup summaries table
## Output: contamination table
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# @sbamin We don't necessarily need contamination table based on variants across all chrs, and 2-3 large chrs are good enough.
# Idea is to query sufficient population level variants to estimate normal contamination of tumor sample.
# @barthf this is fast enough that I'm just gonna leave this for now

rule calculatecontamination:
    input:
        tumortable = lambda wildcards: "results/pileupsummaries/{aliquot_id}.pileupsummaries.txt".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["tumor_aliquot_id"]),
        normaltable = lambda wildcards: "results/pileupsummaries/{aliquot_id}.pileupsummaries.txt".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["normal_aliquot_id"]),
    output:
        "results/contamination/{pair_id}.contamination.txt"
    params:
        mem = CLUSTER_META["calculatecontamination"]["mem"]
    threads:
        CLUSTER_META["calculatecontamination"]["ppn"]
    log:
        "logs/calculatecontamination/{pair_id}.log"
    benchmark:
        "benchmarks/calculatecontamination/{pair_id}.txt"
    message:
        "Computing contamination\n"
        "Pair: {wildcards.pair_id}"
    shell:
        "gatk --java-options -Xmx{params.mem}g CalculateContamination \
            -I {input.tumortable} \
            --matched-normal {input.normaltable} \
            -O {output} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## If variants have not been filtered, filter, else done
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule filtermutect:
    input:
        vcf = "results/m2vcf/{pair_id}.vcf",
        tab = "results/contamination/{pair_id}.contamination.txt"
    output:
        temp("results/m2filter/{pair_id}.filtered.vcf")
    params:
        mem = CLUSTER_META["filtermutect"]["mem"]
    threads:
        CLUSTER_META["filtermutect"]["ppn"]
    log:
        "logs/filtermutect/{pair_id}.log"
    benchmark:
        "benchmarks/filtermutect/{pair_id}.txt"
    message:
        "Filtering Mutect2 calls\n"
        "Pair: {wildcards.pair_id}"
    shell:    
        "gatk --java-options -Xmx{params.mem}g FilterMutectCalls \
            -V {input.vcf} \
            --contamination-table {input.tab} \
            -O {output} \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Collect metrics on sequencing context artifacts
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule collectartifacts:
    input:
        "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam"
    output:
        "results/artifacts/{aliquot_id}.pre_adapter_detail_metrics.txt"
    params:
        prefix = "results/artifacts/{aliquot_id}",
        mem = CLUSTER_META["collectartifacts"]["mem"]
    threads:
        CLUSTER_META["collectartifacts"]["ppn"]
    log:
        "logs/collectartifacts/{aliquot_id}.log"
    benchmark:
        "benchmarks/collectartifacts/{aliquot_id}.txt"
    message:
        "Collecting sequencing artifact metrics\n"
        "Sample: {wildcards.aliquot_id}"
    shell:
        "gatk --java-options -Xmx{params.mem}g CollectSequencingArtifactMetrics \
            -I {input} \
            -O {params.prefix} \
            --FILE_EXTENSION \".txt\" \
            -R {config[reference_fasta]} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Filter by orientation bias
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule filterorientation:
    input:
        art = lambda wildcards: "results/artifacts/{aliquot_id}.pre_adapter_detail_metrics.txt".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["tumor_aliquot_id"]),
        vcf = "results/m2filter/{pair_id}.filtered.vcf"
    output:
        "results/m2filter/{pair_id}.filtered2.vcf"
    params:
        mem = CLUSTER_META["filterorientation"]["mem"]
    threads:
        CLUSTER_META["filterorientation"]["ppn"]
    log:
        "logs/filterorientation/{pair_id}.log"
    benchmark:
        "benchmarks/filterorientation/{pair_id}.txt"
    message:
        "Filtering Mutect2 calls by orientation bias\n"
        "Pair: {wildcards.pair_id}"
    shell:
        "gatk --java-options -Xmx{params.mem}g FilterByOrientationBias \
            -AM \"G/T\" \
            -AM \"C/T\" \
            -V {input.vcf} \
            -P {input.art} \
            -O {output} \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            > {log} 2>&1"

# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# ## VariantFiltration
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# ## Hard filter variants, removes "bad" variants from pre-VEP VCF
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# rule hardfilter:
#     input:
#         "results/m2filter/{pair_id}.filtered2.vcf"
#     output:
#         "results/m2filter/{pair_id}.filtered3.vcf"
#     params:
#         mem = CLUSTER_META["hardfilter"]["mem"]
#     threads:
#         CLUSTER_META["hardfilter"]["ppn"]
#     log:
#         "logs/hardfilter/{pair_id}.log"
#     benchmark:
#         "benchmarks/hardfilter/{pair_id}.txt"
#     message:
#         "Hard-filtering soft-filtered Mutect2 calls\n"
#         "Pair: {wildcards.pair_id}"
#     shell:
#         "gatk --java-options -Xmx{params.mem}g VariantFiltration \
#             -AM \"G/T\" \
#             -AM \"C/T\" \
#             -V {input.vcf} \
#             -P {input.art} \
#             -O {output} \
#             --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
#             > {log} 2>&1"

# @sbamin I suggest that we should put hold at VEP step until we finalize consensus calls. 
# Idea is too get how many filtered calls we have and with what level of confidence (based on consensus calls from other callers).
# To me, VEP with extended anntoations is of little value until vcfs are at freeze level.
# Also, we should have uniform processing past alignment step across both canine and human (adult + pediatric) life history projects.
# A few callers and filtering steps that we already have run on canine are at https://github.com/TheJacksonLaboratory/hourglass/issues/9.
# If we are in need of base annotations, like exonic vs promoter vs non-coding, that should be quick enough from VEP, 
# but I suggest to have a hold on extended annotations using vcf2maf as I guess we will end up doing those iteratively until we freeze vcfs.

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## USE VEP
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule vep:
    input:
        tumor = lambda wildcards: "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["tumor_aliquot_id"]),
        normal = lambda wildcards: "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["normal_aliquot_id"]),
        vcf = "results/m2filter/{pair_id}.filtered2.vcf"
    output:
        "results/vep/{pair_id}.filtered2.anno.maf"
    params:
        mem = CLUSTER_META["vep"]["mem"]
    threads:
        CLUSTER_META["vep"]["ppn"]
    log:
        "logs/vep/{pair_id}.log"
    benchmark:
        "benchmarks/vep/{pair_id}.txt"
    message:
        "Running VEP (variant annotation) on filtered Mutect2 calls\n"
        "Pair: {wildcards.pair_id}"
    shell:
        "TEST_NAM=`samtools view -H {input.tumor} | grep '^@RG' | sed \"s/.*SM:\\([^\\t]*\\).*/\\1/g\" | uniq`;"
        "CTRL_NAM=`samtools view -H {input.normal} | grep '^@RG' | sed \"s/.*SM:\\([^\\t]*\\).*/\\1/g\" | uniq`;"
        "{config[vcf2maf]} \
            --input-vcf {input.vcf} \
            --output-maf {output} \
            --vep-path {config[veppath]} \
            --vep-data {config[vepdata]} \
            --vep-forks 2 \
            --ref-fasta {config[reference_fasta]} \
            --filter-vcf {config[gnomad_vcf]} \
            --tumor-id $TEST_NAM \
            --normal-id $CTRL_NAM \
            --species homo_sapiens \
            --ncbi-build GRCh37 \
            > {log} 2>&1"

## END ##
