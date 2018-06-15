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

## Cluster metadata (memory, CPU, etc)
CLUSTER_META    = json.load(open(config["cluster_json"]))

## JSON data
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
## @barthf : TO-DO

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

## Load modules
include: "snakemake/download.smk"
include: "snakemake/align.smk"
include: "snakemake/mutect2.smk"
include: "snakemake/vep.smk"

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

## END ##
