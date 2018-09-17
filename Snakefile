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

## Set working directory based on configuration file
workdir: config["workdir"]

## GDC token file for authentication
KEYFILE     = config["gdc_token"]

## Cluster metadata (memory, CPU, etc)
CLUSTER_META    = json.load(open(config["cluster_config"]))

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

## Turn an unnamed list of dicts into a nammed list of dicts
## Taken from stackoverflow
## https://stackoverflow.com/questions/4391697/find-the-index-of-a-dict-within-a-list-by-matching-the-dicts-value
def build_dict(seq, key):
    return dict((d[key], dict(d, index=index)) for (index, d) in enumerate(seq))

## CASES should be unique
## CHECK THAT ALL CASE_ID VALUES IN CASES ARE UNIQUE
## IN PROGRESS


## FILES -> FILE_UUID should be unique
## CHECK THAT ALL FILE_UUID VALUES IN FILES ARE UNIQUE
## IN PROGRESS


## PAIR -> PAIR_ID should be unique
## CHECK THAT ALL PAIR_ID VALUES IN PAIR ARE UNIQUE
## IN PROGRESS


## CASES -> DICT
CASES_DICT = build_dict(CASES, "case_id")

## SAMPLES -> DICT
SAMPLES_DICT = build_dict(SAMPLES, "sample_id")

## ALIQUOTS -> DICT
ALIQUOTS_DICT = build_dict(ALIQUOTS, "aliquot_id")

## FILES -> DICT
FILES_DICT = build_dict(FILES, "file_uuid")

## Pair IDs are unique, PAIRS -> DICT
PAIRS_DICT = build_dict(PAIRS, "pair_id")

## Aliquot IDs and BAM files map 1:1
ALIQUOT_TO_BAM_PATH = {}
for file in FILES:
    if file["file_format"] == "BAM":
        ALIQUOT_TO_BAM_PATH[ file["aliquot_id"] ] = file["file_path"]
        
## Case to aliquots
## Dict of aliquots per case
CASE_TO_ALIQUOT = {}
for aliquot in ALIQUOTS:
    aliquot["case_id"] = SAMPLES_DICT[ aliquot["sample_id"] ]["case_id"]  
    if aliquot["case_id"] not in CASE_TO_ALIQUOT:
        CASE_TO_ALIQUOT[ aliquot["case_id"] ] = [ aliquot["aliquot_id"] ]
    elif aliquot["aliquot_id"] not in CASE_TO_ALIQUOT[ aliquot["case_id"] ]:
        CASE_TO_ALIQUOT[ aliquot["case_id"] ].append(aliquot["aliquot_id"])

## Aliquots and RGIDs map 1:many
ALIQUOT_TO_RGID = {}     
ALIQUOT_TO_LEGACY_RGID = {}
for readgroup in READGROUPS:
    if readgroup["aliquot_id"] not in ALIQUOT_TO_RGID:
        ALIQUOT_TO_RGID[ readgroup["aliquot_id"] ] = [ readgroup["readgroup_id"] ]
    else:
        ALIQUOT_TO_RGID[ readgroup["aliquot_id"] ].append(readgroup["readgroup_id"])
    if "legacy_readgroup_id" not in readgroup or len(readgroup["legacy_readgroup_id"]) == 0:
        continue
    if readgroup["aliquot_id"] not in ALIQUOT_TO_LEGACY_RGID:
        ALIQUOT_TO_LEGACY_RGID[ readgroup["aliquot_id"] ] = [ readgroup["legacy_readgroup_id"] ]
    else:
        ALIQUOT_TO_LEGACY_RGID[ readgroup["aliquot_id"] ].append(readgroup["legacy_readgroup_id"])
        
## Readgroup information and 
## Aliquots and RGIDs map 1:many
## RGIDs are unique within an aliquot
## Aliquot IDs and fastQ files map 1:many
## Because FQ files are also seperated by readgroup, create dictionary of FQ files here as well
ALIQUOT_TO_READGROUP = {} 
ALIQUOT_TO_FQ_PATH = {}
for readgroup in READGROUPS:
    if readgroup["aliquot_id"] not in ALIQUOT_TO_READGROUP:
        ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ] = { readgroup["readgroup_id"] : readgroup }
    else:
        ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ][ readgroup["readgroup_id"] ] = readgroup
    ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ][ readgroup["readgroup_id"] ]["file_path"] = FILES_DICT[ ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ][ readgroup["readgroup_id"] ]["file_uuid"] ]["file_path"]
    ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ][ readgroup["readgroup_id"] ]["file_format"] = FILES_DICT[ ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ][ readgroup["readgroup_id"] ]["file_uuid"] ]["file_format"]
    if ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ][ readgroup["readgroup_id"] ]["file_format"] == "FQ":
        if readgroup["aliquot_id"] not in ALIQUOT_TO_FQ_PATH:
            ALIQUOT_TO_FQ_PATH[ readgroup["aliquot_id"] ] = {}
        ALIQUOT_TO_FQ_PATH[ readgroup["aliquot_id"] ][ readgroup["readgroup_id"] ] = ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ][ readgroup["readgroup_id"] ]["file_path"].split(",")

## List of scatterlist items to iterate over
## Each Mutect2 run spawns 50 jobs based on this scatterlist

WGS_SCATTERLIST = ["temp_{num}_of_50".format(num=str(j+1).zfill(4)) for j in range(50)]

## Load modules
include: "snakemake/haplotype-map.smk"
include: "snakemake/download.smk"
include: "snakemake/align.smk"
include: "snakemake/mutect2.smk"
include: "snakemake/vep.smk"
include: "snakemake/lumpy.smk"
include: "snakemake/cnv-gatk.smk"
include: "snakemake/varscan2.smk"
include: "snakemake/fingerprinting.smk"
include: "snakemake/delly.smk"
include: "snakemake/manta.smk"
include: "snakemake/cnvnator.smk"
include: "snakemake/telseq.smk"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Haplotype map creation rule
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule build_haplotype_map:
    input:
        "data/ref/fingerprint.filtered.map"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Alignment rule
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule align:
    input:
        expand("results/align/bqsr/{aliquot_id}.realn.mdup.bqsr.bam", aliquot_id=ALIQUOTS_DICT.keys()),
        expand("results/align/wgsmetrics/{aliquot_id}.WgsMetrics.txt", aliquot_id=ALIQUOT_TO_READGROUP.keys()),
        expand("results/align/validatebam/{aliquot_id}.ValidateSamFile.txt", aliquot_id=ALIQUOT_TO_READGROUP.keys()),
        lambda wildcards: ["results/align/fastqc/{sample}/{sample}.{rg}.unaligned_fastqc.html".format(sample=sample, rg=readgroup)
          for sample, readgroups in ALIQUOT_TO_RGID.items()
          for readgroup in readgroups]

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Download only rule
## Run snakemake with 'snakemake download_only' to activate
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule download_only:
   input: expand("{file}", file=ALIQUOT_TO_BAM_PATH.values())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## QC rule
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule qc:
    input: 
        "results/align/multiqc/multiqc_report.html"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## SNV rule
## Run snakemake with target 'snv'
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule mt2:
    input: expand("results/mutect2/m2filter/{pair_id}.filtered2.vcf", pair_id=PAIRS_DICT.keys())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## SNV rule (VarScan2)
## Run snakemake with target 'snv'
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule vs2:
    input:
        expand("results/varscan2/final/{pair_id}.somatic.hc.filtered.final.vcf.gz", pair_id=PAIRS_DICT.keys())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## CNV calling pipeline
## Run snakemake with target 'svprepare'
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule cnv:
    input:
        expand("results/cnv/callsegments/{pair_id}.called.seg", pair_id=PAIRS_DICT.keys()),
        expand("results/cnv/plotmodeledsegments/{pair_id}/{pair_id}.modeled.png", pair_id=PAIRS_DICT.keys()),
        expand("results/cnv/plotcr/{aliquot_id}/{aliquot_id}.denoised.png", aliquot_id=ALIQUOT_TO_READGROUP.keys())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Call SV using Delly
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule delly:
    input:
        expand("results/delly/filter/{pair_id}.prefilt.bcf", pair_id=PAIRS_DICT.keys())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Call SV using LUMPY-SV
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule lumpy:
    input:
        expand("results/lumpy/svtyper/{pair_id}.dict.svtyper.vcf", pair_id=PAIRS_DICT.keys()),
        expand("results/lumpy/libstat/{pair_id}.libstat.pdf", pair_id=PAIRS_DICT.keys()),
        expand("results/lumpy/filter/{pair_id}.dict.svtyper.filtered.vcf", pair_id=PAIRS_DICT.keys())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Call SV using Manta
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule manta:
    input:
        expand("results/manta/{pair_id}/results/variants/somaticSV.vcf.gz", pair_id=PAIRS_DICT.keys())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Call CNV using Manta
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule cnvnator:
    input:
        expand("results/cnvnator/vcf/{aliquot_id}.call.vcf", aliquot_id=ALIQUOT_TO_READGROUP.keys())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Call SVs using all callers
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule svdetect:
    input:
        expand("results/lumpy/svtyper/{pair_id}.dict.svtyper.vcf", pair_id=PAIRS_DICT.keys()),
        expand("results/lumpy/libstat/{pair_id}.libstat.pdf", pair_id=PAIRS_DICT.keys()),
        expand("results/delly/call/{pair_id}.vcf.gz", pair_id=PAIRS_DICT.keys()),
        expand("results/manta/{pair_id}/results/variants/somaticSV.vcf.gz", pair_id=PAIRS_DICT.keys())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Estimate TL using telseq
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule telseq:
    input:
        expand("results/telseq/{aliquot_id}.telseq.txt", aliquot_id=ALIQUOT_TO_READGROUP.keys())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Fingerprinting pipeline
## Check sample and case fingerprints
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

#rule fingerprint:
#    input:
#        expand("results/fingerprinting/sample/{aliquot_id}.crosscheck_metrics", aliquot_id=ALIQUOT_TO_READGROUP.keys()),
#        expand("results/fingerprinting/case/{case_id}.crosscheck_metrics", case_id=CASES_DICT.keys()),
#        expand("results/fingerprinting/batch/{batch}.crosscheck_metrics", batch=BATCH_TO_ALIQUOT.keys()),
#        "results/fingerprinting/GLASS-WG.crosscheck_metrics",
#        expand("results/fingerprinting/batch/{batch}.clustered.crosscheck_metrics", batch=BATCH_TO_ALIQUOT.keys()),
#        "results/fingerprinting/GLASS-WG.clustered.crosscheck_metrics"
       

## END ##
