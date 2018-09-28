## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Snakefile for GLASS-WG pipeline
## Authors: Floris Barthel, Samir Amin
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

import os
import pandas as pd
import itertools

## Import manifest processing functions
from python.glassfunc import dbconfig, locate
from python.PostgreSQLManifestHandler import PostgreSQLManifestHandler
from python.JSONManifestHandler import JSONManifestHandler

## Connect to database
dbconf = dbconfig(config["db"]["configfile"], config["db"]["configsection"])

## Instantiate manifest
manifest = PostgreSQLManifestHandler(host = dbconf["servername"], port = dbconf["port"], user = dbconf["username"], password = dbconf["password"], database = dbconf["database"])
print(manifest)

## Set working directory based on configuration file
workdir: config["workdir"]

## GDC token file for authentication
KEYFILE = config["gdc_token"]

## Cluster metadata (memory, CPU, etc)
CLUSTER_META = json.load(open(config["cluster_json"]))

## List of scatterlist items to iterate over
## Each Mutect2 run spawns 50 jobs based on this scatterlist

WGS_SCATTERLIST = ["temp_{num}_of_50".format(num=str(j+1).zfill(4)) for j in range(50)]

## Load modules
include: "snakemake/haplotype-map.smk"
include: "snakemake/download.smk"
include: "snakemake/align.smk"
include: "snakemake/fingerprinting.smk"
include: "snakemake/telseq.smk"
include: "snakemake/mutect2.smk"
# include: "snakemake/lumpy.smk"
# include: "snakemake/cnv-gatk.smk"
# include: "snakemake/varscan2.smk"
# include: "snakemake/delly.smk"
# include: "snakemake/manta.smk"
# include: "snakemake/cnvnator.smk"
# include: "snakemake/vep.smk"

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
        expand("results/align/bqsr/{aliquot_id}.realn.mdup.bqsr.bam", aliquot_id = manifest.getSelectedAliquots()),
        expand("results/align/wgsmetrics/{aliquot_id}.WgsMetrics.txt", aliquot_id = manifest.getSelectedAliquots()),
        expand("results/align/validatebam/{aliquot_id}.ValidateSamFile.txt", aliquot_id = manifest.getSelectedAliquots()),
        lambda wildcards: ["results/align/fastqc/{sample}/{sample}.{rg}.unaligned_fastqc.html".format(sample = aliquot_barcode, rg = readgroup)
          for aliquot_barcode, readgroups in manifest.getSelectedReadgroupsByAliquot().items()
          for readgroup in readgroups]

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Download only rule
## Run snakemake with 'snakemake download_only' to activate
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

#rule download_only:
#   input: expand("{file}", file = ALIQUOT_TO_BAM_PATH.values())

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
    input: expand("results/mutect2/final/{pair_id}.final.vcf", pair_id = manifest.getSelectedPairs())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## SNV rule (VarScan2)
## Run snakemake with target 'snv'
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule vs2:
    input:
        expand("results/varscan2/final/{pair_id}.somatic.hc.filtered.final.vcf.gz", pair_id = manifest.getSelectedPairs())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## CNV calling pipeline
## Run snakemake with target 'svprepare'
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule cnv:
    input:
        expand("results/cnv/callsegments/{pair_id}.called.seg", pair_id = manifest.getSelectedPairs()),
        expand("results/cnv/plotmodeledsegments/{pair_id}/{pair_id}.modeled.png", pair_id = manifest.getSelectedPairs()),
        expand("results/cnv/plotcr/{aliquot_id}/{aliquot_id}.denoised.png", aliquot_id = manifest.getSelectedAliquots())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Call SV using Delly
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule delly:
    input:
        expand("results/delly/filter/{pair_id}.prefilt.bcf", pair_id = manifest.getSelectedPairs())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Call SV using LUMPY-SV
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule lumpy:
    input:
        expand("results/lumpy/svtyper/{pair_id}.dict.svtyper.vcf", pair_id = manifest.getSelectedPairs()),
        expand("results/lumpy/libstat/{pair_id}.libstat.pdf", pair_id = manifest.getSelectedPairs()),
        expand("results/lumpy/filter/{pair_id}.dict.svtyper.filtered.vcf", pair_id = manifest.getSelectedPairs())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Call SV using Manta
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule manta:
    input:
        expand("results/manta/{pair_id}/results/variants/somaticSV.vcf.gz", pair_id = manifest.getSelectedPairs())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Call CNV using Manta
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule cnvnator:
    input:
        expand("results/cnvnator/vcf/{aliquot_id}.call.vcf", aliquot_id = manifest.getSelectedAliquots())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Call SVs using all callers
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule svdetect:
    input:
        expand("results/lumpy/svtyper/{pair_id}.dict.svtyper.vcf", pair_id = manifest.getSelectedPairs()),
        expand("results/lumpy/libstat/{pair_id}.libstat.pdf", pair_id = manifest.getSelectedPairs()),
        expand("results/delly/call/{pair_id}.vcf.gz", pair_id = manifest.getSelectedPairs()),
        expand("results/manta/{pair_id}/results/variants/somaticSV.vcf.gz", pair_id = manifest.getSelectedPairs())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Estimate TL using telseq
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule telseq:
    input:
        expand("results/telseq/{aliquot_id}.telseq.txt", aliquot_id = manifest.getSelectedAliquots())

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
