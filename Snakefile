## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Snakefile for GLASS-WG pipeline
## Authors: Floris Barthel, Samir Amin, Frederick Varn
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

import os
import pandas as pd
import itertools

## Import manifest processing functions
from python.glassfunc import dbconfig, locate
from python.PostgreSQLManifestHandler import PostgreSQLManifestHandler
from python.JSONManifestHandler import JSONManifestHandler

## Connect to database
## dbconf = dbconfig("/home/barthf/.odbc.ini", "VerhaakDB")
dbconf = dbconfig(config["db"]["configfile"], config["db"]["configsection"])

#print("Cohort set to ", str(config["cohort"]))
by_cohort = None
if len(str(config["cohort"])) > 0:
    by_cohort = str(config["cohort"]).zfill(2)

## Instantiate manifest
manifest = PostgreSQLManifestHandler(host = dbconf["servername"], port = dbconf["port"], user = dbconf["username"], password = dbconf["password"], database = dbconf["database"],
    source_file_basepath = config["data"]["source_path"], aligned_file_basepath = config["data"]["realn_path"], from_source = config["from_source"], by_cohort = by_cohort)
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

## We do not want the additional DAG processing if not from source
#if(config["from_source"]):
#    include: "snakemake/download.smk"
#include: "snakemake/align.smk"

# include: "snakemake/haplotype-map.smk"
include: "snakemake/fingerprinting.smk"
# include: "snakemake/telseq.smk"
include: "snakemake/mutect2.smk"
include: "snakemake/mutect2-post.smk"
# include: "snakemake/varscan2.smk"
# include: "snakemake/cnvnator.smk"
# include: "snakemake/lumpy.smk"
# include: "snakemake/delly.smk"
# include: "snakemake/manta.smk"
include: "snakemake/cnv.smk"
include: "snakemake/sequenza.smk"
include: "snakemake/optitype.smk"
include: "snakemake/pvacseq.smk"
#include: "snakemake/cnv-post.smk"
include: "snakemake/titan.smk"
#include: "snakemake/pyclone.smk"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Upload coverage to database
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule cov2db:
    input:
        metrics = lambda wildcards: expand("results/align/wgsmetrics/{aliquot_barcode}.WgsMetrics.txt", aliquot_barcode = manifest.getSelectedAliquots())
    output:
        tsv = "results/align/wgsmetrics.merged.tsv"
    params:
        mem = CLUSTER_META["cov2db"]["mem"]
    threads:
        CLUSTER_META["cov2db"]["ppn"]
    #conda:
    #    "../envs/r.yaml"
    log:
        "logs/align/cov2db/cov2db.log"
    benchmark:
        "benchmarks/align/cov2db/cov2db.txt"
    message:
        "Merge coverage file and convert to TSV for easy database upload"
    script:
        "R/snakemake/cov2db.R"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Haplotype map creation rule
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

#rule build_haplotype_map:
#    input:
#        "data/ref/fingerprint.filtered.map"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Alignment rule
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule align:
    input:
        expand("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam", aliquot_barcode = manifest.getSelectedAliquots()),
        expand("results/align/wgsmetrics/{aliquot_barcode}.WgsMetrics.txt", aliquot_barcode = manifest.getSelectedAliquots()),
        expand("results/align/validatebam/{aliquot_barcode}.ValidateSamFile.txt", aliquot_barcode = manifest.getSelectedAliquots()),
        lambda wildcards: ["results/align/fastqc/{sample}/{sample}.{rg}.unaligned_fastqc.html".format(sample = aliquot_barcode, rg = readgroup)
          for aliquot_barcode, readgroups in manifest.getSelectedReadgroupsByAliquot().items()
          for readgroup in readgroups]

rule gencode:
    input:
        expand("results/align/gencode-coverage/{aliquot_barcode}.gencode-coverage.tsv", aliquot_barcode = [f for f in manifest.getSelectedAliquots() if manifest.isExome(f)])

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Get genic coverage given a BAM file
## URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule gencode_coverage:
    input:
        "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam"
    output:
        "results/align/gencode-coverage/{aliquot_barcode}.gencode-coverage.tsv"
    params:
        mem = CLUSTER_META["gencode_coverage"]["mem"]
    #conda:
    #    "envs/align.yaml"
    threads:
        CLUSTER_META["gencode_coverage"]["ppn"]
    log:
        "logs/align/gencode-coverage/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/align/gencode-coverage/{aliquot_barcode}.txt"
    message:
        "Computing coverage using flattened gencode GTF\n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
        "set +o pipefail; /opt/software/helix/samtools/1.8/bin/samtools view -q 10 -b {input} | \
            /opt/software/helix/BEDtools/2.27.0/bin/bedtools coverage -a {config[gencode_gtf_flat]} -b stdin -d -sorted -g {config[bedtools_genome]} | \
            /opt/software/helix/BEDtools/2.27.0/bin/bedtools groupby -i stdin -g 1,2,3,4,5 -c 7 -o sum | \
            sort -k5,5 | \
            /opt/software/helix/BEDtools/2.27.0/bin/bedtools groupby -i stdin -g 5 -c 4,6 -o sum,sum | \
            awk -F\"[+\\t]\" 'BEGIN {{OFS=\"\\t\"}}{{for(i=1;i<(NF-1);i++){{split($i,g,\".\"); print g[1],$(NF-1),$NF}}}}' \
            > {output} 2> {log}"

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
## SNV rule (Mutect2)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule mutect2:
    input:
        expand("results/mutect2/m2filter/{case_barcode}.filtered.vcf.gz", case_barcode = manifest.getSelectedCases()),
        expand("results/mutect2/ssm2filter/{pair_barcode}.filtered.vcf.gz", pair_barcode = manifest.getSelectedPairs())

rule ssmutect2:
    input:
        expand("results/mutect2/ssm2filter/{pair_barcode}.filtered.vcf.gz", pair_barcode = manifest.getSelectedPairs())

rule m2db:
    input:
        #"results/mutect2/consensusvcf/consensus.normalized.sorted.funcotated.tsv",
        expand("results/mutect2/geno2db/{case_barcode}.info.tsv", case_barcode = manifest.getSelectedCases()),
        expand("results/mutect2/geno2db/{case_barcode}.geno.tsv", case_barcode = manifest.getSelectedCases())

rule sequenza:
    input:
        expand("results/sequenza/mergeseqz/{pair_barcode}.small.seqz.gz", pair_barcode = manifest.getSelectedPairs()),
        expand("results/sequenza/seqzR/{pair_barcode}/{pair_barcode}_cellularity.ploidy.txt", pair_barcode = manifest.getSelectedPairs())

rule titancna:
    input:
        expand("results/cnv/titanfinal/seg/{pair_barcode}.seg.txt", pair_barcode = manifest.getSelectedPairs())

# rule mutect2post:
# 	input:
# 		expand("results/mutect2/m2post/{pair_barcode}.normalized.sorted.vcf.gz", pair_barcode = manifest.getSelectedPairs())

# rule genotypefreebayes:
#     input:
#         expand("results/mutect2/freebayes/{aliquot_barcode}.normalized.sorted.vcf.gz", aliquot_barcode = manifest.getSelectedAliquots())

# rule massfreebayes:
#     input:
#         expand("results/mutect2/freebayes/batch{batch}/{aliquot_barcode}.normalized.sorted.vcf.gz", batch = [str(i) for i in range(2,6)], aliquot_barcode = manifest.getSelectedAliquots())

# rule genotypevcf2vcf:
#     input:
#         expand("results/mutect2/genotypes/{aliquot_barcode}.normalized.sorted.vcf.gz", aliquot_barcode = manifest.getSelectedAliquots())

# rule finalfreebayes:
#     input:
#         expand("results/mutect2/batches2db/{aliquot_barcode}.normalized.sorted.tsv", aliquot_barcode = manifest.getSelectedAliquots())

# rule preparem2pon:
#     input:
#         expand("results/mutect2/mergepon/{aliquot_barcode}.pon.vcf", aliquot_barcode = manifest.getPONAliquots())

# rule genodb:
#     input:
#         expand("results/mutect2/geno2db/{aliquot_barcode}.normalized.sorted.tsv", aliquot_barcode = manifest.getSelectedAliquots())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## PON rule (Mutect2)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# rule mutect2pon:
#     input:
#     	"results/mutect2/pon/pon.vcf"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## HLAtyping rule (OptiType)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule call_hla:
	input:
		expand("results/optitype/HLA_calls/{aliquot_barcode}/{aliquot_barcode}_result.tsv", aliquot_barcode = manifest.getSelectedAliquots())
		
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Neoantigen rule (pVACseq)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule call_neoag:
	input:
		expand("results/pvacseq/neoantigens/{case_barcode}/MHC_Class_I/{case_barcode}.final.tsv", case_barcode = manifest.getSelectedCases())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## SNV rule (VarScan2)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule varscan2:
    input:
        expand("results/varscan2/fpfilter/{pair_barcode}.{type}.Somatic.hc.final.vcf", pair_barcode = manifest.getSelectedPairs(), type = ["snp", "indel"])

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## CNV calling pipeline
## Run snakemake with target 'svprepare'
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 



rule cnv:
    input:
        expand("results/cnv/plots/{aliquot_barcode}.pdf", aliquot_barcode = manifest.getSelectedAliquots()),
        expand("results/cnv/callsegments/{aliquot_barcode}.called.seg", aliquot_barcode = manifest.getSelectedAliquots())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Estimate TL using telseq
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule telseq:
    input:
        expand("results/telseq/{aliquot_barcode}.telseq.txt", aliquot_barcode = manifest.getSelectedAliquots())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run PyClone
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule pyclone:
    input:
        lambda wildcards: expand("results/pyclone/run/{pyclone_short_name}/plots/loci/{plot_type}.pdf", pyclone_short_name = manifest.getPyCloneShortNames(), plot_type = ['density','parallel_coordinates','scatter','similarity_matrix','vaf_parallel_coordinates','vaf_scatter']),
        lambda wildcards: expand("results/pyclone/run/{pyclone_short_name}/plots/clusters/{plot_type}.pdf", pyclone_short_name = manifest.getPyCloneShortNames(), plot_type = ['density','parallel_coordinates','scatter']),
        lambda wildcards: expand("results/pyclone/run/{pyclone_short_name}/tables/{table_type}.tsv", pyclone_short_name = manifest.getPyCloneShortNames(), table_type = ['cluster','loci'])

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Fingerprinting pipeline
## Check sample and case fingerprints
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule fingerprint:
   input:
       expand("results/fingerprinting/sample/{aliquot_barcode}.crosscheck_metrics", aliquot_barcode = manifest.getSelectedAliquots()),
       expand("results/fingerprinting/case/{case_barcode}.crosscheck_metrics", case_barcode = manifest.getSelectedCases())
       #"results/fingerprinting/GLASS.crosscheck_metrics",
       

## END ##
