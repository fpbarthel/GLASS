#######################################################
# Create manifest for Yung MD Anderson sample (n = 1, GLASS)
# Date: 2018.08.08
# Author: Kevin J.
#######################################################
# Local directory for github repo.
mybasedir = "/Users/johnsk/Documents/Life-History/GLASS-WG"
setwd(mybasedir)

# Files with information about "Yung" fastq information and barcodes.
# Actual fastq files are stored: /fastscratch/johnsk/GLASS-WG/mdacc/ or on tier2 (long-term).
MDA_batch2_file_path = "data/sequencing-information/MDACC/file_list_C202SC18030593_batch2_n94_20180603.tsv"

# Clinical dataset priovided by Kristin Alfaro-Munoz at MD Anderson.
# The "Yung" sample was included in this dataset.
mda_master_path = "data/clinical-data/MDACC/MDA-Clinical-Dataset/Master Log for WGS_sampletype.20180630.xlsx"

# 2018.07.06 Katie Shao (Novogene) provided the following sheet linking libraryID with submitted sample names in the email:
# "Re:RE: Novogene Project Report - Confirm -C202SC18030593-Davis-MDACC-156-libseq-hWGS-WOBI-NVUS2018022505[J1-C202SC18030593-1000]"
novogene_sample_path = "data/sequencing-information/MDACC/Novogene_SIF_14.xlsx"

# Completed life-history barcodes.
life_history_barcodes = "data/sequencing-information/master_life_history_uniform_naming_complete.txt"

# Create extensions for samples.
json_ext = "json"
text_ext = "tsv"

# Create output files for each metadata set.
cases_file      = "data/manifest/mdacc/cases"
samples_file    = "data/manifest/mdacc/samples"
aliquots_file   = "data/manifest/mdacc/aliquots"
readgroups_file = "data/manifest/mdacc/readgroups"
files_file      = "data/manifest/mdacc/files"
pairs_file      = "data/manifest/mdacc/pairs"

#######################################################

library(tidyverse)
library(openxlsx)
library(rjson)
library(jsonlite)
library(listviewer)
library(stringi)
library(stringr)

#######################################################
# We need to generate the following fields required by the SNV snakemake pipeline:
# aliquots, files, cases, samples, pairs, and readgroups.

### MDA (Novogene cohort) barcode (re)generation ####
# Novogene has processed 120 samples. 
mda_batch2_df <- read.table(MDA_batch2_file_path, col.names="file_path", stringsAsFactors = F)
mda_batch2_df$filename = sapply(strsplit(mda_batch2_df$file_path, "/"), "[[", 3)

# Replace the placeholder for pwd "./" from bash cmd: "find -name ".fq.gz" in parent directory of fastqs. 
mda_batch1_df$file_path <- gsub("\\./RAW/20180405/128.120.88.242/C202SC18030593/raw_data/", "/fastscratch/johnsk/GLASS-WG/mdacc/C202SC18030593_batch1_n14_20180405/", mda_batch1_df$file_path)
mda_batch2_df$file_path <- gsub("^", "/fastscratch/johnsk/GLASS-WG/mdacc/", mda_batch2_df$file_path)
mda_batch3_df$file_path <- gsub("^", "/fastscratch/johnsk/GLASS-WG/mdacc/C202SC18030593_batch2_n13_20180716", mda_batch3_df$file_path)

# Create an old sample.id for these subjects to be linked. No longer *necessary*, but kept in case it's helpful.
mda_batch1_df$old_sample_id <- gsub( "_USPD.*$", "", mda_batch1_df$filename)
mda_batch2_df$old_sample_id <- gsub( "_USPD.*$", "", mda_batch2_df$filename)
mda_batch3_df$old_sample_id <- gsub( "_USPD.*$", "", mda_batch3_df$filename)

# Combine these three sequencing data sets. There should be 552 files.
mda_df <- bind_rows(mda_batch1_df, mda_batch2_df, mda_batch3_df)

# Katie Shao (Novogene provided).
novogene_linker = readWorkbook(novogene_sample_path, sheet = 1, startRow = 18, colNames = TRUE)
novogene_linker_unique <- novogene_linker[1:121, ] # 121st sample represents non-GLASS sample.



