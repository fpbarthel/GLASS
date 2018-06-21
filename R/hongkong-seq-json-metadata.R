#######################################################
# Generate a metadata json file for hong kong sequencing samples
# Date: 2018.06.01
# Author: Kevin J
#######################################################

setwd("/Users/johnsk/Documents/Life-History/GLASS-WG/")

HK_file_path = "data/sequencing-information/HK/tmp-hk-path-list.txt"
HK_meta_path = "data/sequencing-information/HK/hong-kong-sample-maps.txt"
life_history_barcodes = "data/sequencing-information/master_life_history_uniform_naming_incomplete.txt"

#######################################################

library(tidyverse)
library(openxlsx)
library(rjson)
library(jsonlite)
library(listviewer)

#######################################################
# To get the filenames from subsetted fastq in /fastscratch/:
# Generate file map for Hong Kong samples.
# ls *.fq > /fastscratch/johnsk/test/HongKong/HK-all-samples-list.txt

# Append path to samples.
# sed 's,^,'/fastscratch/johnsk/test/HongKong/',' HK-all-samples-list.txt > HK-Path-List.txt

# Move to local computer and process in R.
# rsync -avz johnsk@helix:/fastscratch/johnsk/test/HongKong/HK-Path-List.txt /Users/johnsk/Documents/Life-History/

# Each file needs the following information:
# case_id (replaces patient), case_project (replaces cohort), pairs (grouping normal, primary, and recurrence).
# sample: sample_type, sample_id, sample_type_code.
# file: file_name, file_format, file_uuid.
# readgroups: rg_ID, rg_PL, rg_PU, rg_LB, rg_DT, rg_SM, rg_CN.

# Load in samples from the Hong Kong data set.
hk_df = read.table(HK_file_path, col.names="file_name", stringsAsFactors = F)

# Retrieve the sample name to map to study center provided covariate sheet.
hk_df$sample_id = sub(".*/WG_ *(.*?) *_USPD.*", "\\1", hk_df$file_name)

# Combine with study center provided covariate sheet.
HongKong_covars =  read.table(HK_meta_path, header = T, stringsAsFactors = F)
HongKong_covars$Verhaak_Sample_ID = gsub("-","_", HongKong_covars$Verhaak_Sample_ID)
hk_map_df = HongKong_covars %>% left_join(hk_df, by = c("Verhaak_Sample_ID" = "sample_id"))

# Need to incorporate the life history specific barcodes.
barcodes = read_tsv(life_history_barcodes)
barcodes$Original_ID = gsub("-","_", barcodes$Original_ID)
barcodes_reformatted = barcodes %>% 
  select(uuid, Original_ID, Barcode)
hk_map_barcode_df = hk_map_df %>% left_join(barcodes_reformatted, by = c("Verhaak_Sample_ID" = "Original_ID"))

# Generate some read group information.
hk_rg_df = hk_map_barcode_df %>% 
  mutate(case_id = substr(Barcode, 1, 12)) %>% 
  mutate(case_project = "hong_kong") %>% 
  rename(sample_id = Barcode) %>% 
  mutate(rg_PL = "ILLUMINA") %>% 
  mutate(file_format = "FQ") %>% 
  mutate(rg_DT = strftime(as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%dT%H:%M:%S%z")) %>% 
  mutate(rg_PU = paste(substr(file_name, nchar(file_name)-19, nchar(file_name)-12), 
                       substr(file_name, nchar(file_name)-9, nchar(file_name)-9), sep=".")) %>% 
  mutate(rg_LB = sub(".*_ *(.*?) *_H.*", "\\1", file_name)) %>% 
  mutate(sample_type = recode_factor(Sample_Type, "Primary" = "Primary Tumor", "Blood" = "Blood Derived Normal", "Recurrence" = "First Recurrence Tumor")) %>% 
  mutate(sample_type_code = recode_factor(sample_type, "Primary Tumor" = "TP", "Blood Derived Normal" = "NB", "First Recurrence Tumor" = "R1")) %>% 
  mutate(rg_SM = sample_id) %>% 
  mutate(rg_CN = "NVGN_HK") %>% 
  mutate(rg_ID = paste0(substring(rg_PU, 1, 4), substring(rg_PU, nchar(rg_PU)-1, nchar(rg_PU)), "")) %>% 
  select(-one_of("GT_Sample_ID", "Patient_ID", "Cohort", "Sex", "Age", "Sample_Type", "Seq_Type", "SeqMachine", "WGS_ReadLength", "Verhaak_Sample_ID", "uuid"))

# Nest samples, files, and rg to create JSON.
nested_filtered_files_rgs = hk_rg_df %>%
  select(starts_with("case"), starts_with("sample"), starts_with("file"), starts_with("rg")) %>%
  nest(-starts_with("case"), -starts_with("sample"), -starts_with("file"), .key=readgroups) %>% # starts_with("rg")
  nest(-starts_with("case"), -starts_with("sample"), .key=files) %>%
  nest(-starts_with("case"), .key=samples)

jsonlite::toJSON(nested_filtered_files_rgs, pretty = T)
write(jsonlite::toJSON(nested_filtered_files_rgs, pretty = T), file = "data/sequencing-information/HK/hk-nested-rgs-test.json")

# Need to do the following: 
# 1) FASTQ 1 and 2 on one line with comma seperated
# 2) sample_id and case_id should reflect new sample/case ID, these reflect names going forward
# 3) add a pairs directive at case level that maps all pair of tumor/normal of interest for that case
# 4) consider havin rg_CN and case_project be identical (not sure on this but a thought)
# P.S. rg_SM should also match new sample ID

