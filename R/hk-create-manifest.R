#######################################################
# Create manifest for hong kong samples (GLASS)
# Date: 2018.06.20
# Author: Kevin J.
#######################################################

# Directory for github repo.
mybasedir = "/Users/johnsk/Documents/Life-History/GLASS-WG"
setwd(mybasedir)

HK_file_path = "data/sequencing-information/HK/hk_all_files_test.txt"
HK_meta_path = "data/sequencing-information/HK/hong-kong-sample-maps.txt"
life_history_barcodes = "data/sequencing-information/master_life_history_uniform_naming_incomplete.txt"

# Create extensions for samples.
json_ext = "json"
text_ext = "tsv"

# Create output files for each metadata set.
cases_file      = "data/manifest/hongkong/cases"
samples_file    = "data/manifest/hongkong/samples"
aliquots_file   = "data/manifest/hongkong/aliquots"
readgroups_file = "data/manifest/hongkong/readgroups"
files_file      = "data/manifest/hongkong/files"
pairs_file      = "data/manifest/hongkong/pairs"

#######################################################

library(tidyverse)
library(openxlsx)
library(rjson)
library(jsonlite)
library(listviewer)
library(stringi)
library(stringr)

#######################################################
# We need to generate the following fields required by the SNV snakemake:
# aliquots, cases, files, pairs, readgroups, and samples.

# Generate aliquot tsv.
master_sheet = read.delim("data/sequencing-information/master_life_history_uniform_naming_incomplete.txt", as.is=T)
aliquot_sheet = master_sheet %>% select(aliquot_uuid = uuid, sample_id = Barcode) %>%
  mutate(aliquot_id = sprintf("%s-%s", sample_id, aliquot_uuid),
         analyte_type = "DNA",
         analysis_type = "WGS",
         portion = 1) %>%
  filter(grepl("HK", sample_id))

### aliquots
aliquots = aliquot_sheet %>% select(sample_id, aliquot_uuid, aliquot_id, portion, analyte_type, analysis_type) %>% distinct()
print(sprintf("Exporting manifest as json files"))
write(jsonlite::toJSON(aliquots, pretty = T), file = sprintf("%s.%s", aliquots_file, json_ext))
write.table(aliquots, file = sprintf("%s.%s", aliquots_file, text_ext), sep="\t", row.names = F, col.names = T, quote = F)

# Generate *file* tsv containing: aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format.
# Load in samples from the Hong Kong data set.
hk_df = read.table(HK_file_path, header=T, stringsAsFactors = F)

# Retrieve the sample name to map to study center provided covariate sheet.
hk_df$sample_id = sub(".*/WG_ *(.*?) *_USPD.*", "\\1", hk_df$filenames)

# Create a new identifier on which to group mate pairs onto the same line (i.e., R1 and R2).
hk_df$read_group = paste(hk_df$Library_ID, hk_df$FlowCell_ID, hk_df$Lane_ID, sep = '-')

# Retrieve the file_name from the file_path. 
hk_df$file_name_single = sapply(strsplit(hk_df$filenames, "/"), "[[", 6)

# comma separated file_paths and file_names.
merged_hk_files = hk_df %>% 
  group_by(read_group) %>% 
  mutate(file_path = paste(filenames, collapse=","),
         file_name = paste(file_name_single, collapse=",")) %>% 
  select(-filenames, -file_name_single) %>% 
  distinct()

# Combine with study center provided covariate sheet.
master_sheet$Original_ID = gsub("-","_", master_sheet$Original_ID)
hk_map_df = merged_hk_files %>% 
  inner_join(master_sheet, by = c("sample_id" = "Original_ID")) %>%  
  ungroup()

# Generate uuids for each of the hong kong files.
# 24c6f54a-e7a2-4148-8335-045e3c74096e
set.seed(1)
hk_map_df$file_uuid = paste(stri_rand_strings(dim(hk_map_df)[1], 8, "[a-z0-9]"),
  stri_rand_strings(dim(hk_map_df)[1], 4, "[a-z0-9]"),
  stri_rand_strings(dim(hk_map_df)[1], 4, "[a-z0-9]"),
  stri_rand_strings(dim(hk_map_df)[1], 4, "[a-z0-9]"),
  stri_rand_strings(dim(hk_map_df)[1], 12, "[a-z0-9]"),
  sep = "-")

# Sanity check: make sure each is unique.
n_distinct(hk_map_df$file_uuid)

# Need to record the file_size and file_md5sum for these samples.
hk_map_df = hk_map_df %>% mutate(file_size = "NA",
                     file_md5sum = "NA",
                     aliquot_id = sprintf("%s-%s", Barcode, uuid))
# aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format.

## files
files = hk_map_df %>% select(aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format = File_Type) %>% distinct()
write(jsonlite::toJSON(files, pretty = T), file = sprintf("%s.%s", files_file, json_ext))
write.table(files, file = sprintf("%s.%s", files_file, text_ext), sep="\t", row.names = F, col.names = T, quote = F)


## Samples
samples =  %>% select(case_id, sample_id, legacy_sample_id, sample_type = sample_type_code) %>% distinct()

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
