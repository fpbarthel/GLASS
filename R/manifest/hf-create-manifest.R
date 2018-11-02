#######################################################
# Create manifest for low pass Henry Ford sample (hGBM).
# Date: 2018.09.27
# Author: Kevin J., Floris B.
#######################################################
# Local directory for github repo.
mybasedir = "/Users/johnsk/Documents/Life-History/GLASS-WG"
setwd(mybasedir)

# Files with information about fastq information and barcodes.
hf_file_path = "data/sequencing-information/HenryFord/henry_ford_readgroups.txt"
life_history_filesize = "data/ref/henryford_filesize.txt"
life_history_md5 = "data/ref/henryford_md5.txt"
life_history_barcodes = "data/ref/glass_wg_aliquots_mapping_table.txt"
hf_clinical_info = "data/clinical-data/HF/HGBM-clinical-HF3016.xlsx"

#######################################################

library(tidyverse)
library(openxlsx)
library(rjson)
library(jsonlite)
library(listviewer)
library(stringi)
library(stringr)
library(DBI)

#######################################################
# Establish connection.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

# We need to generate the following fields required by the SNV snakemake pipeline:
# aliquots, files, cases, samples, pairs, and readgroups.

### Aliquots ####
life_history_ids = read.delim(life_history_barcodes, as.is=T)
aliquots_master = life_history_ids %>% mutate(aliquot_barcode = aliquot_id,
                                              sample_barcode = substr(aliquot_barcode, 1, 15),
                                              sample_type = substr(aliquot_barcode, 14, 15),
                                              case_barcode = substr(aliquot_barcode, 1, 12),
                                              aliquot_uuid_short = substr(aliquot_barcode, 25, 30),
                                              aliquot_analyte_type = substr(aliquot_barcode, 19, 19),
                                              aliquot_portion = substr(aliquot_barcode, 17, 18),
                                              aliquot_analysis_type = substr(aliquot_barcode, 21, 23),
                                              aliquot_id_legacy = legacy_aliquot_id) %>%
  select(-legacy_aliquot_id) %>% 
  filter(grepl("GLSS-HF", sample_barcode))

### aliquots
aliquots = aliquots_master %>% 
  select(aliquot_barcode, aliquot_id_legacy, sample_barcode, aliquot_uuid_short,
         aliquot_analyte_type, aliquot_analysis_type, aliquot_portion) %>% 
  distinct()

### Files ####
# Generate *file* tsv containing: aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format.
hf_file_name = read.table(hf_file_path, stringsAsFactors = F)
colnames(hf_file_name) = c("file_path", "@RG" , "RGID", "RGPL", "RGLB", "RGSM", "RGCN")
hf_file_name$file_name <- sapply(strsplit(hf_file_name$file_path, "/"), "[[", 6)
hf_filesize = read_tsv(life_history_filesize, col_names=c("file_size", "file_name"))
hf_md5 = read_tsv(life_history_md5, col_names=c("file_md5sum", "file_name"))
hf_df = hf_file_name %>% 
  inner_join(hf_md5, by="file_name") %>% 
  inner_join(hf_filesize, by="file_name") 

# Remove RG identifier before colon.
drop_prefix = function(x) {gsub(".*:","", x) }
hf_df[3:7] <- lapply(hf_df[3:7], drop_prefix)

# Retrieve the file_name from the file_path. 
hf_map_df = hf_df %>%  
  mutate(legacy_sample_id = sub(".*FORD- *(.*?) *_1.*", "\\1", file_name),
         file_format = "uBAM") %>% 
  inner_join(aliquots_master, by = c("legacy_sample_id" = "aliquot_id_legacy")) 

# Order needs to be: # aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format.
files = hf_map_df %>% 
  select(aliquot_barcode, file_name, file_size, file_md5sum, file_format) %>%
  distinct()


### Cases ####
# Clinical data to extract subject sex.
hf_clinical_data <- readWorkbook(hf_clinical_info, sheet = 1, startRow = 1, colNames = TRUE)

# Prepare data to be merged.
hf_clinical = hf_clinical_data %>% 
  mutate(age = `Age.@.dx`,
         sex = recode(Gender, "M"="male"),
         GLSS_case_id = "GLSS-HF-3016") %>% 
  distinct(GLSS_case_id, .keep_all = TRUE) %>% 
  select(GLSS_case_id, age, sex)

# Merge to lift over `sex` variable.
hf_map_df = hf_map_df %>% 
  mutate(case_id = substring(aliquot_id, 1, 12), 
         case_project = "GLSS") %>% 
  inner_join(hf_clinical, by=c("case_id" = "GLSS_case_id"))

# Select only those relevant fields.
cases = hf_map_df %>% 
  mutate(case_source = substr(aliquot_id,6,7)) %>% 
  select(case_project, case_barcode, case_source, case_age_diagnosis_years=age, case_sex=sex) %>% 
  distinct()


### Samples ####
# Grab last two characrters of barcode.
hf_map_df$sample_type = substring(hf_map_df$aliquot_id, 14, 15)

# Recode variables to match Floris' fields.
samples = hf_map_df %>% select(case_barcode, sample_barcode, sample_type) %>% distinct()


### Readgroups ####
# Necessary information: file_uuid, aliquot_id, RGID, RGPL, RGPU, RGLB, RGPI, RGDT, RGSM, RGCN.
# *** Note: had to generate a new time stamp because it was not found in the RG header for the bams. ****
# Revised: Can't rename RGID, but pipeline is looking for old RGID.
readgroup_df = hf_map_df %>% 
  mutate(readgroup_platform = "ILLUMINA",
         legacy_readgroup_id = RGID,
         readgroup_platform_unit = paste(sub(".*[0-9]{4}_ *(.*?) *_s.*", "\\1", hf_map_df$file_name),
                                         sub(".*s_ *(.*?) *_[A-Za-z]{2}.*", "\\1", hf_map_df$file_name), sep="."),
         readgroup_date = strftime(as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%dT%H:%M:%S%z"),
         readgroup_library = RGLB,
         readgroup_center = RGCN,
         readgroup_sample_id = aliquot_id, 
         readgroup_id = paste0(substring(readgroup_platform_unit, 1, 5), substring(readgroup_platform_unit, nchar(readgroup_platform_unit)-1, nchar(readgroup_platform_unit)), ""))
         
# Finalize readgroup information in predefined order.
readgroups = readgroup_df %>% 
  select(aliquot_barcode, readgroup_idtag=readgroup_id, readgroup_platform, readgroup_platform_unit, readgroup_library,
         readgroup_center) %>% 
  mutate(readgroup_sample_id = aliquot_barcode) %>% distinct()

## Filemap
tmp1 = files %>% select(file_name, aliquot_barcode)
tmp2 = readgroups %>% select(readgroup_idtag, readgroup_sample_id, aliquot_barcode)
files_readgroups = full_join(tmp1, tmp2) %>% select(file_name, readgroup_idtag, readgroup_sample_id)


## Write to database.
dbWriteTable(con, Id(schema="clinical",table="cases"), cases, append=T)
dbWriteTable(con, Id(schema="biospecimen",table="samples"), samples, append=T)
dbWriteTable(con, Id(schema="biospecimen",table="aliquots"), aliquots, append=T)
dbWriteTable(con, Id(schema="biospecimen",table="readgroups"), readgroups, append=T)
dbWriteTable(con, Id(schema="analysis",table="files"), files, append=T)
dbWriteTable(con, Id(schema="analysis",table="files_readgroups"), files_readgroups, append=T)
