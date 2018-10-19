#######################################################
# Create manifest for Hong Kong samples (GLASS).
# Date: 2018.09.27
# Author: Kevin J., Floris B.
#######################################################
# Local directory for github repo.
mybasedir = "/Users/johnsk/Documents/Life-History/GLASS-WG"
setwd(mybasedir)

# Files with information about fastq information and barcodes.
life_history_barcodes = "data/ref/glass_wg_aliquots_mapping_table.txt"
life_history_filesize = "data/ref/hongkong_filesize.txt"
life_history_md5 = "data/ref/hongkong_md5checksums.txt"
hk_file_path = "data/sequencing-information/HK/hong_kong_file_list_complete.20180813.tsv"
hk_clinical_info = "data/clinical-data/HK/recurrent.glioma.HongKong.20161221.xlsx"

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
  filter(grepl("GLSS-HK", sample_barcode))

### aliquots
aliquots = aliquots_master %>% 
  select(aliquot_barcode, aliquot_id_legacy, sample_barcode, aliquot_uuid_short,
         aliquot_analyte_type, aliquot_analysis_type, aliquot_portion) %>% 
  distinct()


### Files ####
# Generate *file* tsv containing: aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format.
hk_file_name = read.table(hk_file_path, col.names="file_name", stringsAsFactors = F)
hk_filesize = read_tsv(life_history_filesize, col_names=c("file_size", "file_name"))
hk_md5 = read_tsv(life_history_md5, col_names=c("file_md5sum", "file_name"))
hk_df = hk_file_name %>% 
  inner_join(hk_md5, by="file_name") %>% 
  inner_join(hk_filesize, by="file_name") 

# Retrieve the original sample name to map to study center provided covariate sheet.
hk_df_meta = hk_df %>% 
  mutate(verhaak_sample_id = sub(".*WG_ *(.*?) *_USPD.*", "\\1", file_name),
         library_id = sub(".*_ *(.*?) *_H.*", "\\1", file_name), 
         flowcell_id = substr(file_name, nchar(file_name)-19, nchar(file_name)-11),
         lane_id = substr(file_name, nchar(file_name)-8, nchar(file_name)-8))

# Reformat so that the aliquot_id_legacy matches between the sample_map and master aliquot file.
hk_df_meta$legacy_aliquot_id = gsub("_","-", hk_df_meta$verhaak_sample_id)
hk_map_df <- hk_df_meta %>% 
  full_join(aliquots_master, by=c("legacy_aliquot_id" = "aliquot_id_legacy")) %>% 
  select(-verhaak_sample_id) %>% 
  ungroup()

# Need to record the file_size and file_md5sum for these samples.
files = hk_map_df %>% 
  mutate(file_format = "FASTQ") %>% 
  select(aliquot_barcode, file_name, file_size, file_md5sum, file_format) %>%
  distinct()

### Cases ####
# Clinical data to extract subject sex.
hk_clinical_data <- readWorkbook(hk_clinical_info, sheet = 1, startRow = 1, colNames = TRUE)

# Prepare data to be merged.
hk_clinical = hk_clinical_data %>% 
  mutate(legacy_sample_id=Name,
         age = Age,
         sex = recode(Sex, "M"="male", "F"="female")) %>% 
  mutate(GLSS_case_id=recode(legacy_sample_id, "CHUK1"="GLSS-HK-0001", "CHUK2"="GLSS-HK-0002", "CHUK3"="GLSS-HK-0003", "CHUK4"="GLSS-HK-0004", "CHUK5"="GLSS-HK-0005")) %>% 
  select(GLSS_case_id, legacy_sample_id, age, sex)

# Merge to lift over `sex` variable.
hk_map_df = hk_map_df %>% 
  mutate(case_id = substring(aliquot_id, 1, 12), 
         case_project = "GLSS") %>% 
  inner_join(hk_clinical, by=c("case_id" = "GLSS_case_id"))

# Select only those relevant fields.
cases = hk_map_df %>% 
  mutate(case_source = substr(aliquot_id,6,7)) %>% 
  select(case_project, case_barcode, case_source, case_age_diagnosis_years=age, case_sex=sex) %>% 
  distinct()

### Samples ####
# Grab last two characrters of barcode.
hk_map_df$sample_type = substring(hk_map_df$aliquot_id, 14, 15)

# Recode variables to match Floris' fields.
samples = hk_map_df %>% 
  select(case_barcode, sample_barcode, sample_type) %>% 
  distinct()

### Readgroups ####
readgroup_df = hk_map_df %>% 
  mutate(readgroup_platform = "ILLUMINA",
         readgroup_platform_unit = paste(substr(file_name, nchar(file_name)-19, nchar(file_name)-11), 
                 substr(file_name, nchar(file_name)-8, nchar(file_name)-8), sep="."),
         readgroup_library = sub(".*_ *(.*?) *_H.*", "\\1", file_name),
         readgroup_date = strftime(as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%dT%H:%M:%S%z"), 
         readgroup_sample_id = aliquot_id,
         readgroup_center = "NVGN_HK",
         readgroup_id = paste0(substring(readgroup_platform_unit, 1, 5), substring(readgroup_platform_unit, nchar(readgroup_platform_unit)-1, nchar(readgroup_platform_unit)), "")) 

# Finalize readgroup information in predefined order.
readgroups = readgroup_df %>% 
  select(aliquot_barcode, readgroup_idtag=readgroup_id, readgroup_platform, readgroup_platform_unit, readgroup_library,
         readgroup_center) %>% 
  mutate(readgroup_sample_id = aliquot_barcode) %>% distinct()

## Filemap
files_readgroups = hk_map_df %>% 
  mutate(readgroup_idtag = paste(substring(hk_map_df$flowcell_id, 1, 5), hk_map_df$lane_id, sep = "."), 
  readgroup_sample_id = aliquot_barcode) %>% 
  select(file_name, readgroup_idtag, readgroup_sample_id)


## Write to database.
dbWriteTable(con, Id(schema="clinical",table="cases"), cases, append=T)
dbWriteTable(con, Id(schema="biospecimen",table="samples"), samples, append=T)
dbWriteTable(con, Id(schema="biospecimen",table="aliquots"), aliquots, append=T)
dbWriteTable(con, Id(schema="biospecimen",table="readgroups"), readgroups, append=T)
dbWriteTable(con, Id(schema="analysis",table="files"), files, append=T)
dbWriteTable(con, Id(schema="analysis",table="files_readgroups"), files_readgroups, append=T)

