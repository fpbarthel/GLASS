#######################################################
# Create manifest for GLASS whole exome sequencing samples.
# Date: 2018.10.19
# Author: Kevin J.
######    Creation of this manifest reflects the data as of: **10.19.2018**   ###########
#######################################################
# Local directory for github repo.
mybasedir = "/Users/johnsk/Documents/Life-History/GLASS-WG/"
setwd(mybasedir)

# Files with information about fastq information and barcodes.
life_history_barcodes = "/Users/johnsk/Documents/Life-History/GLASS-WG/data/ref/glass_wg_aliquots_mapping_table.txt"
wxs_clinical = "data/clinical-data/GLASS/full_glass_clinical_AM_20181019.xlsx"
wxs_bam_map = "data/ref/glass_analysis-paired_bam_map_wes.variant_20181019.txt"
glass_wxs_barcodes = "data/ref/glass_wxs_barcodes.txt"
wxs_revised_readgroups = "data/sequencing-information/glass-wxs/glass_wxs_readgroups_revised_rgid.txt"

#######################################################

# Necessary packages:
library(tidyverse)
library(openxlsx)
library(EnvStats)
library(stringr)
library(stringi)
library(DBI)

#######################################################
# Establish connection with Floris' database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

# We need to generate the following fields required by the SNV snakemake pipeline:
### Aliquots ####
glass_wxs_ids = read.delim(glass_wxs_barcodes, as.is=T)

# Perform a test run where we are only interested in the Columbia wxs samples.
aliquots_master = glass_wxs_ids %>% mutate(sample_barcode = substr(aliquot_barcode, 1, 15),
                                              sample_type = substr(aliquot_barcode, 14, 15),
                                              case_barcode = substr(aliquot_barcode, 1, 12),
                                              aliquot_uuid_short = substr(aliquot_barcode, 25, 30),
                                              aliquot_analyte_type = substr(aliquot_barcode, 19, 19),
                                              aliquot_portion = substr(aliquot_barcode, 17, 18),
                                              aliquot_analysis_type = substr(aliquot_barcode, 21, 23),
                                              aliquot_id_legacy = aliquot_name) %>%
  select(-aliquot_name) 

### aliquots
aliquots = aliquots_master %>% 
  select(aliquot_barcode, aliquot_id_legacy, sample_barcode, aliquot_uuid_short,
         aliquot_analyte_type, aliquot_analysis_type, aliquot_portion) %>% 
  distinct()


### Files ####
# Combine data to get filename to merge with readgroup information:
glass_wxs_bam_map = read.delim(wxs_bam_map, as.is=T)

# Trim off any whitespace that may have been introduced by manipulating in excel. Paranoia strikes deep.
glass_wxs_bam_map = data.frame(lapply(glass_wxs_bam_map, trimws), stringsAsFactors = F)

# Currently, I separate the readgroups so that it's easier to merge with Hoon's "tm_samplename" and "nm_samplename" strategy.
glass_readgroups <- read.delim(wxs_revised_readgroups, as.is=T)
glass_readgroups = data.frame(lapply(glass_readgroups, trimws), stringsAsFactors = F)

# Provide a map to connect the aliquot_id with the legacy_id.
glass_linker = aliquots_master %>% 
  select(aliquot_id_legacy, aliquot_barcode)

# Again, we had issues so making sure whitespace is not an issue...
glass_wxs_bam_map$tm_samplename <- trimws(glass_wxs_bam_map$tm_samplename, "r")
glass_wxs_bam_map$nm_samplename <- trimws(glass_wxs_bam_map$nm_samplename, "r") 
glass_linker$aliquot_id_legacy <- trimws(glass_linker$aliquot_id_legacy, "r")

# Separately merge by tumor and normal because of they are different columns in the bam map.
glass_bam_tm <- glass_wxs_bam_map %>% 
  inner_join(glass_linker, by=c("tm_samplename"="aliquot_id_legacy")) %>% 
  select(starts_with("tm"), aliquot_barcode) %>% 
  distinct()

# These numbers will changed depending on what's included in the analysis (i.e., whether "R1" only samples will be included or tossed).
glass_bam_nm <- glass_wxs_bam_map %>% 
  inner_join(glass_linker, by=c("nm_samplename"="aliquot_id_legacy")) %>% 
  select(starts_with("nm"), aliquot_barcode) %>% 
  distinct()

# Combine the readgroup information.
wxs_tm_df = glass_bam_tm %>% 
  inner_join(glass_readgroups, by=c("tm_bam_rgsm"="bam_rgsm"))
wxs_nm_df = glass_bam_nm %>% 
  inner_join(glass_readgroups, by=c("nm_bam_rgsm"="bam_rgsm"))

# Create function to grab the last element of each vector (varying length vectors) from file path.
last <- function(x) { return( x[length(x)] ) }
wxs_tm_df$file_name <- sapply(strsplit(wxs_tm_df$tm_bam, "/"), last)
wxs_tm_map_df <- wxs_tm_df %>% 
  mutate(file_format = "uBAM",
         file_path = tm_bam)
colnames(wxs_tm_map_df) <- gsub("tm_", "", colnames(wxs_tm_map_df))
wxs_nm_df$file_name <- sapply(strsplit(wxs_nm_df$nm_bam, "/"), last)
wxs_nm_map_df <- wxs_nm_df %>% 
  mutate(file_format = "uBAM",
         file_path =  wxs_nm_df$nm_bam)
colnames(wxs_nm_map_df) <- gsub("nm_", "", colnames(wxs_nm_map_df))

# Create final data.frame for input. This number will be slightly different from RGs because two samples/RGs were removed ("R1" only sample + NB).
wxs_map_df <- bind_rows(wxs_tm_map_df, wxs_nm_map_df)

# Order needs to be: aliquot_barcode, file_path, file_name, file_size, file_md5sum, file_format.
files = wxs_map_df %>% 
  select(aliquot_barcode, file_path, file_name, file_size = bam_size, file_md5sum = bam_md5, file_format) %>%
  distinct()

# Since one of the file sizes was an even number, one of the programs converted to scientific notation.
files$file_size[files$file_size==1.35e+10] = "13500000000"
  
### Cases ####
# Clinical data to extract subject sex and age at first diagnosis.
glass_wxs_clinical = readWorkbook(wxs_clinical, sheet = 1, startRow = 1, colNames = TRUE)
glass_wxs_clinical = data.frame(lapply(glass_wxs_clinical, trimws), stringsAsFactors = F)

# Sanity check because I was having issues with the data manipulated through excel. 
glass_wxs_clinical$patientid = trimws(glass_wxs_clinical$patientid, "r")
glass_wxs_bam_map$patientid  = trimws(glass_wxs_bam_map$patientid, "r")

# Combine data for which we already have bam files generated.
wxs_clinical_data = glass_wxs_clinical %>% 
  select(patientid, age.at.diagnosis, gender) %>% 
  inner_join(glass_wxs_bam_map, by="patientid") %>% 
  inner_join(glass_linker, by=c("tm_samplename"="aliquot_id_legacy"))

# Prepare data to be merged.
wxs_clinical_df = wxs_clinical_data %>% 
  mutate(age = `age.at.diagnosis`,
         sex = recode(gender, "M"="male","F"="female")) %>% 
  select(aliquot_barcode, age, sex)

# Select only those relevant fields. Cases will vary depending on which samples are included.
cases = wxs_clinical_df %>% 
  mutate(case_barcode = substring(aliquot_barcode, 1, 12), 
         case_project = substring(aliquot_barcode, 1, 4),
         case_source = substr(aliquot_barcode, 6, 7)) %>% 
  select(case_project, case_barcode, case_source, case_age_diagnosis_years=age, case_sex=sex) %>% 
  distinct()

### Samples ####
# Grab last two characrters of barcode. These samples have multi-sector sampling. So the 
# number of sample_types will be smaller than all of the tumors.
samples = wxs_map_df %>% 
  mutate(sample_type = substring(wxs_map_df$aliquot_barcode, 14, 15),
         case_barcode =  substring(wxs_map_df$aliquot_barcode, 1, 12),
         sample_barcode = substring(wxs_map_df$aliquot_barcode, 1, 15)) %>% 
  select(case_barcode, sample_barcode, sample_type) %>% 
  distinct()

### Readgroups ####
# Revised: Can't rename RGID, but pipeline is looking for old RGID.
readgroup_df = wxs_map_df %>% 
  mutate(readgroup_platform = bam_rgpl,
         legacy_readgroup_id = bam_rgid,
         readgroup_platform_unit = bam_rgpu,
         readgroup_date = bam_rgdt,
         readgroup_library = bam_rglb,
         readgroup_center = bam_rgcn,
         readgroup_sample_id = aliquot_barcode,
         readgroup_id = glass_rg_id)

# Finalize readgroup information in predefined order.
readgroups = readgroup_df %>% 
  mutate(readgroup_platform = "ILLUMINA") %>% 
  select(aliquot_barcode, readgroup_idtag=readgroup_id, readgroup_platform, readgroup_platform_unit, readgroup_library,
         readgroup_center, readgroup_idtag_legacy = legacy_readgroup_id) %>% 
  mutate(readgroup_sample_id = aliquot_barcode) %>% distinct()

## Filemap
tmp1 = files %>% select(file_name, aliquot_barcode)
tmp2 = readgroups %>% select(readgroup_idtag, readgroup_sample_id, aliquot_barcode)
files_readgroups = full_join(tmp1, tmp2) %>% select(file_name, readgroup_idtag, readgroup_sample_id)

### OUTPUT ####
# Only upload NEW cases. You can query the database and perform an antijoin.
extant_cases <- dbReadTable(con,  Id(schema="clinical",table="cases"))
new_cases = cases %>% 
  anti_join(extant_cases, by="case_barcode")
# New samples.
extant_samples <- dbReadTable(con,  Id(schema="biospecimen",table="samples"))
new_samples = samples %>% 
  anti_join(extant_samples, by="sample_barcode")
# New samples.
extant_aliquots <- dbReadTable(con,  Id(schema="biospecimen",table="aliquots"))
new_aliquots = aliquots %>% 
  anti_join(extant_aliquots, by="aliquot_barcode")

## Write to database.
dbWriteTable(con, Id(schema="clinical",table="cases"), new_cases, append=T)
dbWriteTable(con, Id(schema="biospecimen",table="samples"), new_samples, append=T)
dbWriteTable(con, Id(schema="biospecimen",table="aliquots"), new_aliquots, append=T)
dbWriteTable(con, Id(schema="biospecimen",table="readgroups"), readgroups, append=T)
dbWriteTable(con, Id(schema="analysis",table="files"), files, append=T)
dbWriteTable(con, Id(schema="analysis",table="files_readgroups"), files_readgroups, append=T)

# This is how to remove/delete rows from the database.
# rs <- dbSendStatement(con, "DELETE FROM clinical.cases WHERE case_barcode = ?")
# dbSendStatement(con, "DELETE FROM clinical.cases WHERE case_barcode = ? AND ")
# sapply(test_cases$case_barcode, function(case_barcode) dbBind(rs, list(case_barcode)))
# dbHasCompleted(rs)
# dbGetRowsAffected(rs)
# dbClearResult(rs)
