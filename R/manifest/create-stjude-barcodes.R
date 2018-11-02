#######################################################
# Prepare barcodes for the St. Jude paired primary-recurrent bams.
# Date: 2018.10.31
# Author: Kevin J.
#######################################################

# St. Jude file map derived from file that was shared by Roel via email in Spring 2018.
st_jude_trios = "/Users/johnsk/Documents/Life-History/GLASS-WG/data/st-jude-data/st-jude_glioma_trio_data_wgs.xlsx"

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

# Read in existing aliquots to determine which uuids are already in use.
extant_aliquots <- dbReadTable(con,  Id(schema="biospecimen",table="aliquots"))
extant_aliquots$aliquot_uuid_short

# File map for the St. Jude samples
st_jude_trios = readWorkbook(st_jude_trios, sheet = 1, startRow = 1, colNames = TRUE)

# Filter out SJHGG003 since the two tumor samples came from separte brain regions and only retain bam file rows.
sj_pairs = st_jude_trios %>% 
  filter(!grepl(".bam.bai", file_path)) %>% 
  filter(!grepl("SJHGG003", subject_name)) %>% 
  filter(sj_datasets!="Clinical Pilot")
  
# Create a new barcode for the St. Jude samples.
set.seed(100)
stjude_wgs_barcodes = sj_pairs %>% 
  mutate(patient_number = substr(subject_name, 6, 8),
         portion_id = "01D", 
         glass_uuid = stri_rand_strings(dim(sj_pairs)[1], 6, "[A-Z0-9]"), 
         project = "PCGP",
         cohort_analysis.short = "SJ", 
         analysis_type = "WGS",
         sample_type = recode(sample_type, "Germline"="NB", "Diagnosis"="TP", "Relapse"="R1", "Autopsy"="R1")) %>% 
  mutate_at("patient_number", str_pad, width = 4, side='left', pad = 0) %>% 
  unite(aliquot_barcode, c(project, cohort_analysis.short, patient_number, sample_type, portion_id, analysis_type, glass_uuid), sep = "-", remove = FALSE) %>% 
  select(aliquot_barcode, project, cohort_analysis.short, patient_number, sample_type, portion_id, analysis_type, glass_uuid, legacy_aliquot_id = sample_name, file_type, file_path)

# Check to see whether any of the new uuids were in the GLSS database in case they are ever merged.
sum(extant_aliquots$aliquot_uuid_short%in%stjude_wgs_barcodes$glass_uuid)
n_distinct(stjude_wgs_barcodes$glass_uuid)

# Output table to be used in the generation of a WXS bam manifest.
write.table(stjude_wgs_barcodes, file = "/Users/johnsk/Documents/Life-History/GLASS-WG/data/ref/stjude_wgs_barcodes.txt", sep="\t", row.names = F, col.names = T, quote = F)
