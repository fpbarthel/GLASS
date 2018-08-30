#######################################################
# To create aliquot names for glioma life history.
# Date: 2018.08.30
# Author: Kevin J
#######################################################
# Local directory for github repo.
mybasedir = here::here()
setwd(mybasedir)

# Hong Kong samples sequenced by Novogene.
hong_kong_path = "data/sequencing-information/HK/hong-kong-sample-maps.txt"
# TCGA LGG and GBM datasets.
tcga_path = "data/clinical-data/TCGA/LGG-GBM-samples.tsv"
# Low-pass whole genome sequencing for MD Anderson samples.
rv_jdg_path = "data/clinical-data/Roel-JDG/Roel_JDG-original_bam_map.txt"
# Glioblastoma case from the Henry Ford data set.
hgbm_path = "data/clinical-data/HF/HGBM-original_bam_map.txt"
# The Australian subject (Northern Sydney Cancer Centre) that had multile recurrences and a bone metastasis.
ns_path ="data/sequencing-information/NS/Sample_Identifiers.txt"
# Different tissue source sites from TCGA.
tss_path = "data/sequencing-information/tissueSourceSite.txt"

# MDA: Files with information about fastq information and barcodes, but also aliquot names.
mda_batch1_file_path = "data/sequencing-information/MDACC/file_list_C202SC18030593_batch1_n14_20180405.tsv"
mda_batch2_file_path = "data/sequencing-information/MDACC/file_list_C202SC18030593_batch2_n94_20180603.tsv"
mda_batch3_file_path = "data/sequencing-information/MDACC/file_list_C202SC18030593_batch2_n13_20180716.tsv"

# Clinical dataset priovided by Kristin Alfaro-Munoz at MD Anderson.
# Modified to generate a sample_type (TP, R1, R2, R3) using the clinical information. See github issue #16.
mda_master_path = "data/clinical-data/MDACC/MDA-Clinical-Dataset/Master Log for WGS_sampletype.20180630.xlsx"

# 2018.07.06 Katie Shao (Novogene) provided the following sheet linking libraryID with submitted sample names in the email:
# "Re:RE: Novogene Project Report - Confirm -C202SC18030593-Davis-MDACC-156-libseq-hWGS-WOBI-NVUS2018022505[J1-C202SC18030593-1000]"
novogene_sample_path = "data/sequencing-information/MDACC/Novogene_SIF_14.xlsx"

#######################################################
# Necessary packages.
library(tidyverse)
library(openxlsx)
library(stringr)
library(stringi)

#######################################################
# Inspect how each cohort currently barcodes each file. 
hk_df = read.table(hong_kong_path, header=T)
tcga_df = read_tsv(tcga_path)
rv_jdg_df = read_tsv(rv_jdg_path)
hgbm_df = read_tsv(hgbm_path)
ns_df = read_tsv(ns_path)

# The sample-specific barcode should include the following:
# project_id, tissue source site, subject_id, sample_type, portion, analysis_type, and legacy_aliquot_id
hk_sample_map = hk_df %>% 
  mutate(project_id = "GLSS",
         tss = "HK",
         subject_id = recode_factor(Patient_ID, "CHUK5" = "0005", "CHUK4" = "0004", "CHUK3" = "0003", "CHUK2" = "0002", "CHUK1" = "0001"),
         sample_type = recode_factor(Sample_Type, "Blood" = "NB", "Primary" = "TP", "Recurrence" = "R1"),
         portion = "01D",
         analysis_type = "WGS",
         legacy_aliquot_id = Verhaak_Sample_ID) %>%
  unite(glass_aliquot_id, c(project_id, tss, subject_id, sample_type, portion, analysis_type), sep = "-", remove = FALSE) %>% 
  select(glass_aliquot_id, legacy_aliquot_id)


### MDA (Novogene cohort) barcode generation ####
# Novogene has processed 120 samples. In our master sheet normal blood samples were not provided. We need them for barcodes.
# Inspect both batch1, batch2, and batch 3 sequencing data for tumor/normal samples.
mda_batch1_df <- read.table(mda_batch1_file_path, col.names="relative_file_path", stringsAsFactors = F)
mda_batch1_df$filename = sapply(strsplit(mda_batch1_df$relative_file_path, "/"), "[[", 3)
mda_batch2_df <- read.table(mda_batch2_file_path, col.names="relative_file_path", stringsAsFactors = F)
mda_batch2_df$filename = sapply(strsplit(mda_batch2_df$relative_file_path, "/"), "[[", 3)
mda_batch3_df <- read.table(mda_batch3_file_path, col.names="relative_file_path", stringsAsFactors = F)
mda_batch3_df$filename = sapply(strsplit(mda_batch3_df$relative_file_path, "/"), "[[", 3)

# Replace the placeholder for pwd "./" from unix cmd: "find -name ".fq.gz" in parent directory of fastqs. 
mda_batch1_df$working_file_path <- gsub("^", "/fastscratch/johnsk/GLASS-WG/mdacc/C202SC18030593_batch1_n14_20180405", mda_batch1_df$relative_file_path)
mda_batch2_df$working_file_path <- gsub("^", "/fastscratch/johnsk/GLASS-WG/mdacc/C202SC18030593_batch1_n94_20180603", mda_batch2_df$relative_file_path)
mda_batch3_df$working_file_path <- gsub("^", "/fastscratch/johnsk/GLASS-WG/mdacc/C202SC18030593_batch2_n13_20180716", mda_batch3_df$relative_file_path)

# Combine these three sequencing data sets. There should be 552 files.
mda_df <- bind_rows(mda_batch1_df, mda_batch2_df, mda_batch3_df)

# Katie Shao (Novogene) provided the filename and sample linker file.
novogene_linker = readWorkbook(novogene_sample_path, sheet = 1, startRow = 18, colNames = TRUE)
novogene_linker_unique <- novogene_linker[1:121, ]

# Retrieve only the blood samples ("10D"). Blood samples not contained in master sheet.
normal_blood_samples = novogene_linker_unique[grep("[b|B]lood", novogene_linker_unique$`*SampleName`), ]

# Retrieve the subject ID in same format as tumor samples.
normal_blood_samples$SubjectID = sapply(strsplit(normal_blood_samples$`*SampleName`, "-"), "[", 3)

# Create a barcode for available blood samples (n=39).
mda_normal_blood_map <- normal_blood_samples %>% 
  mutate(project_id = "GLSS",
         tss = "MD",
         subject_id = SubjectID,
         portion = "01D",
         analysis_type = "WGS",
         sample_type = "NB",
         legacy_aliquot_id = `*SampleName`) %>%   
  unite(glass_aliquot_id, c(project_id, tss, subject_id, sample_type, portion, analysis_type), sep = "-", remove = FALSE) %>% 
  select(glass_aliquot_id, legacy_aliquot_id)

### Use master sheet containing tumor sample information, necessary to build barcodes.####
# Kristin Alfaro-Munoz kindly pointed me to the sequencing identifier link to these samples. 
mda_master = readWorkbook(mda_master_path, sheet = 2, startRow = 1, colNames = TRUE)

# The master sheet from MDA only contains normal samples. Sum should equal 81 (tumor samples).
sum(novogene_linker_unique$'*SampleName'%in%mda_master$Jax.Lib.Prep.Customer.Sample.Name)

# Extract the 4-digit SUBJECT identifier.
mda_master$SubjectID = sapply(strsplit(mda_master$Jax.Lib.Prep.Customer.Sample.Name, "-"), "[", 3)

# Create a new MD Anderson to map the old to new sample names.
mda_tumor_sample_map = mda_master[1:81,] %>% 
  mutate(project_id = "GLSS",
         tss = "MD",
         legacy_aliquot_id = Jax.Lib.Prep.Customer.Sample.Name,
         subject_id = SubjectID,
         portion = "01D",
         analysis_type = "WGS") %>% 
  unite(glass_aliquot_id, c(project_id, tss, subject_id, sample_type, portion, analysis_type), sep = "-", remove = FALSE) %>% 
  select(glass_aliquot_id, legacy_aliquot_id)

# Combine all tumor and normal samples together.
mda_all_samples_map <- bind_rows(mda_normal_blood_map, mda_tumor_sample_map)



##### TCGA #####
# Create a new TCGA to map the old to new sample names.
tcga_sample_map = tcga_df %>% 
  rename(legacy_aliquot_id = aliquot_id) %>% 
  separate(legacy_aliquot_id,  c("project_id", "tss", "subject_id", "sample_type", "portion", "plate_id", "center"), remove = FALSE) %>% 
  mutate(glass_sample_type = recode_factor(sample_type, "01A" = "TP", "01B" = "TP", "02A" = "R1", "02B" = "R2", "10A" = "NB", "10B" = "NB", "10D" = "NB"),
         analysis_type = "WGS") %>% 
unite(glass_aliquot_id, c(project_id, tss, subject_id, glass_sample_type, portion, analysis_type), sep = "-", remove = FALSE) %>% 
  select(glass_aliquot_id, legacy_aliquot_id)
  

##### RV-JDG #####
# Modify barcode to make it the same as the other studies with only 4 characters.
# Replacing Neuro-oncologist with "LP" identifier to indicate "Low Pass" and disambiguate from other MD set.
rv_jdg_df$AnalysisID = gsub("JDG", "LP", rv_jdg_df$AnalysisID)

# Create a new RV_JDG sample map to connect with the old sample names.
rvjdg_sample_map = rv_jdg_df %>%   
  separate(col=AnalysisID, into = c("Study", "PtIdentifier", "WrongTissueID", "Analyte"), sep="-") %>% 
  separate(PtIdentifier, into = c("Text", "Num"), sep = "(?<=[A-Za-z])(?=[0-9])") %>% 
  mutate_at("Num", str_pad, width = 2, side='left', pad = 0) %>% 
  unite(PtIdentifier, c(Text,Num), sep = "") %>% 
  mutate(sample_type = recode_factor(Tissue.Type.Revised, "first" = "TP", "second" = "R1", "normal" = "NB"),
         project_id = "GLSS",
         tss = "MD",
         subject_id = PtIdentifier,
         portion = Analyte,
         analysis_type = "WGS",
         legacy_aliquot_id = samplename) %>% 
  group_by(PatientID) %>% 
  filter(n() >= 3)  %>%  # Extract only trios.
  unite(glass_aliquot_id, c(project_id, tss, subject_id, sample_type, portion, analysis_type), sep = "-", remove = FALSE) %>% 
  ungroup() %>% 
  select(glass_aliquot_id, legacy_aliquot_id)


##### HGBM ####
# Create a new hGBM sample map to connect with the old sample names.
hgbm_sample_map = hgbm_df %>% 
  filter(samplename!="HF-3177-10-01D") %>% # Removing the second normal blood sample from this same patient.
  separate(samplename,  c("tss", "subject_id", "sample_type", "portion"), remove = FALSE) %>% 
  mutate(project_id = "GLSS", 
         tss = "HF",
         analysis_type = "WGS",
         subject_id = substring(PatientID, 4, 7),
         sample_type = recode_factor(SampleType2, "P" = "TP", "R" = "R1", "N" = "NB"),
         legacy_aliquot_id = samplename) %>% 
unite(glass_aliquot_id, c(project_id, tss, subject_id, sample_type, portion, analysis_type), sep = "-", remove = FALSE) %>% 
  select(glass_aliquot_id, legacy_aliquot_id)


###### NS subject ##### 
ns_sample_map = ns_df %>% 
  separate(Sample, c("original_id", "lane_id"), sep="_L[0-9]{1}_", extra = "merge") %>% 
  distinct(original_id) %>% 
  mutate(project_id = "GLSS",
         subject_id = "0001",
         tss ="NS",
         portion = "01D",
         sample_type = c("TP", "R1", "R2", "M1", "NB"),
         analysis_type = "WGS",
         legacy_aliquot_id = original_id) %>% 
  unite(glass_aliquot_id, c(project_id, tss, subject_id, sample_type, portion, analysis_type), sep = "-", remove = FALSE) %>% 
  select(glass_aliquot_id, legacy_aliquot_id)


#############################################
# Combine 6 datasets into master table with each linker
#############################################
life_history_barcodes = bind_rows(hk_sample_map, mda_tumor_sample_map, tcga_sample_map, hgbm_sample_map, rvjdg_sample_map, mda_normal_blood_map, ns_sample_map)

## Addition of random string to use as a shortened identifier.
# Create a random string identifier for each SAMPLE in the Life History study.
set.seed(1)
uuid_life_history = as.data.frame(stri_rand_strings(dim(life_history_barcodes)[1], 6, "[A-Z0-9]"))
colnames(uuid_life_history) = "uuid"

# Sanity check: make sure each is unique.
ifelse(n_distinct(uuid_life_history)==dim(uuid_life_history)[1], 
       message("UUIDs are unique."), 
       message("WARNING! Not unique"))

# Write final combined dataset set for a data freeze.
life_history_all_samples = bind_cols(life_history_barcodes, uuid_life_history) %>% 
  unite(aliquot_id, c(glass_aliquot_id, uuid), sep = "-", remove = TRUE)

# Write file to be uploaded to GLASS-WG github page.
write.table(life_history_all_samples, file='data/ref/glass_wg_sample_mapping_table.txt', quote=FALSE, sep='\t', row.names = F)

### Create dictionaries for different tissue centers. ####
# Subset list of TissueSourceSites from TCGA plus add new rows for other datasets.
tss_info = read_tsv(tss_path)

# Retrieve list of TissueSourceSites for future use.
tcga_map = tcga_df %>% 
  rename(legacy_aliquot_id = aliquot_id) %>% 
  separate(legacy_aliquot_id,  c("project_id", "tss", "subject_id", "sample_type", "portion", "plate_id", "center"), remove = FALSE) %>% 
  distinct(tss, as.character)
tcga_tss <- unique(tcga_map$tss)

# Collapse space into "." separated. 
colnames(tss_info) = gsub( " ", ".", colnames(tss_info))

# Create more formal definitions.
tss_glass = data.frame(tss_code   = c("HK", "MD", "HF", "NS"),
                source_site  = c("Chinese University of Hong Kong", "MD Anderson Cancer Center", "Henry Ford Hospital", "Northern Sydney Cancer Centre"),
                study_name = c("Lower Grade Glioma / Glioblastoma", "Lower Grade Glioma / Glioblastoma", "Glioblastoma", "Metastatic Gliosarcoma"),
                bcr      = rep("GLASS", 4))
# Only those TCGA centers that are present in the NB, TP, R1/R2 samples.
tss_tcga_filtered = tss_info %>% 
  filter(TSS.Code %in% tcga_tss)

# Create new tissue source site dictionary.
tissue_source_site_dict = bind_rows(tss_tcga_filtered, tss_glass)
write.table(tissue_source_site_dict, file='data/ref/life-history-tissue-source-sites.txt', quote=FALSE, sep='\t', row.names = F)

# Create sample type dictionary. 
sample_type_dict = data.frame(Code   = c("NB", "TP", "R1", "R2", "R3", "M1"),
                       Definition  = c("Blood Derived Normal", "Primary Tumor", "First Recurrence Tumor", "Second Recurrence Tumor", "Third Reccurence Tumor", "Bone Metastasis"))
write.table(sample_type_dict, file='data/ref/life-history-sample-type-dict.txt', quote=FALSE, sep='\t', row.names = F)
