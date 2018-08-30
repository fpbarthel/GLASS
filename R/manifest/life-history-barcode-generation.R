#######################################################
# To create sample barcodes for glioma life history.
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
# Glioblastoma from the Henry Ford data set.
hgbm_path = "data/clinical-data/HF/HGBM-original_bam_map.txt"
# The Australian subject (Northern Sydney Cancer Centre) that had multile recurrences and a bone metastasis.
ns_path ="data/sequencing-information/NS/Sample_Identifiers.txt"
# Different tissue source sites from TCGA.
tss_path = "data/sequencing-information/tissueSourceSite.txt"

# MDA: Files with information about fastq information and barcodes.
# Actual fastq files are stored: /fastscratch/johnsk/GLASS-WG/mdacc/
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
HK_df = read.table(hong_kong_path, header=T)
TCGA_df = read_tsv(tcga_path)
RV_JDG_df = read_tsv(rv_jdg_path)
hGBM_df = read_tsv(hgbm_path)
ns_df = read_tsv(ns_path)

# The sample-specific barcode should include the following:
# Project (TCGA/GLSS), Tissue Source Site (##), Subject Code (4), Sample Code (8), Tissue Code (01/02/10)

# Need to recode the Hong Kong sample names to a [:alnum:]{8}. 
HK_df$Verhaak_Sample_ID = as.character(HK_df$Verhaak_Sample_ID)
HK_df$New_Sample_ID = gsub("-", "", (gsub(".*S", "S", HK_df$Verhaak_Sample_ID)))
HK_df$New_Sample_ID = gsub("bld", "NB", HK_df$New_Sample_ID)
HK_df$New_Sample_ID = gsub("tis", "R1", HK_df$New_Sample_ID)
HK_df$New_Sample_ID = gsub("T3", "R1", HK_df$New_Sample_ID)
HK_df$New_Sample_ID = str_pad(HK_df$New_Sample_ID, width = 8, side = "right", pad = 0)

# Create a new data frame to be able to map old SampleID to new sample ID.
HongKong_sample_map = HK_df %>% 
  mutate(SeqID = "GLSS") %>% 
  mutate(TSS = "HK") %>% 
  mutate(SubjectCode = recode_factor(Patient_ID, "CHUK5" = "0005", "CHUK4" = "0004", "CHUK3" = "0003", "CHUK2" = "0002", "CHUK1" = "0001")) %>%
  rename(SampleCode = New_Sample_ID) %>%
  mutate(TissueCode = recode_factor(Sample_Type, "Blood" = "NB", "Primary" = "TP", "Recurrence" = "R1")) %>%
  unite(Barcode, SeqID:TissueCode, sep = "-", remove = FALSE)


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

# Create an old sample.id for these subjects to be linked. No longer necessary, but kept in case it's helpful.
mda_batch1_df$old_sample_id <- gsub( "_USPD.*$", "", mda_batch1_df$filename)
mda_batch2_df$old_sample_id <- gsub( "_USPD.*$", "", mda_batch2_df$filename)
mda_batch3_df$old_sample_id <- gsub( "_USPD.*$", "", mda_batch3_df$filename)

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
  mutate(SeqID = "GLSS") %>% 
  mutate(TSS = "MD") %>% 
  rename(SubjectCode = SubjectID) %>% 
  mutate(TissueCode = "NB") %>% 
  mutate(Original_ID = `*SampleName`) %>%   
  unite(Barcode, c(SeqID, TSS, SubjectCode, TissueCode), sep = "-", remove = FALSE) %>% 
  select(Original_ID, Barcode, SeqID, TSS, SubjectCode, TissueCode) %>% 
  distinct()

### Use master sheet containing tumor sample information, necessary to build barcodes.####
# Kristin Alfaro-Munoz kindly pointed me to the sequencing identifier link to these samples. 
mda_master = readWorkbook(mda_master_path, sheet = 2, startRow = 1, colNames = TRUE)

# The master sheet from MDA only contains normal samples. Sum should equal 81 (tumor samples).
sum(novogene_linker_unique$'*SampleName'%in%mda_master$Jax.Lib.Prep.Customer.Sample.Name)

# Extract the 4-digit SUBJECT identifier.
mda_master$SubjectID = sapply(strsplit(mda_master$Jax.Lib.Prep.Customer.Sample.Name, "-"), "[", 3)

# Create a new MD Anderson to map the old to new sample names.
mda_tumor_sample_map = mda_master[1:81,] %>% 
  mutate(SeqID = "GLSS") %>% 
  mutate(TSS = "MD") %>% 
  rename(SubjectCode = SubjectID) %>% 
  rename(TissueCode = sample_type) %>% 
  mutate(Original_ID = Jax.Lib.Prep.Customer.Sample.Name) %>% 
  unite(Barcode, c(SeqID, TSS, SubjectCode, TissueCode), sep = "-", remove = FALSE) %>% 
  select(Original_ID, Barcode, SeqID, TSS, SubjectCode, TissueCode)

# Combine all tumor and normal samples together.
MDA_all_samples_map <- bind_rows(mda_normal_blood_map, mda_tumor_sample_map)



##### TCGA #####
# TCGA tissue source sites.
unique(sapply(strsplit(TCGA_df$aliquot_id, "-"), "[", 2))

# Create a new TCGA to map the old to new sample names.
TCGA_sample_map = TCGA_df %>% 
  mutate(SeqID = "TCGA") %>% 
  mutate(TSS = sapply(strsplit(TCGA_df$aliquot_id, "-"), "[", 2)) %>% 
  mutate(SubjectCode = substring(aliquot_id, 9, 12)) %>% 
  mutate(TissueCode = recode_factor(substring(aliquot_id, 14, 16), "01A" = "TP", "01B" = "TP", "02A" = "R1", "02B" = "R2", "10A" = "NB", "10B" = "NB", "10D" = "NB")) %>% 
unite(Barcode, c(SeqID, TSS, SubjectCode, TissueCode), sep = "-", remove = FALSE)
  
##### RV-JDG #####
# This dataset contains more samples, including some unmatched recurrences. 
glimpse(RV_JDG_df)   

# Modify barcode to make it the same as the other studies with only 4 characters.
# Replacing Neuro-oncologist with "LP" identifier to indicate "Low Pass" and disambiguate from other MD set.
RV_JDG_df$AnalysisID = gsub("JDG", "LP", RV_JDG_df$AnalysisID)

# Create a new RV_JDG sample map to connect with the old sample names.
RVJDG_sample_map = RV_JDG_df %>%   
  mutate(SeqID = "GLSS") %>% 
  mutate(TSS = "MD") %>% 
  separate(col=AnalysisID, into = c("Study", "PtIdentifier", "WrongTissueID", "Analyte"), sep="-") %>% 
  separate(PtIdentifier, into = c("Text", "Num"), sep = "(?<=[A-Za-z])(?=[0-9])") %>% 
  mutate_at("Num", str_pad, width = 2, side='left', pad = 0) %>% 
  unite(PtIdentifier, c(Text,Num), sep = "") %>% 
  rename(SubjectCode = PtIdentifier) %>% 
  mutate(TissueCode = recode_factor(Tissue.Type.Revised, "first" = "TP", "second" = "R1", "normal" = "NB")) %>% 
  group_by(PatientID) %>% 
  filter(n() >= 3)  %>%  # Extract only trios.
  unite(Barcode, c(SeqID, TSS, SubjectCode, TissueCode), sep = "-", remove = TRUE) %>% 
  ungroup()

##### HGBM ####
# Create a new hGBM sample map to connect with the old sample names.
HGBM_sample_map = hGBM_df %>% 
  filter(samplename!="HF-3177-10-01D") %>% # Removing the second normal blood sample from this same patient.
  mutate(SeqID = "GLSS") %>% 
  mutate(TSS = "HF") %>% 
  mutate(SubjectCode = substring(PatientID, 4, 7)) %>% 
  mutate(TissueCode = recode_factor(SampleType2, "P" = "TP", "R" = "R1", "N" = "NB")) %>% 
  unite(Barcode, c(SeqID, TSS, SubjectCode, TissueCode), sep = "-", remove = TRUE)

###### AU subject ##### 
ns_sample_map = ns_df %>% 
  separate(Sample, c("Original_ID", "lane_id"), sep="_L[0-9]{1}_", extra = "merge") %>% 
  distinct(Original_ID) %>% 
  mutate(Barcode = sprintf("GLSS-NS-0001-%s", c("TP", "R1", "R2", "M1", "NB")),
         SubjectCode = "GLSS-NS-0001",
         TissueSourceSite ="NS",
         TissueCode = c("TP", "R1", "R2", "M1", "NB"))


#############################################
# Combine 6 datasets into master table with each linker
#############################################
HongKong_merge = HongKong_sample_map %>% select(Verhaak_Sample_ID, Barcode) %>% rename(Original_ID = Verhaak_Sample_ID)
MDA_tumor_merge = MDA_all_samples_map %>% filter(TissueCode!="NB") %>% select(Original_ID, Barcode) %>% filter()
TCGA_merge = TCGA_sample_map %>%  select(aliquot_id, Barcode)  %>% rename(Original_ID = aliquot_id)
HGBM_merge = HGBM_sample_map %>%  select(samplename, Barcode)  %>% rename(Original_ID = samplename)
RVJDG_merge = RVJDG_sample_map %>%  select(samplename, Barcode)  %>% rename(Original_ID = samplename)
MDA_normal_merge = MDA_all_samples_map %>% filter(TissueCode=="NB") %>% select(Original_ID, Barcode) 
NS_merge = ns_sample_map %>% select(Original_ID, Barcode) 

# Combine all 6 cohorts into one master table. 
LifeHistory_barcodes = bind_rows(HongKong_merge, MDA_tumor_merge, TCGA_merge, HGBM_merge, RVJDG_merge, MDA_normal_merge, NS_merge)

## Addition of random string to use as a shortened identifier.
# Create a random string identifier for each SAMPLE in the Life History study.
set.seed(1)
uuid_LifeHistory = as.data.frame(stri_rand_strings(dim(LifeHistory_barcodes)[1], 6, "[A-Z0-9]"))
colnames(uuid_LifeHistory) = "uuid"

# Sanity check: make sure each is unique.
ifelse(n_distinct(uuid_LifeHistory)==dim(LifeHistory_barcodes)[1], 
       message("UUIDs are unique"), 
       message("WARNING! Not unique"))

# Write final combined dataset set for a data freeze.
Life_history_all_samples = bind_cols(uuid_LifeHistory, LifeHistory_barcodes)
Final_Life_History = Life_history_all_samples %>% 
  separate(Barcode, c("ProjectID", "TissueSourceSite", "SubjectID", "TissueType"), remove = FALSE) %>% 
  mutate(portion = "01",
         analyte_type = "WGS") %>% 
  rename(legacy_sample_id = Original_ID, sample_id = Barcode) %>% 
  select(sample_id, uuid, legacy_sample_id, portion, analyte_type)


# Write file to be uploaded to GLASS-WG github page.
write.table(Final_Life_History, file='data/ref/glass_wg_sample_mapping_table.txt', quote=FALSE, sep='\t', row.names = F)

### Create dictionaries for different tissue centers. ####
# Subset list of TissueSourceSites from TCGA plus add new rows for other datasets.
TSS_info = read_tsv(tss_path)

# Retrieve list of TissueSourceSites for future use.
TCGA_TSS = unique(TCGA_sample_map$TSS)

# Collapse space into "." separated. 
colnames(TSS_info) = gsub( " ", ".", colnames(TSS_info))

# Create more formal definitions.
TSS_GLASS = data.frame(TSS.Code   = c("HK", "MD", "HF", "NS"),
                Source.Site  = c("Chinese University of Hong Kong", "MD Anderson Cancer Center", "Henry Ford Hospital", "Northern Sydney Cancer Centre"),
                Study.Name = c("Lower Grade Glioma / Glioblastoma", "Lower Grade Glioma / Glioblastoma", "Glioblastoma", "Metastatic Gliosarcoma"),
                BCR      = rep("GLASS", 4))
# Only those TCGA centers that are present in the NB, TP, R1/R2 samples.
TSS_TCGA_filtered = TSS_info %>% 
  filter(TSS.Code %in% TCGA_TSS)

# Create new tissue source site dictionary.
TissuSourceSite_dict = bind_rows(TSS_TCGA_filtered, TSS_GLASS)
write.table(TissuSourceSite_dict, file='data/ref/life-history-tissue-source-sites.txt', quote=FALSE, sep='\t', row.names = F)

# Create sample type dictionary. 
SampleType_dict = data.frame(Code   = c("NB", "TP", "R1", "R2", "R3", "M1"),
                       Definition  = c("Blood Derived Normal", "Primary Tumor", "First Recurrence Tumor", "Second Recurrence Tumor", "Third Reccurence Tumor", "Bone Metastasis"))
write.table(SampleType_dict, file='data/ref/life-history-sample-type-dict.txt', quote=FALSE, sep='\t', row.names = F)