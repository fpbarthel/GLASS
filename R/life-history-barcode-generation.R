#######################################################
# To create sample barcodes for glioma life history
# Date: 2018.05.30
# Author: Kevin J
#######################################################

############################################################
###   WARNING: Incomplete. Need MDACC Sample Type Info   ###
############################################################

setwd("/Users/johnsk/Documents/Life-History/GLASS-WG/")

HongKong_path = "data/sequencing-information/HK/hong-kong-sample-maps.txt"
MDA_master_path = "data/clinical-data/MDACC/MDA-Clinical-Dataset/Master Log for WGS_26Apr2018wData.xlsx"
TCGA_path = "data/clinical-data/TCGA/LGG-GBM-samples.tsv"
RV_JDG_path = "data/clinical-data/Roel-JDG/Roel_JDG-original_bam_map.txt"
hGBM_path = "data/clinical-data/HF/HGBM-original_bam_map.txt"
TSS_path = "data/sequencing-information/tissueSourceSite.txt"

#######################################################

library(tidyverse)
library(openxlsx)
library(stringr)
library(stringi)

#######################################################
# Inspect how each cohort currently barcodes each file. 
HK_df = read.table(HongKong_path, header=T)
MDA_df = readWorkbook(MDA_master_path, sheet = 2, startRow = 1, colNames = TRUE)
TCGA_df = read_tsv(TCGA_path)
RV_JDG_df = read_tsv(RV_JDG_path)
hGBM_df = read_tsv(hGBM_path)

# The sample-specific barcode should include the following:
# Project (TCGA/GLSS), Tissue Source Site (##), Subject Code (4), Sample Code (8), Tissue Code (01/02/10)
## Hong Kong:
glimpse(HK_df) # GLSS - HK - 000# - 

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

# MD Anderson:
glimpse(MDA_df)

# Generate a sample type column.
MDA_df$Sample_Type <- "NA"
MDA_df$Sample_Type[MDA_df$`Order.of.Sx-1`=="first"] <- "first_surgery"
MDA_df$Sample_Type[MDA_df$`Order.of.Sx-2`=="second"] <- "second_surgery"
MDA_df$Sample_Type[MDA_df$`Order.of.Sx-3`=="third"] <- "third_surgery"
MDA_df$Sample_Type[MDA_df$`Order.of.Sx-4`=="fourth"] <- "fourth_surgery"





# Simplify the SAMPLE identifier to 8-digits.
MDA_df$New_Sample_ID = gsub("S", "", (gsub("-", "", MDA_df$mdacc_sx_acc)))
# Extract the 4-digit SUBJECT identifier.
MDA_df$SubjectID = sapply(strsplit(MDA_df$Jax.Lib.Prep.Customer.Sample.Name, "-"), "[", 3)
# We still don't have any information about which samples are primary, recurrent, or blood.

# Create a new MD Anderson to map the old to new sample names.
MDA_sample_map = MDA_df %>% 
  mutate(SeqID = "GLSS") %>% 
  mutate(TSS = "MD") %>% 
  rename(SubjectCode = SubjectID) %>% 
  rename(SampleCode = New_Sample_ID) %>%
  mutate(TissueCode = "??") %>%
  unite(Barcode, c(SeqID, TSS, SubjectCode, TissueCode), sep = "-", remove = FALSE)

# TCGA:
glimpse(TCGA_df) # TCGA - 0# - [:alnum:]{4} -
# TCGA tissue source sites.
unique(sapply(strsplit(TCGA_df$aliquot_id, "-"), "[", 2))

# Create a new TCGA to map the old to new sample names.
TCGA_sample_map = TCGA_df %>% 
  mutate(SeqID = "TCGA") %>% 
  mutate(TSS = sapply(strsplit(TCGA_df$aliquot_id, "-"), "[", 2)) %>% 
  mutate(SubjectCode = substring(aliquot_id, 9, 12)) %>% 
  mutate(TissueCode = recode_factor(substring(aliquot_id, 14, 16), "01A" = "TP", "01B" = "TP", "02A" = "R1", "02B" = "R2", "10A" = "NB", "10B" = "NB", "10D" = "NB")) %>% 
unite(Barcode, c(SeqID, TSS, SubjectCode, TissueCode), sep = "-", remove = FALSE)
  
# RV-JDG:
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

# HGBM:
glimpse(hGBM_df)

# Create a new RVJDG sample map to connect with the old sample names.
HGBM_sample_map = hGBM_df %>% 
  filter(samplename!="HF-3177-10-01D") %>% 
  mutate(SeqID = "GLSS") %>% 
  mutate(TSS = "HF") %>% 
  mutate(SubjectCode = substring(PatientID, 4, 7)) %>% 
  mutate(TissueCode = recode_factor(SampleType2, "P" = "TP", "R" = "R1", "N" = "NB")) %>% 
  unite(Barcode, c(SeqID, TSS, SubjectCode, TissueCode), sep = "-", remove = TRUE)

#############################################
# Combine 5 datasets into master table with each linker
#############################################
HongKong_merge = HongKong_sample_map %>% select(Verhaak_Sample_ID, Barcode) %>% rename(Original_ID = Verhaak_Sample_ID)
MDA_merge = MDA_sample_map %>%  select(Jax.Lib.Prep.Customer.Sample.Name, Barcode)  %>% rename(Original_ID = Jax.Lib.Prep.Customer.Sample.Name)
TCGA_merge = TCGA_sample_map %>%  select(aliquot_id, Barcode)  %>% rename(Original_ID = aliquot_id)
HGBM_merge = HGBM_sample_map %>%  select(samplename, Barcode)  %>% rename(Original_ID = samplename)
RVJDG_merge = RVJDG_sample_map %>%  select(samplename, Barcode)  %>% rename(Original_ID = samplename)

# Combine all 5 cohorts into one master table. 
LifeHistory_barcodes = bind_rows(HongKong_merge, MDA_merge, TCGA_merge, HGBM_merge, RVJDG_merge)

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
  separate(Barcode, c("ProjectID", "TissueSourceSite", "SubjectID", "TissueType"), remove = FALSE)

# Write file to be uploaded to GLASS-WG github page.
write.table(Final_Life_History, file='data/sequencing-information/master_life_history_uniform_naming_incomplete.txt', quote=FALSE, sep='\t', row.names = F)


### Create dictionaries for different tissue centers.
# Subset list of TissueSourceSites from TCGA plus add new rows for other datasets.
TSS_info = read_tsv(TSS_path)

# Retrieve list of TissueSourceSites for future use.
TCGA_TSS = unique(TCGA_sample_map$TSS)

# Collapse space into "." separated. 
colnames(TSS_info) = gsub( " ", ".", colnames(TSS_info))

# Create more formal definitions.
TSS_GLASS = data.frame(TSS.Code   = c("HK", "MD", "HF"),
                Source.Site  = c("Chinese University of Hong Kong", "MD Anderson Cancer Center", "Henry Ford Hospital"),
                Study.Name = c("Lower Grade Glioma / Glioblastoma", "Lower Grade Glioma / Glioblastoma", "Glioblastoma"),
                BCR      = rep("GLASS", 3))
# Only those TCGA centers that are present in the NB, TP, R1/R2 samples.
TSS_TCGA_filtered = TSS_info %>% 
  filter(TSS.Code %in% TCGA_TSS)

# Create new tissue source site dictionary.
TissuSourceSite_dict = bind_rows(TSS_TCGA_filtered, TSS_GLASS)
write.table(TissuSourceSite_dict, file='data/sequencing-information/life-history-tissue-source-sites.txt', quote=FALSE, sep='\t', row.names = F)

# Create sample type dictionary. 
SampleType_dict = data.frame(Code   = c("NB", "TP", "R1", "R2"),
                       Definition  = c("Blood Derived Normal", "Primary Tumor", "First Recurrence Tumor", "Second Recurrence Tumor"))
write.table(SampleType_dict, file='data/sequencing-information/life-history-sample-type-dict.txt', quote=FALSE, sep='\t', row.names = F)

