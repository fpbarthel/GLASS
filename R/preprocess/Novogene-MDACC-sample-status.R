#######################################################
# Determine which samples have been sequenced by Novogene at this time.
# Date: 2018.05.29
#######################################################

# This file was provided by Katio Shao @ Novogene regarding sequencing of 107 MDACC samples. 
NVGN_seq_path = "/Users/johnsk/Documents/Life-History/GLASS-WG/data/sequencing-information/MDACC/C202SC18030593_qc_20180423/src/tables/qc.summary.xlsx"
# Novogene shipped MDACC 156 samples to JAX for library preparation. For JAX GT QC reasons, only 121 have been sequenced.
NVGN_JAX_shipment_path = "/Users/johnsk/Documents/Life-History/GLASS-WG/data/sequencing-information/MDACC/Novogene-sample-return-C202SC17090320.xlsx"
# Samples submitted to GT from Novogene, accessed through LIMS submission. Hoon submitted them under 18-verhaak-005.
GT_submitted_path = "/Users/johnsk/Documents/Life-History/GLASS-WG/data/sequencing-information/MDACC/JAX_GT_LIMS_SUBMISSION-SIF-C202SC17090320-Davis-US-MDACC-150-human-WGS-WOBI-20180205-20180206222105S-Submitted.xlsx"
# Master sheet provided Kristin (Program coordinator of Neuro-oncology at MDACC). Contains sample information for sequenced samples (n=39 subjects).
master_log_path = "/Users/johnsk/Documents/Life-History/GLASS-WG/data/clinical-data/MDACC/MDA-Clinical-Dataset/Master Log for WGS_26Apr2018wData.xlsx"

#######################################################

library(tidyr)
library(dplyr)
library(openxlsx)

#######################################################
#### Master ###
# Kristin kindly pointed me to the sequencing identifier link to these samples. 
MDA_master = readWorkbook(master_log_path, sheet = 2, startRow = 1, colNames = TRUE)

# "mdacc_sx_acc" or "MRN" can link to other files.
# "Jax.Lib.Prep.Customer.Sample.Name" contains the file name submitted to Novogene (Tumor samples only).
n_distinct(substring(MDA_master$Jax.Lib.Prep.Customer.Sample.Name, 1, 17)) # All 39 samples are present in the master file.
n_distinct(MDA_master$mdacc_sx_acc) # 81 distinct tumor samples are present in the master file.

# Not all of the samples are separated by the same characters. Revise for consistency.
MDA_master$LinkerID = gsub("____", "__", MDA_master$Jax.Lib.Prep.Customer.Sample.Name)

# Isolate parts of long identifier to easily distinguish a subject and its samples.
MDA_master_df = MDA_master %>% 
  separate(LinkerID, c("FastqID", "TubeID"), sep = "---", remove = FALSE) %>% 
  separate(TubeID, c("ID", "SampleID", "ShortenedID"), sep="__", remove= FALSE)  


#### Novogene Sequencing Batch 2 ####
# We are trying to figure how many triplets have already been sequenced by Novogene for MDA samples.
# From 107 samples made available (not delivered) to us on April 26.
Novogene_sequenced = readWorkbook(NVGN_seq_path, sheet = 1, startRow = 1, colNames = TRUE)

# Extract unique subject ID.
Novogene_sequenced$SubjectID = sapply(strsplit(Novogene_sequenced$Sample.name, "_"), "[", 4)

# All of the blood samples (minus Dr. Yung sample) were sequenced in this batch. Both "blood" and "Blood" used in the name.
length(grep("lood", unique(Novogene_sequenced$Sample.name)))

# Flag the samples where there are not at least three samples.
IDs_missing_data = Novogene_sequenced %>% 
  filter(!duplicated(Sample.name)) %>% 
  group_by(SubjectID) %>% 
  summarise(count = n()) %>% 
  filter(count < 3)

# Isolate subjects without three unique samples.
Incomplete_data_samples = Novogene_sequenced %>% 
  filter(SubjectID%in%IDs_missing_data$SubjectID)

# There's a sample from Alfred Yung without other accompanying tissues.
# This pushes the number of unique subjects to 40, but should be 39.
n_distinct(Novogene_sequenced$SubjectID)

# The samples that were sequenced at an earlier date (2018.04.05)
# contain the third sample for each of these subjects plus a fourth sample
# for GLASS-01-0027. File path:
# /verhaak-temp/kimh/kimh_JAX_SEQ/MDACC_NOVOGENE_BATCH1/RAW/20180405/128.120.88.242/C202SC18030593/raw_data

# So what was originally shipped back to JAX? Is there a discrepancy between what was shipped versus what JAX processed?
# There were 232 samples (DNA isolation) totalling 156 tissue samples in 51 subjects.
Novogene_shipped = readWorkbook(NVGN_JAX_shipment_path, sheet = 1, startRow = 1, colNames = TRUE)
Novogene_shipped$Remove.Space = gsub("ID ", "ID", Novogene_shipped$Sample.name)
Novogene_shipped$Add.Sep = gsub(" ", "_", Novogene_shipped$Remove.Space)

# Reformat Master MDA sample identifier to match that of the shipped sample.
MDA_master_df$ReformattedID = gsub("__", "_", MDA_master_df$TubeID)

# Extract samples that do not match between what was shipped and what was in the master sheet.
Samples_Absent_Master = dplyr::setdiff(Novogene_shipped$Add.Sep, MDA_master_df$ReformattedID)
# There were no blood samples in the master list. Count and then remove those.
to_drop = grep("lood",  Samples_Absent_Master)

# These are the tumor samples for which we have samples, but have not been prepared for sequencing. 
Tumor_Samples_Shipped_Master_Diff = Samples_Absent_Master[-to_drop]

# Files in LIMS that were submitted to JAX Genome Technologies.
Samples_submitted_GT = readWorkbook(GT_submitted_path, sheet = 1, startRow = 1, colNames = TRUE)

# Not all of the samples are separated by the same characters. Revise for consistency.
Samples_submitted_GT$LinkerID = gsub("____", "__", Samples_submitted_GT$Customer.Sample.Name)

# Isolate parts of long identifier to easily distinguish a subject and its samples. Blood sample throws Warning, ignore.
GT_samples_master = Samples_submitted_GT %>% 
  separate(LinkerID, c("FastqID", "TubeID"), sep = "---", remove = FALSE) %>% 
  separate(TubeID, c("ID", "SampleID", "ShortenedID"), sep="__", remove= FALSE)  

# Prepare list of sample names that were not sequenced.
total_sample = GT_samples_master %>% 
  group_by(ID) %>% 
  summarise(count = n())

# Quickly visualize the number of biological samples per subject.
hist(total_sample$count)

# There's no real difference between these two samples sets as there was
# "ID185_Blood_958398" and "ID185_blood_958398". Perhaps copy error.
# Correct some inconsistencies with TubeID.
GT_samples_master$ReformattedID = gsub("__", "_", GT_samples_master$TubeID)
dplyr::setdiff(Novogene_shipped$Add.Sep, GT_samples_master$ReformattedID)

# The "missing" samples failed one of two QC steps, either insufficient material for part of the sample trio
# or low-input FFPE for one of the tumor tissues. The result was the removal of 4 subjects based on a failed library preparation
# and 11 subjects who were excluded on the basis of low-input samples where a library preparation was not attempted.

