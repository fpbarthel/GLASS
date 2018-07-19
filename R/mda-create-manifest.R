#######################################################
# Create manifest for MD Anderson samples (GLASS)
# Date: 2018.07.13
# Author: Kevin J.
#######################################################

# Local directory for github repo.
mybasedir = "/Users/johnsk/Documents/Life-History/GLASS-WG"
setwd(mybasedir)

# Files with information about fastq information and barcodes.
MDA_batch1_file_path = "data/sequencing-information/MDACC/file_list_C202SC18030593_batch1_n14_20180405.tsv"
MDA_batch2_file_path = "data/sequencing-information/MDACC/file_list_C202SC18030593_batch2_n94_20180603.tsv"
mda_master_path = "data/clinical-data/MDACC/MDA-Clinical-Dataset/Master Log for WGS_sampletype.20180630.xlsx"
novogene_sample_path = "data/sequencing-information/MDACC/Novogene_SIF_14.xlsx"

# Create extensions for samples.
json_ext = "json"
text_ext = "tsv"

# Create output files for each metadata set.
cases_file      = "data/manifest/mdacc/cases"
samples_file    = "data/manifest/mdacc/samples"
aliquots_file   = "data/manifest/mdacc/aliquots"
readgroups_file = "data/manifest/mdacc/readgroups"
files_file      = "data/manifest/mdacc/files"
pairs_file      = "data/manifest/mdacc/pairs"

#######################################################

library(tidyverse)
library(openxlsx)
library(rjson)
library(jsonlite)
library(listviewer)
library(stringi)
library(stringr)

#######################################################
# We need to generate the following fields required by the SNV snakemake pipeline:
# aliquots, files, cases, samples, pairs, and readgroups.


#### NEED TO GENERATE NEW CODES FOR MDACC NVGN SAMPLES ####

### Aliquot ####
# barcode_sheet = read.delim("data/sequencing-information/master_life_history_uniform_naming_incomplete.txt", as.is=T)
# aliquot_sheet = master_sheet %>% select(aliquot_uuid = uuid, sample_id = Barcode) %>%
#  mutate(aliquot_id = sprintf("%s-%s", sample_id, aliquot_uuid),
#         analyte_type = "DNA",
#         analysis_type = "WGS",
#         portion = 1) %>%
#  filter(grepl("HK", sample_id))

# Aliquot file to be written.
# aliquots = aliquot_sheet %>% select(sample_id, aliquot_uuid, aliquot_id, portion, analyte_type, analysis_type) %>% distinct()


### MDA (Novogene cohort) barcode generation ####
# Novogene has processed 114 samples. In our master sheet normal blood samples were not provided. We need them for barcodes.
# Inspect both batch1 and batch2 sequencing data for 
mda_batch1_df <- read.table(MDA_batch1_file_path, col.names="file_path", stringsAsFactors = F)
mda_batch1_df$filename = sapply(strsplit(mda_batch1_df$file_path, "/"), "[[", 8)
mda_batch2_df <- read.table(MDA_batch2_file_path, col.names="file_path", stringsAsFactors = F)
mda_batch2_df$filename = sapply(strsplit(mda_batch2_df$file_path, "/"), "[[", 3)

# Replace the placeholder for pwd "./" from unix cmd: "find -name ".fq.gz" in parent directory of fastqs. 
mda_batch1_df$file_path <- gsub("\\./", "/verhaak-temp/kimh/kimh_JAX_SEQ/MDACC_NOVOGENE_BATCH1/", mda_batch1_df$file_path)
mda_batch2_df$file_path <- gsub("^", "/verhaak-temp/GLASS/life_history/l1_seq/mdacc/", mda_batch2_df$file_path)
mda_batch1_df$old_sample_id <- gsub( "_USPD.*$", "", mda_batch1_df$filename)
mda_batch2_df$old_sample_id <- gsub( "_USPD.*$", "", mda_batch2_df$filename)

# Combine both sequencing data sets.
mda_df <- bind_rows(mda_batch1_df, mda_batch2_df)

# Retrieve only the blood samples ("10D"). Blood samples not contained in master sheet.
normal_blood_samples = mda_df[grep("[b|B]lood", mda_df$filename), ]

# Sequencing results stripped the name, changing it from "GLASS_01_00" to "G01_" in the output. Reverse that modification.
normal_blood_samples$revised_id <- gsub("G01_", "GLASS_01_00", normal_blood_samples$old_sample_id)

# Retrieve the subject ID in same format as tumor samples.
normal_blood_samples$SubjectID = sapply(strsplit(normal_blood_samples$revised_id, "_"), "[", 3)

# Create a barcode for available blood samples (n=33).
mda_normal_blood_map <- normal_blood_samples %>% 
  mutate(SeqID = "GLSS") %>% 
  mutate(TSS = "MD") %>% 
  rename(SubjectCode = SubjectID) %>% 
  mutate(TissueCode = "NB") %>% 
  mutate(Original_ID = revised_id) %>%   
  unite(Barcode, c(SeqID, TSS, SubjectCode, TissueCode), sep = "-", remove = FALSE) %>% 
  select(Original_ID, Barcode, SeqID, TSS, SubjectCode, TissueCode) %>% 
  distinct()
  

### Use master sheet containing tumor sample information, necessary to build barcodes.####
# Kristin Alfaro-Munoz kindly pointed me to the sequencing identifier link to these samples. 
mda_master = readWorkbook(mda_master_path, sheet = 2, startRow = 1, colNames = TRUE)

novogene_linker = readWorkbook(novogene_sample_path, sheet = 1, startRow = 18, colNames = TRUE)
novogene_linker_unique <- novogene_linker[1:121, ]

sum(novogene_linker_unique$'*SampleName'%in%mda_master$Jax.Lib.Prep.Customer.Sample.Name)

# Extract the 4-digit SUBJECT identifier.
mda_master$SubjectID = sapply(strsplit(mda_master$Jax.Lib.Prep.Customer.Sample.Name, "-"), "[", 3)

# Create a new MD Anderson to map the old to new sample names.
mda_tumor_sample_map = mda_master %>% 
  mutate(SeqID = "GLSS") %>% 
  mutate(TSS = "MD") %>% 
  rename(SubjectCode = SubjectID) %>% 
  rename(TissueCode = sample_type) %>% 
  mutate(Original_ID = Jax.Lib.Prep.Customer.Sample.Name) %>% 
  unite(Barcode, c(SeqID, TSS, SubjectCode, TissueCode), sep = "-", remove = FALSE) %>% 
  select(Original_ID, Barcode, SeqID, TSS, SubjectCode, TissueCode)

mda_all_samples_df <- bind_rows(mda_normal_blood_map, mda_tumor_sample_map)

## Addition of random string to use as a shortened identifier.
# Create a random string identifier for each SAMPLE in the Life History study.
set.seed(39)
mda_life_history_barcodes <- mda_all_samples_df %>% 
  mutate(uuid = stri_rand_strings(dim(mda_all_samples_df)[1], 6, "[A-Z0-9]")) %>% 
  rename(TissueSourceSite = TSS, ProjectID = SeqID, TissueType = TissueCode, SubjectID = SubjectCode) %>% 
  select(uuid, Original_ID, Barcode, ProjectID, TissueSourceSite, SubjectID, TissueType) 

# Check to make sure there is no overlap with previous barcodes.
# sum(barcode_sheet$uuid%in%mda_life_history_barcodes) # No matches. I was concerned because I already generated the uuid
# for the previous datasets that Floris used in his pipeline.
# Write file to be uploaded to GLASS-WG github page. ***Batch1 and Batch2 only. Missing samples. ***
# write.table(mda_life_history_barcodes, file='data/sequencing-information/mda_life_history_uniform_naming.txt', quote=FALSE, sep='\t', row.names = F)


### Aliquot ####
mda_barcode_sheet = read.delim("data/sequencing-information/mda_life_history_uniform_naming.txt", as.is=T)
aliquot_sheet = mda_barcode_sheet %>% select(aliquot_uuid = uuid, sample_id = Barcode) %>%
  mutate(aliquot_id = sprintf("%s-%s", sample_id, aliquot_uuid),
         analyte_type = "DNA",
         analysis_type = "WGS",
         portion = 1) 

# Aliquot file to be written.
aliquots = aliquot_sheet %>% select(sample_id, aliquot_uuid, aliquot_id, portion, analyte_type, analysis_type) %>% distinct()

# write(jsonlite::toJSON(aliquots, pretty = T), file = sprintf("%s.%s", aliquots_file, json_ext))
# write.table(aliquots, file = sprintf("%s.%s", aliquots_file, text_ext), sep="\t", row.names = F, col.names = T, quote = F)


### Files ####
# Generate *file* tsv containing: aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format.
# From above ^: mda_df = bind_rows(mda_batch1_df, mda_batch2_df)

# Retrieve the library, flowcell, and lane id from the filename.
mda_df_meta  = mda_df %>% 
  mutate(library_id = sub(".*_ *(.*?) *_H.*", "\\1", filename), 
         flowcell_id = substr(filename, nchar(filename)-19, nchar(filename)-12),
         lane_id = substr(filename, nchar(filename)-8, nchar(filename)-8))

# Create a new identifier on which to group mate pairs onto the same line (i.e., R1 and R2).
mda_df_meta$read_group = paste(mda_df_meta$library_id, mda_df_meta$flowcell_id, mda_df_meta$lane_id, sep = '-')

# Comma separated file_paths and file_names.
merged_mda_files = mda_df_meta %>% 
  group_by(read_group) %>% 
  mutate(file_path = paste(file_path, collapse=","), 
         file_name = paste(filename, collapse=",")) %>% 
  select(-filename) %>% 
  distinct()

# Combine with study center provided covariate sheet.
# Not all of the samples are separated by the same characters. Revise for consistency.
# 1. Fixing "GLASS_01_00" >> "G01_".
merged_mda_files$revised_id <- gsub("G01_", "GLASS_01_00", merged_mda_files$old_sample_id)

test_map_df = merged_mda_files %>% 
  inner_join(novogene_linker_unique, by = c("library_id" = "Novogene.ID")) %>%  
  ungroup()


# Create a temporary file to determine overlap of IDs.
tmp1 = unique(merged_mda_files$revised_id)

# Seq batch ID:  GLASS_01_0003_99D_131077_S_11 ***However, all samples don't follow the same pattern or edits
# Master ID:     GLASS-01-0003-99D---131077__S-11-087625__495670
mda_master$linker_id = gsub( "____", "__", mda_master$Jax.Lib.Prep.Customer.Sample.Name)
mda_master$revised_id = gsub( "---", "-", mda_master$linker_id)

tmp = mda_master %>% 
  separate(revised_id, c("SUBJECT_ID", "SAMPLE_ID", "TUBE_LABEL"), sep = "__") %>%  
  separate(SAMPLE_ID, c("S_code", "N_code", "Acc_code"), sep = "-") %>% 
  mutate(SUBJECT_ID = gsub("-", "_", SUBJECT_ID)) %>% 
  mutate(TO_MERGE = paste(SUBJECT_ID, S_code, N_code, sep = "_"))

# Create temporary file on which you would like to merge.
tmp2 = tmp$TO_MERGE

# Remove the blood samples from tmp1.
tumor_no_blood = tmp1[-grep("(B|b)lood", tmp1)]
blood_ids = tmp1[grep("(B|b)lood", tmp1)]
sum(tumor_no_blood[15:75]%in%tmp2)

# List of sample names as they appear in the fastq names.
sequenced_samples <- c(tumor_no_blood, blood_ids)
write.table(sequenced_samples, file='/Users/johnsk/Documents/Life-History/life-history-batch1-batch2-sequenced-sample-names.txt', quote=FALSE, sep='\t', row.names = F)
write.table(blood_ids, file='/Users/johnsk/Documents/Life-History/life-history-batch1-batch2-sequenced-blood-names.txt', quote=FALSE, sep='\t', row.names = F)

# Eventually, we want this data.frame.
test_map_df = merged_mda_files %>% 
  inner_join(tmp, by = c("revised_id" = "TO_MERGE")) %>%  
  ungroup()


# Generate uuids for each of the hong kong files.
# Example from TCGA: 24c6f54a-e7a2-4148-8335-045e3c74096e
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
                                 aliquot_id = sprintf("%s-%s", Barcode, uuid), 
                                 file_format = "FQ")

# Order needs to be: # aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format.
files = hk_map_df %>% select(aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format) %>% distinct()




# 1. Batch 1 samples do not possess sufficient data to merge directly with the master sheet. We'll work by a process of
# elimination. Merge with the batch2 samples then only merge with the remaining samples in batch1.
# 2. The sample names in the sequencing files and the master don't match up. Missing leading zeros for some, and very different
# character formats for other samples.


# Convert to "-" separator and pad single digits with one zero.
sprintf("%02d", gsub(".*_", "", mda_batch2_df$old_sample_id))

# mda_df <- bind_rows(mda_batch1_df, mda_batch2_df)
# mda_df$old_sample_id <- gsub( "_USPD.*$", "", mda_df$filename)
# Create a new data.frame with just the old sample_id.
# mda_all_samples <- as.data.frame(unique(mda_df$old_sample_id))
# colnames(mda_all_samples) <- "mda_ids"
# mda_all_samples$mda_ids <- as.character(mda_all_samples$mda_ids)
# mda_all_samples$mda_ids[1:14] <- sapply(strsplit(mda_all_samples$mda_ids[1:14], "_"), "[[", 4)
# Create consistent 17 character string for the "subject" indentifier with tumor (99D) and normal (10D)
# designation.
# mda_all_samples$mda_ids[15:108] <-gsub("G01_", "GLASS_01_00", mda_all_samples$mda_ids[15:108])
# mda_all_samples$subject_id = substring(mda_all_samples$mda_ids, 1, 17)
# Create tumor  
# mda_all_samples$sample_type <- "NA"
# mda_all_samples$sample_type <- sapply(strsplit(mda_all_samples$mda_ids, "_"), "[[", 4)



# Not all of the samples are separated by the same characters. Revise for consistency.
mda_master$linker_id = gsub( "---.*$", "", mda_master$Jax.Lib.Prep.Customer.Sample.Name)

# Generate a sample type column.
mda_master$surgery_type <- "NA"
mda_master$surgery_type[mda_master$`Order.of.Sx-1`=="first"] <- "1"
mda_master$surgery_type[mda_master$`Order.of.Sx-2`=="second"] <- "2"
mda_master$surgery_type[mda_master$`Order.of.Sx-3`=="third"] <- "3"
mda_master$surgery_type[mda_master$`Order.of.Sx-4`=="fourth"] <- "4"

# Based on consensus from the life-history group samples were labelled as "TP", "R1", "R2", or "R3"
# depending on whether the subject had missing first surgery samples.
table(mda_master$sample_type)

# Combine datasets based on linker files. 
mda_batch2_merge = mda_all_samples[15:108, ]



mda_map_df = mda_master %>% 
  inner_join(mda_all_samples, by = c("linker_id" = "mda_ids"))
  




# Need to generate a complete barcode as the MD Anderson normals weren't included in the original.
# Create a new MD Anderson to map the old to new sample names.
MDA_sample_map = mda_all_samples_id %>% 
  mutate(SeqID = "GLSS") %>% 
  mutate(TSS = "MD") %>% 
  rename(SubjectCode = SubjectID) %>% 
  rename(SampleCode = New_Sample_ID) %>%
  mutate(TissueCode = recode_factor(sample_type), "99D" = "TU", "10D" = "NB") %>%
  unite(Barcode, c(SeqID, TSS, SubjectCode, TissueCode), sep = "-", remove = FALSE)



# Normal blood samples in this dataset defined by "B|blood" and 10D (Tumor = 99D).
normal_blood_samples = mda_df[grep("[b|B]lood", mda_df$filename), ]
putative_normal_blood_samples = mda_df[grep("10D", mda_df$filename), ]
putative_tumor_samples = mda_df[grep("99D", mda_df$filename), ]

# Retrieve the original sample name to map to study center provided covariate sheet.
hk_df_meta = mda_df %>% 
  mutate(verhaak_sample_id = sub(".*/WG_ *(.*?) *_USPD.*", "\\1", filenames),
         library_id = sub(".*_ *(.*?) *_H.*", "\\1", filenames), 
         flowcell_id = substr(filenames, nchar(filenames)-19, nchar(filenames)-12),
         lane_id = substr(filenames, nchar(filenames)-8, nchar(filenames)-8))

# Create a new identifier on which to group mate pairs onto the same line (i.e., R1 and R2).
hk_df_meta$read_group = paste(hk_df_meta$library_id, hk_df_meta$flowcell_id, hk_df_meta$lane_id, sep = '-')

# Retrieve the file_name from the file_path. 
hk_df_meta$file_name_single = sapply(strsplit(hk_df_meta$filenames, "/"), "[[", 8)

# Comma separated file_paths and file_names.
merged_hk_files = hk_df_meta %>% 
  group_by(read_group) %>% 
  mutate(file_path = paste(filenames, collapse=","),
         file_name = paste(file_name_single, collapse=",")) %>% 
  select(-filenames, -file_name_single) %>% 
  distinct()

# Combine with study center provided covariate sheet.
master_sheet$Original_ID = gsub("-","_", master_sheet$Original_ID)
hk_map_df = merged_hk_files %>% 
  inner_join(master_sheet, by = c("verhaak_sample_id" = "Original_ID")) %>%  
  ungroup()

# Generate uuids for each of the hong kong files.
# Example from TCGA: 24c6f54a-e7a2-4148-8335-045e3c74096e
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
                                 aliquot_id = sprintf("%s-%s", Barcode, uuid), 
                                 file_format = "FQ")

# Order needs to be: # aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format.
files = hk_map_df %>% select(aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format) %>% distinct()




### Cases ####
hk_map_df = hk_map_df %>% 
  mutate(case_id = substring(Barcode, 1, 12), 
         project_id = "GLSS-HK")
# Select only those relevant fields.
cases = hk_map_df %>% select(case_id, project_id)


