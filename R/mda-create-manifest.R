#######################################################
# Create manifest for MD Anderson samples (n = 120, GLASS)
# Date: 2018.08.07
# Author: Kevin J.
#######################################################
# Local directory for github repo.
mybasedir = "/Users/johnsk/Documents/Life-History/GLASS-WG"
setwd(mybasedir)

# Files with information about fastq information and barcodes.
# Actual fastq files are stored: /fastscratch/johnsk/GLASS-WG/mdacc/ or on tier2 (long-term).
MDA_batch1_file_path = "data/sequencing-information/MDACC/file_list_C202SC18030593_batch1_n14_20180405.tsv"
MDA_batch2_file_path = "data/sequencing-information/MDACC/file_list_C202SC18030593_batch2_n94_20180603.tsv"
MDA_batch3_file_path = "data/sequencing-information/MDACC/file_list_C202SC18030593_batch2_n13_20180716.tsv"

# Clinical dataset priovided by Kristin Alfaro-Munoz at MD Anderson.
# Modified to generate a sample_type (TP, R1, R2, R3) using the clinical information. See github issue #16.
mda_master_path = "data/clinical-data/MDACC/MDA-Clinical-Dataset/Master Log for WGS_sampletype.20180630.xlsx"

# 2018.07.06 Katie Shao (Novogene) provided the following sheet linking libraryID with submitted sample names in the email:
# "Re:RE: Novogene Project Report - Confirm -C202SC18030593-Davis-MDACC-156-libseq-hWGS-WOBI-NVUS2018022505[J1-C202SC18030593-1000]"
novogene_sample_path = "data/sequencing-information/MDACC/Novogene_SIF_14.xlsx"

# Completed life-history barcodes.
life_history_barcodes = "data/sequencing-information/master_life_history_uniform_naming_complete.txt"

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

### MDA (Novogene cohort) barcode (re)generation ####
# Novogene has processed 120 samples. In our master sheet normal blood samples were not provided. We need them for barcodes.
# Inspect both batch1 and batch2 sequencing data for all blood/tumor data.
mda_batch1_df <- read.table(MDA_batch1_file_path, col.names="file_path", stringsAsFactors = F)
mda_batch1_df$filename = sapply(strsplit(mda_batch1_df$file_path, "/"), "[[", 8)
mda_batch2_df <- read.table(MDA_batch2_file_path, col.names="file_path", stringsAsFactors = F)
mda_batch2_df$filename = sapply(strsplit(mda_batch2_df$file_path, "/"), "[[", 3)
mda_batch3_df <- read.table(MDA_batch3_file_path, col.names="file_path", stringsAsFactors = F)
mda_batch3_df$filename = sapply(strsplit(mda_batch3_df$file_path, "/"), "[[", 3)

# Replace the placeholder for pwd "./" from bash cmd: "find -name ".fq.gz" in parent directory of fastqs. 
mda_batch1_df$file_path <- gsub("\\./RAW/20180405/128.120.88.242/C202SC18030593/raw_data/", "/fastscratch/johnsk/GLASS-WG/mdacc/C202SC18030593_batch1_n14_20180405/", mda_batch1_df$file_path)
mda_batch2_df$file_path <- gsub("^", "/fastscratch/johnsk/GLASS-WG/mdacc/", mda_batch2_df$file_path)
mda_batch3_df$file_path <- gsub("^", "/fastscratch/johnsk/GLASS-WG/mdacc/C202SC18030593_batch2_n13_20180716", mda_batch3_df$file_path)

# Create an old sample.id for these subjects to be linked. No longer *necessary*, but kept in case it's helpful.
mda_batch1_df$old_sample_id <- gsub( "_USPD.*$", "", mda_batch1_df$filename)
mda_batch2_df$old_sample_id <- gsub( "_USPD.*$", "", mda_batch2_df$filename)
mda_batch3_df$old_sample_id <- gsub( "_USPD.*$", "", mda_batch3_df$filename)

# Combine these three sequencing data sets. There should be 552 files.
mda_df <- bind_rows(mda_batch1_df, mda_batch2_df, mda_batch3_df)

# Katie Shao (Novogene provided).
novogene_linker = readWorkbook(novogene_sample_path, sheet = 1, startRow = 18, colNames = TRUE)
novogene_linker_unique <- novogene_linker[1:121, ] # 121st sample represents non-GLASS sample.

# Retrieve only the blood samples ("10D"). Blood samples not presented in master sheet.
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

# The master sheet from MDA only contains tumor samples. Sum should equal 81 (tumor samples).
sum(novogene_linker_unique$'*SampleName'%in%mda_master$Jax.Lib.Prep.Customer.Sample.Name)

# Extract the 4-digit SUBJECT identifier.
mda_master$SubjectID = sapply(strsplit(mda_master$Jax.Lib.Prep.Customer.Sample.Name, "-"), "[", 3)

# Map the old to new sample name for the tumors.
mda_tumor_sample_map = mda_master[1:81,] %>% 
  mutate(SeqID = "GLSS") %>% 
  mutate(TSS = "MD") %>% 
  rename(SubjectCode = SubjectID) %>% 
  rename(TissueCode = sample_type) %>% 
  mutate(Original_ID = Jax.Lib.Prep.Customer.Sample.Name) %>% 
  unite(Barcode, c(SeqID, TSS, SubjectCode, TissueCode), sep = "-", remove = FALSE) %>% 
  select(Original_ID, Barcode, SeqID, TSS, SubjectCode, TissueCode)

# Combine all tumor and normal samples together. Should equal 120 for this GLASS dataset.
mda_all_samples_df <- bind_rows(mda_normal_blood_map, mda_tumor_sample_map)


### Aliquot ####
mda_barcode_sheet = read.delim("data/sequencing-information/master_life_history_uniform_naming_complete.txt", as.is=T)
aliquot_sheet = mda_barcode_sheet %>% select(aliquot_uuid = uuid, sample_id = Barcode, Original_ID)  %>%
  mutate(aliquot_id = sprintf("%s-%s", sample_id, aliquot_uuid),
         analyte_type = "DNA",
         analysis_type = "WGS",
         portion = 1) %>% 
  filter(grepl("GLASS", Original_ID))

# Aliquot file to be written.
aliquots = aliquot_sheet %>% select(sample_id, aliquot_uuid, aliquot_id, portion, analyte_type, analysis_type) %>% distinct()


### Files ####
# Generate *file* tsv containing: aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format.
# From above ^: mda_df = bind_rows(mda_batch1_df, mda_batch2_df, mda_batch3_df)

# Retrieve the library, flowcell, and lane id from the filename.
mda_df_meta  = mda_df %>% 
  mutate(library_id = sub(".*_ *(.*?) *_H.*", "\\1", filename), 
         flowcell_id = substr(filename, nchar(filename)-19, nchar(filename)-12),
         lane_id = substr(filename, nchar(filename)-8, nchar(filename)-8))

# Create a new identifier on which to group mate pairs onto the same line (i.e., R1 and R2).
mda_df_meta$read_group = paste(mda_df_meta$library_id, mda_df_meta$flowcell_id, mda_df_meta$lane_id, sep = '-')

# Comma separated file_paths and file_names. "L#_1" and "L#_2" can be any order. Doesn't yet matter for Snakemake.
merged_mda_files = mda_df_meta %>% 
  group_by(read_group) %>% 
  mutate(file_path = paste(file_path, collapse=","), 
         file_name = paste(filename, collapse=",")) %>% 
  select(-filename) %>% 
  distinct()

# Not all of the samples are separated by the same characters. Revise for consistency.
# 1. Fixing  "G01_" >> GLASS_01_00".
merged_mda_files$revised_id <- gsub("G01_", "GLASS_01_00", merged_mda_files$old_sample_id)

# Combine with the new files with the sample name provided by Novogene.
# This action will remove the Yung sample.
mda_map_df = merged_mda_files %>% 
  inner_join(novogene_linker_unique, by = c("library_id" = "Novogene.ID")) %>%  
  inner_join(mda_barcode_sheet, by = c("*SampleName" = "Original_ID")) %>%  
  ungroup()

# Generate uuids for each of the mdacc files.
# Example from TCGA: 24c6f54a-e7a2-4148-8335-045e3c74096e
set.seed(39)
mda_map_df$file_uuid = paste(stri_rand_strings(dim(mda_map_df)[1], 8, "[a-z0-9]"),
                            stri_rand_strings(dim(mda_map_df)[1], 4, "[a-z0-9]"),
                            stri_rand_strings(dim(mda_map_df)[1], 4, "[a-z0-9]"),
                            stri_rand_strings(dim(mda_map_df)[1], 4, "[a-z0-9]"),
                            stri_rand_strings(dim(mda_map_df)[1], 12, "[a-z0-9]"),
                            sep = "-")

# Sanity check: make sure each is unique. Ought to be 274 at this point.
n_distinct(mda_map_df$file_uuid)

# We noticed that several fastq files are empty causing errors in Snakemake.
empty_fastqs <- read.table("data/sequencing-information/MDACC/empty_fastqs_to_remove_from_json.txt", stringsAsFactors = F, skip = 1)
empty_fastqs$fastqname <- sapply(strsplit(empty_fastqs$V1, "/"), "[[", 3)
empty_fastqs$read_group <- substr(empty_fastqs$fastqname, nchar(empty_fastqs$fastqname)-32, nchar(empty_fastqs$fastqname)-8)

# Stitched the files together as was done for the rest of the dataset.
files_to_remove = empty_fastqs %>% 
  group_by(read_group) %>% 
  mutate(filename = paste(fastqname, collapse=","),
         filename_flipped = paste(rev(fastqname), collapse=",")) %>% 
  select(-V1, -fastqname) %>% 
  distinct()

# Since, the order of the paired fastqs was not ordered I needed to use either forward or reverse orientation.
rows_to_remove <- which(mda_map_df$file_name%in%files_to_remove$filename==TRUE | mda_map_df$file_name%in%files_to_remove$filename_flipped==TRUE)

# Remove the empty files. Should be 270 at this point.
mda_map_df <- mda_map_df[-rows_to_remove, ]


# Need to record the file_size and file_md5sum for these samples.
mda_map_df = mda_map_df %>% mutate(file_size = "NA",
                                 file_md5sum = "NA",
                                 aliquot_id = sprintf("%s-%s", Barcode, uuid), 
                                 file_format = "FQ")

# Order needs to be: # aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format.
files = mda_map_df %>% select(aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format) %>% distinct()




### Cases ####
mda_map_df = mda_map_df %>% 
  mutate(case_id = substring(Barcode, 1, 12), 
         project_id = "GLSS-MD")
# Select only those relevant fields.
cases = mda_map_df %>% select(case_id, project_id) %>% 
  distinct()



### Samples ####
# Grab last two characrters of barcode.
mda_map_df$sample_type = substring(mda_map_df$Barcode, 14, 15)
# Recode variables to match Floris' fields.
samples = mda_map_df %>% select(case_id, sample_id = Barcode, legacy_sample_id = `*SampleName`, sample_type) %>% distinct()



### Pairs ####
p1 = samples %>% 
  left_join(aliquots) %>%
  select(sample_type, aliquot_id, case_id) %>%
  filter(sample_type %in% c("TP", "NB")) %>% 
  spread(sample_type, aliquot_id) %>%
  mutate(pair_id = TP) %>%
  select(case_id, pair_id, tumor_aliquot_id = TP, normal_aliquot_id = NB)

p2 = samples %>% 
  left_join(aliquots) %>%
  select(sample_type, aliquot_id, case_id) %>%
  filter(sample_type %in% c("R1", "NB")) %>%
  spread(sample_type, aliquot_id) %>%
  mutate(pair_id = R1) %>%
  select(case_id, pair_id, tumor_aliquot_id = R1, normal_aliquot_id = NB)

p3 = samples %>% 
  left_join(aliquots) %>%
  select(sample_type, aliquot_id, case_id) %>%
  filter(sample_type %in% c("R2", "NB")) %>%
  spread(sample_type, aliquot_id) %>%
  mutate(pair_id = R2) %>%
  select(case_id, pair_id, tumor_aliquot_id = R2, normal_aliquot_id = NB)

p4 = samples %>% 
  left_join(aliquots) %>%
  select(sample_type, aliquot_id, case_id) %>%
  filter(sample_type %in% c("R3", "NB")) %>%
  spread(sample_type, aliquot_id) %>%
  mutate(pair_id = R3) %>%
  select(case_id, pair_id, tumor_aliquot_id = R3, normal_aliquot_id = NB)

# Note: May need to revise once we gather additional information from Kristin @ MD Anderson.
pairs = rbind(p1, p2, p3, p4) %>% filter(complete.cases(tumor_aliquot_id, normal_aliquot_id))


### Readgroups ####
# Necessary information: file_uuid, aliquot_id, RGID, RGPL, RGPU, RGLB, RGPI, RGDT, RGSM, RGCN.
readgroup_df = mda_map_df %>% 
  mutate(RGPL = "ILLUMINA",
         RGPU = paste(substr(file_name, nchar(file_name)-19, nchar(file_name)-12), 
                      substr(file_name, nchar(file_name)-8, nchar(file_name)-8), sep="."),
         RGLB = sub(".*_ *(.*?) *_H.*", "\\1", file_name),
         RGPI = 0,
         RGDT = strftime(as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%dT%H:%M:%S%z"), 
         RGSM = Barcode,
         RGCN = "NVGN_MD",
         RGID = paste0(substring(RGPU, 1, 4), substring(RGPU, nchar(RGPU)-1, nchar(RGPU)), ""))

# Finalize readgroup information in predefined order.
readgroups = readgroup_df %>% select(file_uuid, aliquot_id, RGID, RGPL, RGPU, RGLB, RGPI, RGDT, RGSM, RGCN) %>% distinct()

### OUTPUT ####
# Output the json and .tsv files.
print(sprintf("Exporting manifest as json files for snakemake use."))
write(jsonlite::toJSON(aliquots, pretty = T), file = sprintf("%s.%s", aliquots_file, json_ext))
write(jsonlite::toJSON(files, pretty = T), file = sprintf("%s.%s", files_file, json_ext))
write(jsonlite::toJSON(cases, pretty = T), file = sprintf("%s.%s", cases_file, json_ext))
write(jsonlite::toJSON(pairs, pretty = T), file = sprintf("%s.%s", pairs_file, json_ext))
write(jsonlite::toJSON(readgroups, pretty = T), file = sprintf("%s.%s", readgroups_file, json_ext))
write(jsonlite::toJSON(samples, pretty = T), file = sprintf("%s.%s", samples_file, json_ext))

print(sprintf("Exporting manifest as tsv files for visualization ease."))
write.table(aliquots, file = sprintf("%s.%s", aliquots_file, text_ext), sep="\t", row.names = F, col.names = T, quote = F)
write.table(files, file = sprintf("%s.%s", files_file, text_ext), sep="\t", row.names = F, col.names = T, quote = F)
write.table(cases, file = sprintf("%s.%s", cases_file, text_ext), sep="\t", row.names = F, col.names = T, quote = F)
write.table(pairs, file = sprintf("%s.%s", pairs_file, text_ext), sep="\t", row.names = F, col.names = T, quote = F)
write.table(readgroups, file = sprintf("%s.%s", readgroups_file, text_ext), sep="\t", row.names = F, col.names = T, quote = F)
write.table(samples, file = sprintf("%s.%s", samples_file, text_ext), sep="\t", row.names = F, col.names = T, quote = F)


# Output RData object and timestamp along with package versions.
mysession_info <- devtools::session_info()
timetag = make.names(format(Sys.time(),"t%d_%b_%y_%H%M%S%Z"))
save.image(file.path(sprintf("R/RData/mda-create-manifest_%s.RData", timetag)))






