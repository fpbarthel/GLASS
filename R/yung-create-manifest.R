#######################################################
# Create manifest for Yung MD Anderson sample (n = 1, GLASS)
# Date: 2018.08.13
# Author: Kevin J.
#######################################################
# Local directory for github repo.
mybasedir = "/Users/johnsk/Documents/Life-History/GLASS-WG"
setwd(mybasedir)

# Files with information about "Yung" fastq information and barcodes.
# Actual fastq files are stored: /fastscratch/johnsk/GLASS-WG/mdacc/ or on tier2 (long-term).
MDA_batch2_file_path = "data/sequencing-information/MDACC/file_list_C202SC18030593_batch2_n94_20180603.tsv"

# 2018.07.06 Katie Shao (Novogene) provided the following sheet linking libraryID with submitted sample names in the email:
# "Re:RE: Novogene Project Report - Confirm -C202SC18030593-Davis-MDACC-156-libseq-hWGS-WOBI-NVUS2018022505[J1-C202SC18030593-1000]"
novogene_sample_path = "data/sequencing-information/MDACC/Novogene_SIF_14.xlsx"

# Create extensions for samples.
json_ext = "json"
text_ext = "tsv"

# Create output files for each metadata set.
cases_file      = "data/manifest/mdacc/yung/cases"
samples_file    = "data/manifest/mdacc/yung/samples"
aliquots_file   = "data/manifest/mdacc/yung/aliquots"
readgroups_file = "data/manifest/mdacc/yung/readgroups"
files_file      = "data/manifest/mdacc/yung/files"

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

### MDA (Novogene cohort) fastq file list ####
mda_batch2_df <- read.table(MDA_batch2_file_path, col.names="file_path", stringsAsFactors = F)
mda_batch2_df$filename = sapply(strsplit(mda_batch2_df$file_path, "/"), "[[", 3)

# Replace the placeholder for pwd "./" from bash cmd: "find -name ".fq.gz" in parent directory of fastqs. 
mda_batch2_df$file_path <- gsub("^", "/fastscratch/johnsk/GLASS-WG/mdacc/", mda_batch2_df$file_path)

# Create an old sample.id for these subjects to be linked. No longer *necessary*, but kept in case it's helpful.
mda_batch2_df$old_sample_id <- gsub( "_USPD.*$", "", mda_batch2_df$filename)

# Find the "Yung" files. This will be used below in the *files* secyion.
yung_mda_df <- mda_batch2_df[grepl("Yung", mda_batch2_df$filename),]

# Katie Shao (Novogene provided) to line up sample names with files.
novogene_linker = readWorkbook(novogene_sample_path, sheet = 1, startRow = 18, colNames = TRUE)
novogene_linker_unique <- novogene_linker[1:121, ] # 121st sample represents non-GLASS sample.
yung_sample_df <- novogene_linker_unique[grepl("Yung", novogene_linker_unique$`*SampleName`), ]

# Extract the 4-digit SUBJECT identifier.
yung_sample_df$SubjectID = sapply(strsplit(yung_sample_df$`*SampleName`, "-"), "[", 3)

# Map the old to new sample name for the tumors.
yung_tumor_sample_map = yung_sample_df %>% 
  mutate(SeqID = "GLSS") %>% 
  mutate(TSS = "MD") %>% 
  rename(SubjectCode = SubjectID) %>% 
  mutate(TissueCode = "TP") %>% 
  mutate(Original_ID = yung_sample_df$`*SampleName`) %>% 
  unite(Barcode, c(SeqID, TSS, SubjectCode, TissueCode), sep = "-", remove = FALSE) %>% 
  select(Original_ID, Barcode, SeqID, TSS, SubjectCode, TissueCode)

set.seed(42)
yung_tumor_sample_map$uuid = stri_rand_strings(1, 6, "[A-Z0-9]")

### Aliquot ####
aliquot_sheet = yung_tumor_sample_map %>% select(aliquot_uuid = uuid, sample_id = Barcode, Original_ID)  %>%
  mutate(aliquot_id = sprintf("%s-%s", sample_id, aliquot_uuid),
         analyte_type = "DNA",
         analysis_type = "WGS",
         portion = 1)

# Aliquot file to be written.
aliquots = aliquot_sheet %>% select(sample_id, aliquot_uuid, aliquot_id, portion, analyte_type, analysis_type) %>% distinct()

### Files ####
# Generate *file* tsv containing: aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format.

# Retrieve the library, flowcell, and lane id from the filename.
yung_df_meta  = yung_mda_df %>% 
  mutate(library_id = sub(".*g_ *(.*?) *_H.*", "\\1", filename), 
         flowcell_id = substr(filename, nchar(filename)-19, nchar(filename)-12),
         lane_id = substr(filename, nchar(filename)-8, nchar(filename)-8))

# Create a new identifier on which to group mate pairs onto the same line (i.e., R1 and R2).
yung_df_meta$read_group = paste(yung_df_meta$library_id, yung_df_meta$flowcell_id, yung_df_meta$lane_id, sep = '-')

# Comma separated file_paths and file_names. "L#_1" and "L#_2" can be any order. Doesn't yet matter for Snakemake.
merged_yung_files = yung_df_meta %>% 
  group_by(read_group) %>% 
  mutate(file_path = paste(file_path, collapse=","), 
         file_name = paste(filename, collapse=",")) %>% 
  select(-filename) %>% 
  distinct()

# Combine with the new files with the sample name provided by Novogene.
# This action will remove the Yung sample.
yung_map_df = merged_yung_files %>% 
  inner_join(novogene_linker_unique, by = c("library_id" = "Novogene.ID")) %>%  
  ungroup()

# Generate uuids for each of the mdacc files.
# Example from TCGA: 24c6f54a-e7a2-4148-8335-045e3c74096e
set.seed(39)
yung_map_df$file_uuid = paste(stri_rand_strings(dim(yung_map_df)[1], 8, "[a-z0-9]"),
                             stri_rand_strings(dim(yung_map_df)[1], 4, "[a-z0-9]"),
                             stri_rand_strings(dim(yung_map_df)[1], 4, "[a-z0-9]"),
                             stri_rand_strings(dim(yung_map_df)[1], 4, "[a-z0-9]"),
                             stri_rand_strings(dim(yung_map_df)[1], 12, "[a-z0-9]"),
                             sep = "-")

# Sanity check: make sure each is unique. Ought to be 274 at this point.
n_distinct(yung_map_df$file_uuid)

# Need to record the file_size and file_md5sum for these samples.
yung_map_df = yung_map_df %>% mutate(file_size = "NA",
                                   file_md5sum = "NA",
                                   aliquot_id = sprintf("%s-%s", aliquots$sample_id, aliquots$aliquot_uuid), 
                                   file_format = "FQ",
                                   Barcode = substring(yung_map_df$aliquot_id,1,15))

# Order needs to be: # aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format.
files = yung_map_df %>% select(aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format) %>% distinct()


### Cases ####
yung_map_df = yung_map_df %>% 
  mutate(case_id = substring(yung_map_df$aliquot_id, 1, 12), 
         project_id = "GLSS-MD")
# Select only those relevant fields.
cases = yung_map_df %>% select(case_id, project_id) %>% 
  distinct()

### Samples ####
# Grab last two characrters of barcode.
yung_map_df$sample_type = substring(yung_map_df$aliquot_id, 14, 15)
# The Yung sample was not included in the GLASS cohort metadata that we received from Kristin.
# It was only available from the sequencing data release.
yung_map_df$sequencing_submitted_id = substring(yung_map_df$file_name, 1, 27)
# Recode variables to match Floris' fields.
samples = yung_map_df %>% select(case_id, sample_id = Barcode, legacy_sample_id = sequencing_submitted_id, sample_type) %>% distinct()


### Pairs ####


### This sample did not have a tumor-normal pair ###


### Readgroups ####
# Necessary information: file_uuid, aliquot_id, RGID, RGPL, RGPU, RGLB, RGPI, RGDT, RGSM, RGCN.
readgroup_df = yung_map_df %>% 
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
write(jsonlite::toJSON(readgroups, pretty = T), file = sprintf("%s.%s", readgroups_file, json_ext))
write(jsonlite::toJSON(samples, pretty = T), file = sprintf("%s.%s", samples_file, json_ext))

print(sprintf("Exporting manifest as tsv files for visualization ease."))
write.table(aliquots, file = sprintf("%s.%s", aliquots_file, text_ext), sep="\t", row.names = F, col.names = T, quote = F)
write.table(files, file = sprintf("%s.%s", files_file, text_ext), sep="\t", row.names = F, col.names = T, quote = F)
write.table(cases, file = sprintf("%s.%s", cases_file, text_ext), sep="\t", row.names = F, col.names = T, quote = F)
write.table(readgroups, file = sprintf("%s.%s", readgroups_file, text_ext), sep="\t", row.names = F, col.names = T, quote = F)
write.table(samples, file = sprintf("%s.%s", samples_file, text_ext), sep="\t", row.names = F, col.names = T, quote = F)


# Output RData object and timestamp along with package versions.
mysession_info <- devtools::session_info()
timetag = make.names(format(Sys.time(),"t%d_%b_%y_%H%M%S%Z"))
save.image(file.path(sprintf("R/RData/yung-create-manifest_%s.RData", timetag)))


