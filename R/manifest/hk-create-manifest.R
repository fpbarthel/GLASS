#######################################################
# Create manifest for hong kong samples (GLASS)
# Date: 2018.09.06
# Author: Kevin J., Floris B.
#######################################################
# Local directory for github repo.
mybasedir = "/Users/johnsk/Documents/Life-History/GLASS-WG"
setwd(mybasedir)

# Files with information about fastq information and barcodes.
hk_file_path = "data/sequencing-information/HK/hong_kong_file_list_complete.20180813.tsv"
life_history_barcodes = "data/ref/glass_wg_aliquots_mapping_table.txt"
hk_clinical_info = "data/clinical-data/HK/recurrent.glioma.HongKong.20161221.xlsx"

# Create extensions for samples.
json_ext = "json"
text_ext = "tsv"

# Create output files for each metadata set.
cases_file      = "data/manifest/hongkong/cases"
samples_file    = "data/manifest/hongkong/samples"
aliquots_file   = "data/manifest/hongkong/aliquots"
readgroups_file = "data/manifest/hongkong/readgroups"
files_file      = "data/manifest/hongkong/files"
pairs_file      = "data/manifest/hongkong/pairs"

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



### Aliquot ####
master_sheet = read.delim(life_history_barcodes, as.is=T)
aliquot_sheet = master_sheet %>% select(aliquot_uuid = uuid, sample_id = Barcode) %>%
  mutate(aliquot_id = sprintf("%s-%s", sample_id, aliquot_uuid),
         analyte_type = "DNA",
         analysis_type = "WGS",
         portion = 1) %>%
  filter(grepl("HK", sample_id))

# Aliquot file to be written.
aliquots = aliquot_sheet %>% select(sample_id, aliquot_uuid, aliquot_id, portion, analyte_type, analysis_type) %>% distinct()



### Files ####
# Generate *file* tsv containing: aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format.
hk_df = read.table(hk_file_path, col.names="filenames", stringsAsFactors = F)

# Retrieve the original sample name to map to study center provided covariate sheet.
hk_df_meta = hk_df %>% 
  mutate(verhaak_sample_id = sub(".*WG_ *(.*?) *_USPD.*", "\\1", filenames),
         library_id = sub(".*_ *(.*?) *_H.*", "\\1", filenames), 
         flowcell_id = substr(filenames, nchar(filenames)-19, nchar(filenames)-12),
         lane_id = substr(filenames, nchar(filenames)-8, nchar(filenames)-8))

# Create a new identifier on which to group mate pairs onto the same line (i.e., R1 and R2).
hk_df_meta$read_group = paste(hk_df_meta$library_id, hk_df_meta$flowcell_id, hk_df_meta$lane_id, sep = '-')

# Create a new identifier on which to group mate pairs onto the same line (i.e., R1 and R2).
hk_df_meta$read_group = paste(hk_df_meta$library_id, hk_df_meta$flowcell_id, hk_df_meta$lane_id, sep = '-')

# Provide path to files temporarily stored on /fastscratch/.
hk_df_meta$fullpath_filenames <- paste("/fastscratch/johnsk/GLASS-WG/hongkong", hk_df_meta$filenames, sep="/")

# Comma separated file_paths and file_names.
merged_hk_files = hk_df_meta %>% 
  group_by(read_group) %>% 
  mutate(file_path = paste(fullpath_filenames, collapse=","),
         file_name = paste(filenames, collapse=",")) %>% 
  select(-filenames, -fullpath_filenames) %>% 
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
# Clinical data to extract subject sex.
hk_clinical_data <- readWorkbook(hk_clinical_info, sheet = 1, startRow = 1, colNames = TRUE)
hk_clinical_data$Sex[1:5]

hk_map_df = hk_map_df %>% 
  mutate(case_id = substring(Barcode, 1, 12), 
         project_id = "GLSS-HK")
# Select only those relevant fields.
cases = hk_map_df %>% select(case_id, project_id) %>% distinct()




### Samples ####
# Grab last two characrters of barcode.
hk_map_df$sample_type = substring(hk_map_df$Barcode, 14, 15)
# Recode variables to match Floris' fields.
samples = hk_map_df %>% select(case_id, sample_id = Barcode, legacy_sample_id = verhaak_sample_id, sample_type) %>% distinct()



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

# Note: Hong Kong samples do not have second recurrences.
pairs = rbind(p1,p2) %>% filter(complete.cases(tumor_aliquot_id, normal_aliquot_id))



### Readgroups ####
# Necessary information: file_uuid, aliquot_id, RGID, RGPL, RGPU, RGLB, RGPI, RGDT, RGSM, RGCN.
readgroup_df = hk_map_df %>% 
  mutate(RGPL = "ILLUMINA",
    RGPU = paste(substr(file_name, nchar(file_name)-19, nchar(file_name)-12), 
                 substr(file_name, nchar(file_name)-8, nchar(file_name)-8), sep="."),
    RGLB = sub(".*_ *(.*?) *_H.*", "\\1", file_name),
    RGPI = 0,
    RGDT = strftime(as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%dT%H:%M:%S%z"), 
    RGSM = Barcode,
    RGCN = "NVGN_HK",
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
save.image(file.path(sprintf("R/RData/hk-create-manifest_%s.RData", timetag)))

