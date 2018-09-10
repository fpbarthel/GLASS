#######################################################
# Create manifest for low pass Henry Ford sample (hGBM)
# Date: 2018.09.10
# Author: Kevin J., Floris B.
#######################################################
# Local directory for github repo.
mybasedir = "/Users/johnsk/Documents/Life-History/GLASS-WG"
setwd(mybasedir)

# Files with information about fastq information and barcodes.
hf_file_path = "data/sequencing-information/HenryFord/henry_ford_readgroups.txt"
life_history_barcodes = "data/ref/glass_wg_aliquots_mapping_table.txt"
hf_clinical_info = "data/clinical-data/HF/HGBM-clinical-HF3016.xlsx"

# Create extensions for samples.
json_ext = "json"
text_ext = "tsv"

# Create output files for each metadata set.
cases_file      = "data/manifest/henryford/cases"
samples_file    = "data/manifest/henryford/samples"
aliquots_file   = "data/manifest/henryford/aliquots"
readgroups_file = "data/manifest/henryford/readgroups"
files_file      = "data/manifest/henryford/files"
pairs_file      = "data/manifest/henryford/pairs"

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

### Aliquots ####
life_history_ids = read.delim(life_history_barcodes, as.is=T)
aliquots_master = life_history_ids %>% mutate(sample_id = substr(aliquot_id, 1, 15),
                                              sample_type_code = substr(aliquot_id, 14, 15),
                                              case_id = substr(aliquot_id, 1, 12),
                                              aliquot_uuid = substr(aliquot_id, 25, 30),
                                              analyte = substr(aliquot_id, 19, 19),
                                              portion = substr(aliquot_id, 17, 18),
                                              analysis_type = substr(aliquot_id, 21, 23)) %>%
  filter(grepl("GLSS-HF", sample_id))

### Aliquots ###
aliquots = aliquots_master %>% 
  select(aliquot_id, legacy_aliquot_id, sample_id, case_id, aliquot_uuid, analyte, portion, analysis_type) %>% 
  distinct()

### Files ####
# Generate *file* tsv containing: aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format.
hf_df = read.table(hf_file_path, stringsAsFactors = F)
colnames(hf_df) = c("file_path", "@RG" , "RGID", "RGPL", "RGLB", "RGSM", "RGCN")
# Remove RG identifier before colon.
drop_prefix = function(x) {gsub(".*:","", x) }
hf_df[3:7] <- lapply(hf_df[3:7], drop_prefix)

# Retrieve the file_name from the file_path. 
hf_map_df = hf_df %>%  
  mutate(file_name = sapply(strsplit(hf_df$file_path, "/"), "[[", 6),
  legacy_sample_id = sub(".*FORD- *(.*?) *_1.*", "\\1", file_name)) %>% 
  inner_join(aliquots_master, by = c("legacy_sample_id" = "legacy_aliquot_id")) 

# Generate uuids for each of the henry ford files.
# Example from TCGA: 24c6f54a-e7a2-4148-8335-045e3c74096e
hf_map_df$file_uuid = paste(stri_rand_strings(dim(hf_map_df)[1], 8, "[a-z0-9]"),
                            stri_rand_strings(dim(hf_map_df)[1], 4, "[a-z0-9]"),
                            stri_rand_strings(dim(hf_map_df)[1], 4, "[a-z0-9]"),
                            stri_rand_strings(dim(hf_map_df)[1], 4, "[a-z0-9]"),
                            stri_rand_strings(dim(hf_map_df)[1], 12, "[a-z0-9]"),
                            sep = "-")

# Sanity check: make sure each is unique.
n_distinct(hf_map_df$file_uuid)

# Need to record the file_size and file_md5sum for these samples. Floris says not necessary.
hf_map_df = hf_map_df %>% mutate(file_size = "NA",
                                 file_md5sum = "NA",
                                 file_format = "BAM")

# Order needs to be: # aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format.
files = hf_map_df %>% select(aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format) %>% distinct()


### Cases ####
# Clinical data to extract subject sex.
hf_clinical_data <- readWorkbook(hf_clinical_info, sheet = 1, startRow = 1, colNames = TRUE)

# Prepare data to be merged.
hf_clinical = hf_clinical_data %>% 
  mutate(age = `Age.@.dx`,
         sex = recode(Gender, "M"="male"),
         GLSS_case_id = "GLSS-HF-3016") %>% 
  distinct(GLSS_case_id, .keep_all = TRUE) %>% 
  select(GLSS_case_id, age, sex)

# Merge to lift over `sex` variable.
hf_map_df = hf_map_df %>% 
  mutate(case_id = substring(aliquot_id, 1, 12), 
         case_project = "GLSS-HF") %>% 
  inner_join(hf_clinical, by=c("case_id" = "GLSS_case_id"))

# Select only those relevant fields.
cases = hf_map_df %>% 
  select(case_id, case_project, age, sex) %>% 
  distinct()



### Samples ####
# Grab last two characrters of barcode.
hf_map_df$sample_type = substring(hf_map_df$aliquot_id, 14, 15)

# Recode variables to match Floris' fields.
samples = hf_map_df %>% select(case_id, sample_id, sample_type) %>% distinct()


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

# Note: The Henry Ford sample HF-3016 originally had two NB, we selected HF-3016-10-28D because
# it was presumed to be taken at first timepoint.
pairs = rbind(p1,p2) %>% filter(complete.cases(tumor_aliquot_id, normal_aliquot_id))



### Readgroups ####
# Necessary information: file_uuid, aliquot_id, RGID, RGPL, RGPU, RGLB, RGPI, RGDT, RGSM, RGCN.
# *** Note: had to generate a new time stamp because it was not found in the RG header for the bams. ****
# Revised: Can't rename RGID, but pipeline is looking for old RGID.
readgroup_df = hf_map_df %>% 
  mutate(readgroup_platform = "ILLUMINA",
         legacy_readgroup_id = RGID,
         readgroup_platform_unit = paste(sub(".*[0-9]{4}_ *(.*?) *_s.*", "\\1", hf_map_df$file_name),
                                         sub(".*s_ *(.*?) *_[A-Za-z]{2}.*", "\\1", hf_map_df$file_name), sep="."),
         readgroup_date = strftime(as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%dT%H:%M:%S%z"),
         readgroup_library = RGLB,
         readgroup_center = RGCN,
         readgroup_sample_id = RGSM, 
         readgroup_id = paste0(substring(readgroup_platform_unit, 1, 5), substring(readgroup_platform_unit, nchar(readgroup_platform_unit)-1, nchar(readgroup_platform_unit)), ""))
         

# Finalize readgroup information in predefined order.
readgroups = readgroup_df %>% select(file_uuid, aliquot_id, readgroup_id, legacy_readgroup_id, readgroup_platform, readgroup_platform_unit, readgroup_library, readgroup_date, readgroup_center, readgroup_sample_id) %>% distinct()


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
save.image(file.path(sprintf("R/RData/hf-create-manifest_%s.RData", timetag)))

