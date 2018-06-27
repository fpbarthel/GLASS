#######################################################
# Create manifest for low pass MD Anderson samples (Roel-JDG)
# Date: 2018.06.26
# Author: Kevin J., Floris B.
#######################################################

# Local directory for github repo.
mybasedir = "/Users/johnsk/Documents/Life-History/GLASS-WG"
setwd(mybasedir)

# Files with information about fastq information and barcodes.
JDG_file_path = "data/sequencing-information/Roel-JDG/Roel_JDG_readgroups.txt"
life_history_barcodes = "data/sequencing-information/master_life_history_uniform_naming_incomplete.txt"

# Create extensions for samples.
json_ext = "json"
text_ext = "tsv"

# Create output files for each metadata set.
cases_file      = "data/manifest/roel-jdg/cases"
samples_file    = "data/manifest/roel-jdg/samples"
aliquots_file   = "data/manifest/roel-jdg/aliquots"
readgroups_file = "data/manifest/roel-jdg/readgroups"
files_file      = "data/manifest/roel-jdg/files"
pairs_file      = "data/manifest/roel-jdg/pairs"

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
master_sheet = read.delim("data/sequencing-information/master_life_history_uniform_naming_incomplete.txt", as.is=T)
# grepl low-pass (LP) samples.
aliquot_sheet = master_sheet %>% select(aliquot_uuid = uuid, sample_id = Barcode) %>%
  mutate(aliquot_id = sprintf("%s-%s", sample_id, aliquot_uuid),
         analyte_type = "DNA",
         analysis_type = "WGS",
         portion = 1) %>%
  filter(grepl("LP", sample_id))

# Aliquot file to be written.
aliquots = aliquot_sheet %>% select(sample_id, aliquot_uuid, aliquot_id, portion, analyte_type, analysis_type) %>% distinct()

### Files ####
# Generate *file* tsv containing: aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format.
# hk_df = read.table(HK_file_path, header=T, stringsAsFactors = F)
jdg_df = read.table(JDG_file_path, stringsAsFactors = F)
colnames(jdg_df) = c("file_path", "@RG" , "RGID", "RGPL", "RGLB", "RGSM", "RGCN")
# Remove RG identifier before colon.
drop_prefix = function(x) {gsub(".*:","", x) }
jdg_df[3:7] <- lapply(jdg_df[3:7], drop_prefix)

# Retrieve the file_name from the file_path. 
jdf_map_df = jdg_df %>%  
  mutate(file_name = sapply(strsplit(jdg_df$file_path, "/"), "[[", 6),
         legacy_sample_id = sub(".*ROEL-JDG- *(.*?) *_.*", "\\1", file_name)) %>% 
  inner_join(master_sheet, by = c("legacy_sample_id" = "Original_ID")) 

# Generate uuids for each of the roel-jdg  files.
# Example from TCGA: 24c6f54a-e7a2-4148-8335-045e3c74096e
set.seed(1)
jdf_map_df$file_uuid = paste(stri_rand_strings(dim(jdf_map_df)[1], 8, "[a-z0-9]"),
                            stri_rand_strings(dim(jdf_map_df)[1], 4, "[a-z0-9]"),
                            stri_rand_strings(dim(jdf_map_df)[1], 4, "[a-z0-9]"),
                            stri_rand_strings(dim(jdf_map_df)[1], 4, "[a-z0-9]"),
                            stri_rand_strings(dim(jdf_map_df)[1], 12, "[a-z0-9]"),
                            sep = "-")

# Sanity check: make sure each is unique.
n_distinct(jdf_map_df$file_uuid)

# Need to record the file_size and file_md5sum for these samples. Floris says not necessary.
jdf_map_df = jdf_map_df %>% mutate(file_size = "NA",
                                 file_md5sum = "NA",
                                 aliquot_id = sprintf("%s-%s", Barcode, uuid), 
                                 file_format = "BAM")

# Order needs to be: # aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format.
files = jdf_map_df %>% select(aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format) %>% distinct()


### Cases ####
jdf_map_df = jdf_map_df %>% 
  mutate(case_id = substring(Barcode, 1, 12), 
         project_id = "GLSS-MD")
# Select only those relevant fields.
cases = jdf_map_df %>% select(case_id, project_id)

### Samples ####
# Grab last two characrters of barcode.
jdf_map_df$sample_type = substring(jdf_map_df$Barcode, 14, 15)
# Recode variables to match Floris' fields.
samples = jdf_map_df %>% select(case_id, sample_id = Barcode, legacy_sample_id, sample_type) %>% distinct()


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
readgroup_df = jdf_map_df %>% 
  mutate(RGPU = paste(sub(".*[0-9]{4}_ *(.*?) *_s.*", "\\1", jdf_map_df$file_name), 
                      sub(".*s_ *(.*?) *_[A-Za-z]{2}.*", "\\1", jdf_map_df$file_name), sep="."),
         RGPI = 0,
         RGDT = strftime(as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%dT%H:%M:%S%z"))

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
save.image(file.path(sprintf("R/RData/jdg-create-manifest_%s.RData", timetag)))









