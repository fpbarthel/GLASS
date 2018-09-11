
library(tidyverse)
library(stringi)

cases_file      = "data/manifest/K2/cases"
samples_file    = "data/manifest/K2/samples"
aliquots_file   = "data/manifest/K2/aliquots"
readgroups_file = "data/manifest/K2/readgroups"
files_file      = "data/manifest/K2/files"
pairs_file      = "data/manifest/K2/pairs"

json_ext = "json"
text_ext = "tsv"

rginfo = read.delim("/projects/verhaak-lab/metastatic_gliosarcoma/Sample_Identifiers.txt", as.is=T)
fqfiles = list.files("/projects/verhaak-lab/metastatic_gliosarcoma/N06057_MK_MJ1801151/raw_data", pattern = "[1-2].fq.gz$", full.names = T)

aliquots_master = read.delim("data/ref/glass_wg_aliquots_mapping_table.txt", as.is=T)
aliquots_master = aliquots_master %>% mutate(sample_id = substr(aliquot_id, 1, 15),
                                             sample_type_code = substr(aliquot_id, 14, 15),
                                             case_id = substr(aliquot_id, 1, 12),
                                             aliquot_uuid = substr(aliquot_id, 25, 30),
                                             analyte = substr(aliquot_id, 19, 19),
                                             portion = substr(aliquot_id, 17, 18),
                                             analysis_type = substr(aliquot_id, 21, 23)) %>%
  filter(grepl("GLSS-NS", sample_id))

files = data.frame(legacy_aliquot_id = unique(gsub("_L6_\\d.fq.gz","",basename(fqfiles))),
                   file_path = paste(fqfiles[seq(1,9,2)], fqfiles[seq(2,10,2)], sep=","),
                   file_name = paste(basename(fqfiles)[seq(1,9,2)], basename(fqfiles)[seq(2,10,2)], sep=","),
                   file_size = "NA",
                   file_md5sum = "NA",
                   file_format = "FQ",
                   stringsAsFactors = F)

files = files %>%
  mutate(file_uuid=paste(stri_rand_strings(n(), 8, "[a-z0-9]"),
                         stri_rand_strings(n(), 4, "[a-z0-9]"),
                         stri_rand_strings(n(), 4, "[a-z0-9]"),
                         stri_rand_strings(n(), 4, "[a-z0-9]"),
                         stri_rand_strings(n(), 12, "[a-z0-9]"),
                         sep = "-")) %>% 
  left_join(aliquots_master) %>%
  select(aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format) %>%
  distinct()

readgroups = data.frame(readgroup_id = sprintf("%s.%s", substr(rginfo$Flowcell.ID[seq(1,9,2)],1,5), rginfo$Flowcell.Lane[seq(1,9,2)]),
                        readgroup_platform = "ILLUMINA",
                        readgroup_platform_unit = sprintf("%s.%s", rginfo$Flowcell.ID[seq(1,9,2)], rginfo$Flowcell.Lane[seq(1,9,2)]),
                        readgroup_library = unique(gsub("_L6_\\d.fq.gz","",basename(rginfo$Sample))),
                        readgroup_date = "NA",
                        legacy_aliquot_id = unique(gsub("_L6_\\d.fq.gz","",basename(rginfo$Sample))),
                        readgroup_center = "GENEWIZ",
                        stringsAsFactors = F)

readgroups = readgroups %>% 
  left_join(select(aliquots_master, aliquot_id, legacy_aliquot_id)) %>%
  left_join(select(files,file_uuid,aliquot_id)) %>%
  select(file_uuid, aliquot_id, readgroup_id, everything(), -legacy_aliquot_id) %>%
  mutate(readgroup_sample_id = aliquot_id)


### aliquots
aliquots = aliquots_master %>% 
  select(aliquot_id, legacy_aliquot_id, sample_id, case_id, aliquot_uuid, analyte, portion, analysis_type) %>% 
  distinct()

## Samples
samples = aliquots_master %>% 
  select(case_id, sample_id, sample_type = sample_type_code) %>% 
  distinct()

### Cases
cases = aliquots_master %>% 
  mutate(case_project = "GLSS", age = NA, sex = "female") %>%
  select(case_id, case_project, age, sex) %>% distinct()
# 
# ### aliquots
# aliquots = data.frame(sample_id = sprintf("GLSS-K2-0001-%s", c("TP", "R1", "R2", "M1", "NB")),
#                       aliquot_uuid = uuids,
#                       aliquot_id = sprintf("GLSS-K2-0001-%s-%s", c("TP", "R1", "R2", "M1", "NB"), uuid = uuids),
#                       portion = 1,
#                       analyte_type = "DNA",
#                       analysis_type = "WGS",
#                       stringsAsFactors = F)
# 
# ## Samples
# samples = data.frame(case_id = "GLSS-K2-0001",
#                      sample_id = sprintf("GLSS-K2-0001-%s", c("TP", "R1", "R2", "M1", "NB")),
#                      legacy_sample_id = gsub(".fq.gz", "", basename(fqfiles)[seq(1,9,2)]),
#                      sample_type = c("TP", "R1", "R2", "M1", "NB"),
#                      stringsAsFactors = F)
# 
# ### Cases
# cases = data.frame(case_id = "GLSS-K2-0001",
#                    project_id = "K2",
#                    stringsAsFactors = F)

## Pairs
pairs = data.frame(case_id = "GLSS-NS-0001",
                   pair_id = sprintf("%s-%s-%s", substr(aliquots$aliquot_id[1:4],1,18),substr(aliquots$aliquot_id[5],14,18),"WGS"),
                   tumor_aliquot_id = aliquots$aliquot_id[1:4],
                   normal_aliquot_id = aliquots$aliquot_id[5])

## Make manifest
write(jsonlite::toJSON(files, pretty = T), file = sprintf("%s.%s", files_file, json_ext))
write(jsonlite::toJSON(cases, pretty = T), file = sprintf("%s.%s", cases_file, json_ext))
write(jsonlite::toJSON(samples, pretty = T), file = sprintf("%s.%s", samples_file, json_ext))
write(jsonlite::toJSON(aliquots, pretty = T), file = sprintf("%s.%s", aliquots_file, json_ext))
write(jsonlite::toJSON(readgroups, pretty = T), file = sprintf("%s.%s", readgroups_file, json_ext))
write(jsonlite::toJSON(pairs, pretty = T), file = sprintf("%s.%s", pairs_file, json_ext))

write.table(files, file = sprintf("%s.%s", files_file, text_ext), sep="\t", row.names = F, col.names = T, quote = F)
write.table(cases, file = sprintf("%s.%s", cases_file, text_ext), sep="\t", row.names = F, col.names = T, quote = F)
write.table(samples, file = sprintf("%s.%s", samples_file, text_ext), sep="\t", row.names = F, col.names = T, quote = F)
write.table(aliquots, file = sprintf("%s.%s", aliquots_file, text_ext), sep="\t", row.names = F, col.names = T, quote = F)
write.table(readgroups, file = sprintf("%s.%s", readgroups_file, text_ext), sep="\t", row.names = F, col.names = T, quote = F)
write.table(pairs, file = sprintf("%s.%s", pairs_file, text_ext), sep="\t", row.names = F, col.names = T, quote = F)
