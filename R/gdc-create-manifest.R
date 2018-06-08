## Create a manifest of TCGA whole genome (WGS) files from LGG and GBM cohorts
## Limit to primary-recurrent triplets and 2nd recurrences
## @Author Floris Barthel

library(GenomicDataCommons)
library(listviewer)
library(tidyverse)

## record base directory, typically ./.git dir is present
mybasedir = here::here()
setwd(mybasedir)
getwd()

cases_file      = "data/manifest/tcga/cases"
samples_file    = "data/manifest/tcga/samples"
aliquots_file   = "data/manifest/tcga/aliquots"
readgroups_file = "data/manifest/tcga/readgroups"
files_file      = "data/manifest/tcga/files"
pairs_file      = "data/manifest/tcga/pairs"

json_ext = "json"
text_ext = "tsv"

## Aliquots
df2 = read.delim("data/sequencing-information/master_life_history_uniform_naming_incomplete.txt", as.is=T)
df2 = df2 %>% select(aliquot_uuid = uuid, sample_id = Barcode) %>%
  mutate(aliquot_id = sprintf("%s-%s", sample_id, aliquot_uuid),
         analyte_type = "DNA",
         analysis_type = "WGS",
         portion = 1) %>%
  filter(grepl("TCGA", sample_id))

tmp = jsonlite::read_json("data/ref/TCGA_WGS_GDC_legacy_UUIDs.json", simplifyVector=T)
df = tmp %>% unnest(samples) %>%
  unnest(files) %>%
  unnest(readgroups) %>%
  mutate(legacy_sample_id = sample_id,
         sample_id = sprintf("%s-%s", case_id, sample_type_code)) %>%
  left_join(df2)

files = df %>% 
  mutate(file_path = sprintf("/fastscratch/barthf/GLASS-WG/download/%s/%s", file_uuid, file_name)) %>%
  select(aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format) %>%
  distinct()

## Readgroups
rgs = readLines('data/ref/TCGA_BAM_readgroups.txt')
readgroups = data.frame(file_uuid = basename(dirname(gsub("\\t@RG.*$","",rgs))),
                        RGID = gsub("^.*ID:([\\w\\.\\-\\:]+).*$", "\\1", rgs, perl=T),
                        RGPL = gsub("^.*PL:([\\w\\.\\-\\:]+).*$", "\\1", rgs, perl=T),
                        RGPU = gsub("^.*PU:([\\w\\.\\-\\:]+).*$", "\\1", rgs, perl=T),
                        RGLB = gsub("^.*LB:([\\w\\.\\-\\:]+).*$", "\\1", rgs, perl=T),
                        RGPI = gsub("^.*PI:([\\w\\.\\-\\:]+).*$", "\\1", rgs, perl=T),
                        RGDT = gsub("^.*DT:([\\w\\.\\-\\:]+).*$", "\\1", rgs, perl=T),
                        RGSM = gsub("^.*SM:([\\w\\.\\-\\:]+).*$", "\\1", rgs, perl=T),
                        RGCN = gsub("^.*CN:([\\w\\.\\-\\:]+).*$", "\\1", rgs, perl=T),
                        stringsAsFactors = F) %>%
  mutate(RGPI = ifelse(RGPI == 0, 0, NA)) %>%
  left_join(select(files, file_uuid, aliquot_id)) %>%
  select(file_uuid, aliquot_id, everything())

### aliquots
aliquots = df %>% select(sample_id, aliquot_uuid, aliquot_id, portion, analyte_type, analysis_type) %>% distinct()

## Samples
samples = df %>% select(case_id, sample_id, legacy_sample_id, sample_type = sample_type_code) %>% distinct()

### Cases
cases = df %>% select(case_id, project_id = case_project)

## Pairs
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

pairs = rbind(p1,p2,p3) %>% filter(complete.cases(tumor_aliquot_id, normal_aliquot_id))

print(sprintf("Exporting manifest as json files"))
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

## Old - new name mapping table
## hack that works bc for now samples and aliquots map 1:1
renamemap = samples %>% left_join(aliquots) %>%
  mutate(old_bam = sprintf("%s.realn.dedup.bqsr.bam",legacy_sample_id),
         old_bai = sprintf("%s.realn.dedup.bqsr.bai",legacy_sample_id),
         old_md5 = sprintf("%s.realn.dedup.bqsr.bam.md5",legacy_sample_id),
         old_bqsr = sprintf("%s.bqsr.txt",legacy_sample_id),
         new_bam = sprintf("%s.realn.mdup.bqsr.bam",aliquot_id),
         new_bai = sprintf("%s.realn.mdup.bqsr.bai",aliquot_id),
         new_md5 = sprintf("%s.realn.mdup.bqsr.bam.md5",aliquot_id),
         new_bqsr = sprintf("%s.bqsr.txt",aliquot_id)) %>%
  select(aliquot_id, starts_with("new"), starts_with("old")) %>%
  gather(starts_with("new"), starts_with("old"), key="var", value="filename") %>%
  separate(var, into=c("a","file"), sep="_") %>%
  spread(a, filename) %>%
  select(old,new)

write.table(renamemap, file="results_bqsr_renamemap.tsv", quote=F, row.names=F, sep="\t")

mysession_info <- devtools::session_info()

timetag = make.names(format(Sys.time(),"t%d_%b_%y_%H%M%S%Z"))
save.image(file.path(sprintf("R/RData/gdc-create-manifest_%s.RData", timetag)))
print(sprintf("Done! Manifest RData file saved to R/RData/gdc-create-manifest_%s.RData\nManifest files in json formats are at data/manifest/tcga/", timetag))

## end ##
