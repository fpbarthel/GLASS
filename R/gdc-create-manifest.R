## Create a manifest of TCGA whole genome (WGS) files from LGG and GBM cohorts
## Limit to primary-recurrent triplets and 2nd recurrences
## @Author Floris Barthel

library(GenomicDataCommons)
library(listviewer)
library(tidyverse)

setwd("/Volumes/Helix-Projects/GLASS-WG")

cases_file      = "data/manifest/tcga/cases"
samples_file    = "data/manifest/tcga/samples"
aliquots_file   = "data/manifest/tcga/aliquots"
readgroups_file = "data/manifest/tcga/readgroups"
files_file      = "data/manifest/tcga/files"
pairs_file      = "data/manifest/tcga/pairs"

json_ext = "json"
text_ext = "tsv"

## Make sure to include case ids
# "cases.project.project_id" = project (eg. TCGA-LGG)
# "cases.samples.sample_type" = sample type (eg. Primary Tumor)
# "cases.samples.portions.analytes.aliquots.submitter_id" = TCGA barcode

add_fields = c("cases.project.project_id",
               "cases.samples.sample_type", 
               "cases.samples.portions.analytes.aliquots.submitter_id")

## Inspect data types
files(legacy = TRUE) %>% facet(c('data_type')) %>% aggregations()

## Get a list of all WXS/RNA-Seq aligned BAM files from primary tumors
fq = files(legacy = TRUE) %>% 
  GenomicDataCommons::filter( ~ cases.project.project_id %in% c("TCGA-LGG", "TCGA-GBM") & 
                                experimental_strategy == "WGS" & 
                                data_type == "Aligned reads") %>% 
  GenomicDataCommons::select(c(default_fields(files()), add_fields))

message(sprintf("Found %s hits", GenomicDataCommons::count(fq)))

## Extract results
fres = results_all(fq)

## Inspect list using listviewer::jsonedit
jsonedit(fres)

## Flatten nested variables
fres$project = map(fres$cases, "project") %>% map_chr("project_id")
fres$sample_type = map(fres$cases, "samples") %>% map(unlist) %>% map_chr("sample_type")
fres$aliquot_id = map(fres$cases, "samples") %>% map(unlist) %>% map_chr("portions.analytes.aliquots.submitter_id")

## Convert to dataframe (finally!)
df = as.data.frame(fres[-which(names(fres) %in% c("cases", "acl", "analysis"))], stringsAsFactors = F) %>%
  select(id, aliquot_id, project, sample_type, experimental_strategy, file_size, md5sum, file_name, created_datetime, updated_datetime) %>%
  filter(grepl("TCGA", project)) %>%
  mutate(sample_id = substr(aliquot_id,1,16),
         case_id = substr(aliquot_id,1,12),
         format = "BAM")

## Pair primary-recurrent-2ndrecurrence samples
filtered_files = df %>% 
  group_by(sample_id) %>%
  mutate(p = order(file_size, decreasing = T)) %>%
  ungroup() %>%
  group_by(case_id) %>%
  mutate(hasRec = any(sample_type == "Recurrent Tumor")) %>%
  ungroup() %>%
  filter(hasRec, p == 1) %>%
  select(-hasRec, -p, -created_datetime, -updated_datetime, -experimental_strategy) %>%
  mutate(file_size_readable = gdata::humanReadable(file_size, standard="Unix"),
         sample_type_code = recode_factor(substr(aliquot_id, 14,16), "01A" = "TP", "01B" = "TP", "02A" = "R1", "02B" = "R2", "10A" = "NB", "10B" = "NB", "10D" = "NB"))

## Full-join
filtered_files_rgs = filtered_files %>% right_join(rgs)

## Nest
nested_filtered_files_rgs = filtered_files_rgs %>%
  mutate(file_uuid = id, case_project = project, file_md5sum = md5sum, file_format = format) %>%
  select(starts_with("case"), starts_with("sample"), starts_with("file"), starts_with("rg")) %>%
  nest(-starts_with("case"), -starts_with("sample"), -starts_with("file"), .key=readgroups) %>% # starts_with("RG")
  nest(-starts_with("case"), -starts_with("sample"),.key=files) %>%
  nest(-starts_with("case"),.key=samples)

jsonlite::toJSON(nested_filtered_files_rgs, pretty = T)
write(jsonlite::toJSON(nested_filtered_files_rgs, pretty = T), file = unpaired_json)

####################################################################################

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


