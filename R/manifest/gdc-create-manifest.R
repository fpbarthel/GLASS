## Create a manifest of TCGA whole genome (WGS) files from LGG and GBM cohorts
## Limit to primary-recurrent triplets and 2nd recurrences
## @Author Floris Barthel

library(tidyverse)

## Flowcell ID parsing taken from
## https://github.com/10XGenomics/supernova/blob/master/tenkit/lib/python/tenkit/illumina_instrument.py
# platform_unit_to_model <- function()

# FCIDs = c("C[A-Z,0-9]{4}ANXX$" = "HiSeq 1500/2000/2500",
#           "C[A-Z,0-9]{4}ACXX$" = "HiSeq 1000/1500/2000/2500",
#           "H[A-Z,0-9]{4}ADXX$" = "HiSeq 1500/2500",
#           "H[A-Z,0-9]{4}BCXX$" = "HiSeq 1500/2500",
#           "H[A-Z,0-9]{4}BCXY$" = "HiSeq 1500/2500",
#           "H[A-Z,0-9]{4}BBXX$" = "HiSeq 4000",
#           "H[A-Z,0-9]{4}BBXY$" = "HiSeq 4000",
#           "H[A-Z,0-9]{4}CCXX$" = "HiSeq X",
#           "H[A-Z,0-9]{4}CCXY$" = "HiSeq X",
#           "H[A-Z,0-9]{4}ALXX$" = "HiSeq X",
#           "H[A-Z,0-9]{4}BGXX$" = "NextSeq",
#           "H[A-Z,0-9]{4}BGXY$" = "NextSeq",
#           "H[A-Z,0-9]{4}BGX2$" = "NextSeq",
#           "H[A-Z,0-9]{4}AFXX$" = "NextSeq",
#           "A[A-Z,0-9]{8}$" = "MiSeq",
#           "B[A-Z,0-9]{8}$" = "MiSeq",
#           "D[A-Z,0-9]{8}$" = "MiSeq",
#           "G[A-Z,0-9]{8}$" = "MiSeq",
#           "H[A-Z,0-9]{4}DMXX$" = "NovaSeq")

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
aliquots_master = read.delim("data/ref/glass_wg_aliquots_mapping_table.txt", as.is=T)
aliquots_master = aliquots_master %>% mutate(sample_id = substr(aliquot_id, 1, 15),
                                             sample_type_code = substr(aliquot_id, 14, 15),
                                             case_id = substr(aliquot_id, 1, 12),
                                             aliquot_uuid = substr(aliquot_id, 25, 30),
                                             analyte = substr(aliquot_id, 19, 19),
                                             portion = substr(aliquot_id, 17, 18),
                                             analysis_type = substr(aliquot_id, 21, 23)) %>%
  filter(grepl("TCGA", sample_id))

## Files
files_master = jsonlite::read_json("data/ref/TCGA_WGS_GDC_legacy_UUIDs.json", simplifyVector=T)
files_master = files_master %>% #unnest(samples) %>% unnest(files) %>% unnest(readgroups) %>% mutate(sample_id = sprintf("%s-%s", case_id, sample_type_code)) %>%
  select(case_project=project,file_uuid=id,file_md5sum=md5sum,file_name,file_size,legacy_aliquot_id=aliquot_id,age,sex) %>%
  mutate(file_format="BAM") %>%
  left_join(aliquots_master)

files = files_master %>% 
  mutate(file_path = sprintf("/fastscratch/verhaak-lab/GLASS-WG/data/bam/%s/%s", file_uuid, file_name)) %>%
  select(aliquot_id, file_path, file_name, file_uuid, file_size, file_md5sum, file_format) %>%
  distinct()

## Readgroups
rgs = readLines('data/ref/TCGA_BAM_readgroups.txt')
readgroups = data.frame(file_uuid = basename(dirname(gsub("\\t@RG.*$","",rgs))),
                        legacy_readgroup_id = gsub("^.*ID:([\\w\\.\\-\\:]+).*$", "\\1", rgs, perl=T),
                        readgroup_platform = toupper(gsub("^.*PL:([\\w\\.\\-\\:]+).*$", "\\1", rgs, perl=T)),
                        readgroup_platform_unit = substr(gsub("^.*PU:([\\w\\.\\-\\:]+).*$", "\\1", rgs, perl=T),1,17),
                        readgroup_library = gsub("^.*LB:([\\w\\.\\-\\:]+).*$", "\\1", rgs, perl=T),
                        readgroup_date = gsub("^.*DT:([\\w\\.\\-\\:]+).*$", "\\1", rgs, perl=T),
                        readgroup_center = gsub("^.*CN:([\\w\\.\\-\\:]+).*$", "\\1", rgs, perl=T),
                        stringsAsFactors = F) %>%
  left_join(select(files, file_uuid, aliquot_id)) %>%
  mutate(readgroup_id = sprintf("%s.%s", substr(readgroup_platform_unit,1,5), substr(readgroup_platform_unit,17,17)),
         readgroup_sample_id = aliquot_id,
         readgroup_center = ifelse(readgroup_center == "broad", "BI", readgroup_center)) %>%
  select(file_uuid, aliquot_id, readgroup_id, everything())

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
  left_join(select(files_master,legacy_aliquot_id,age,sex)) %>%
  mutate(case_project = "TCGA") %>%
  select(case_id, case_project, age, sex) %>% distinct()

## Pairs
p1 = samples %>% 
  left_join(aliquots) %>%
  select(sample_type, aliquot_id, case_id) %>%
  filter(sample_type %in% c("TP", "NB")) %>% 
  spread(sample_type, aliquot_id) %>%
  mutate(pair_id = sprintf("%s-%s-%s-%s", case_id, substr(TP, 14, 18), substr(NB, 14, 18), substr(TP, 21, 23))) %>%
  select(case_id, pair_id, tumor_aliquot_id = TP, normal_aliquot_id = NB)

p2 = samples %>% 
  left_join(aliquots) %>%
  select(sample_type, aliquot_id, case_id) %>%
  filter(sample_type %in% c("R1", "NB")) %>%
  spread(sample_type, aliquot_id) %>%
  mutate(pair_id = sprintf("%s-%s-%s-%s", case_id, substr(R1, 14, 18), substr(NB, 14, 18), substr(R1, 21, 23))) %>%
  select(case_id, pair_id, tumor_aliquot_id = R1, normal_aliquot_id = NB)

p3 = samples %>% 
  left_join(aliquots) %>%
  select(sample_type, aliquot_id, case_id) %>%
  filter(sample_type %in% c("R2", "NB")) %>%
  spread(sample_type, aliquot_id) %>%
  mutate(pair_id = sprintf("%s-%s-%s-%s", case_id, substr(R2, 14, 18), substr(NB, 14, 18), substr(R2, 21, 23))) %>%
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
# renamemap = samples %>% left_join(aliquots) %>%
#   mutate(old_bam = sprintf("%s.realn.dedup.bqsr.bam",legacy_sample_id),
#          old_bai = sprintf("%s.realn.dedup.bqsr.bai",legacy_sample_id),
#          old_md5 = sprintf("%s.realn.dedup.bqsr.bam.md5",legacy_sample_id),
#          old_bqsr = sprintf("%s.bqsr.txt",legacy_sample_id),
#          new_bam = sprintf("%s.realn.mdup.bqsr.bam",aliquot_id),
#          new_bai = sprintf("%s.realn.mdup.bqsr.bai",aliquot_id),
#          new_md5 = sprintf("%s.realn.mdup.bqsr.bam.md5",aliquot_id),
#          new_bqsr = sprintf("%s.bqsr.txt",aliquot_id)) %>%
#   select(aliquot_id, starts_with("new"), starts_with("old")) %>%
#   gather(starts_with("new"), starts_with("old"), key="var", value="filename") %>%
#   separate(var, into=c("a","file"), sep="_") %>%
#   spread(a, filename) %>%
#   select(old,new)
# 
# write.table(renamemap, file="results_bqsr_renamemap.tsv", quote=F, row.names=F, sep="\t")

mysession_info <- devtools::session_info()

timetag = make.names(format(Sys.time(),"t%d_%b_%y_%H%M%S%Z"))
save.image(file.path(sprintf("R/RData/gdc-create-manifest_%s.RData", timetag)))
print(sprintf("Done! Manifest RData file saved to R/RData/gdc-create-manifest_%s.RData\nManifest files in json formats are at data/manifest/tcga/", timetag))

## end ##
