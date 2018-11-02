
library(tidyverse)
library(stringi)
library(DBI)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

md5sums = read.delim("/projects/verhaak-lab/metastatic_gliosarcoma/N06057_MK_MJ1801151/raw_data/md5_postupload.md5", as.is=T, header=F, sep=" ")
rginfo = read.delim("/projects/verhaak-lab/metastatic_gliosarcoma/Sample_Identifiers.txt", as.is=T)
fqfiles = list.files("/projects/verhaak-lab/metastatic_gliosarcoma/N06057_MK_MJ1801151/raw_data", pattern = ".fq.gz$", full.names = T)
fqfiles = fqfiles[-c(2,4)]

aliquots_master = read.delim("data/ref/glass_wg_aliquots_mapping_table.txt", as.is=T)
aliquots_master = aliquots_master %>% mutate(aliquot_barcode = aliquot_id,
                                             sample_barcode = substr(aliquot_barcode, 1, 15),
                                             sample_type = substr(aliquot_barcode, 14, 15),
                                             case_barcode = substr(aliquot_barcode, 1, 12),
                                             aliquot_uuid_short = substr(aliquot_barcode, 25, 30),
                                             aliquot_analyte_type = substr(aliquot_barcode, 19, 19),
                                             aliquot_portion = substr(aliquot_barcode, 17, 18),
                                             aliquot_analysis_type = substr(aliquot_barcode, 21, 23),
                                             aliquot_id_legacy = legacy_aliquot_id) %>%
  filter(grepl("GLSS-NS", sample_barcode))

md5sums = md5sums %>% select(file_name = V3, file_md5sum = V1)

files = data.frame(aliquot_id_legacy = gsub("_L6_\\d(?:.fixed)*.fq.gz","",basename(fqfiles)),
                   file_name = basename(fqfiles),
                   file_size = file.info(fqfiles)$size,
                   file_format = "FASTQ",
                   stringsAsFactors = F)

files = files %>%
  left_join(md5sums) %>%
  left_join(aliquots_master) %>%
  select(aliquot_barcode=aliquot_id, file_name, file_size, file_md5sum, file_format) %>%
  distinct()

readgroups = data.frame(readgroup_id = sprintf("%s.%s", substr(rginfo$Flowcell.ID[seq(1,9,2)],1,5), rginfo$Flowcell.Lane[seq(1,9,2)]),
                        readgroup_platform = "ILLUMINA",
                        readgroup_platform_unit = sprintf("%s.%s", rginfo$Flowcell.ID[seq(1,9,2)], rginfo$Flowcell.Lane[seq(1,9,2)]),
                        readgroup_library = unique(gsub("_L6_\\d.fq.gz","",basename(rginfo$Sample))),
                        legacy_aliquot_id = unique(gsub("_L6_\\d.fq.gz","",basename(rginfo$Sample))),
                        readgroup_center = "GENEWIZ",
                        stringsAsFactors = F)

readgroups = readgroups %>% 
  left_join(select(aliquots_master, aliquot_id, legacy_aliquot_id)) %>%
  select(aliquot_barcode=aliquot_id, readgroup_idtag=readgroup_id, everything(), -legacy_aliquot_id) %>%
  mutate(readgroup_sample_id = aliquot_barcode)


### aliquots
aliquots = aliquots_master %>% 
  select(aliquot_barcode, aliquot_id_legacy, sample_barcode, aliquot_uuid_short,
         aliquot_analyte_type, aliquot_portion, aliquot_analysis_type) %>% 
  distinct()

## Samples
samples = aliquots_master %>% 
  select(case_barcode, sample_barcode, sample_type) %>% 
  distinct()

### Cases
cases = aliquots_master %>% 
  mutate(case_project = "GLSS", age = NA, sex = "female", case_barcode = substr(aliquot_id,1,12), case_source = substr(aliquot_id,6,7)) %>%
  select(case_project, case_barcode, case_source, case_age_diagnosis_years=age, case_sex=sex) %>% distinct()

## Filemap
tmp1 = files %>% select(file_name, aliquot_barcode)
tmp2 = readgroups %>% select(readgroup_idtag, readgroup_sample_id, aliquot_barcode)
files_readgroups = full_join(tmp1, tmp2) %>% select(file_name, readgroup_idtag, readgroup_sample_id)

## Write to database
dbWriteTable(con, Id(schema="clinical",table="cases"), cases, append=T)
dbWriteTable(con, Id(schema="biospecimen",table="samples"), samples, append=T)
dbWriteTable(con, Id(schema="biospecimen",table="aliquots"), aliquots, append=T)
dbWriteTable(con, Id(schema="biospecimen",table="readgroups"), readgroups, append=T)
dbWriteTable(con, Id(schema="analysis",table="files"), files, append=T)
dbWriteTable(con, Id(schema="analysis",table="files_readgroups"), files_readgroups, append=T)

## Close connection
dbDisconnect(con)
