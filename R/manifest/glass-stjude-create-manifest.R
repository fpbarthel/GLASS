#######################################################
# Create manifest for GLASS whole exome sequencing samples.
# Date: 2018.10.22
# Author: Kevin J., modified by FP BARTHEL
######    Creation of this manifest reflects the data as of: **10.19.2018**   ###########
#######################################################
# Local directory for github repo.
mybasedir = "/fastscratch/verhaak-lab/GLASS-WG/"
setwd(mybasedir)

#######################################################

# Necessary packages:
library(tidyverse)
#library(DBI)

#######################################################
# Establish connection with Floris' database.
#con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

md5sums = read.delim("data/ref/stjude_md5.txt", as.is = T, header = F, sep=" ") %>%
  select(file_name = V3, file_md5sum = V1) %>% mutate(file_name=gsub("./", "", file_name))

# We need to generate the following fields required by the SNV snakemake pipeline:
### Aliquots ####
stjude_barcodes = "data/ref/stjude_wgs_barcodes.txt"
stjude_ids = read.delim("data/ref/stjude_wgs_barcodes.txt", as.is=T)

# Perform a test run where we are only interested in the Columbia wxs samples.
aliquots_master = stjude_ids %>% mutate(sample_barcode = substr(aliquot_barcode, 1, 15),
                                              sample_type = substr(aliquot_barcode, 14, 15),
                                              case_barcode = substr(aliquot_barcode, 1, 12),
                                              aliquot_uuid_short = substr(aliquot_barcode, 25, 30),
                                              aliquot_analyte_type = substr(aliquot_barcode, 19, 19),
                                              aliquot_portion = substr(aliquot_barcode, 17, 18),
                                              aliquot_analysis_type = substr(aliquot_barcode, 21, 23),
                                              aliquot_id_legacy = legacy_aliquot_id) %>%
  select(-file_path)

### aliquots
aliquots = aliquots_master %>% 
  select(aliquot_barcode, aliquot_id_legacy, sample_barcode, aliquot_uuid_short,
         aliquot_analyte_type, aliquot_analysis_type, aliquot_portion) %>% 
  distinct()

bamfiles = list.files('/fastscratch/johnsk/st_jude_data', pattern = "WholeGenome.bam$", full.names = T)
files = data.frame(aliquot_id_legacy = gsub("(\\w*_[A-Z]{1}).WholeGenome.bam","\\1",basename(bamfiles)),
                   file_name = basename(bamfiles),
                   file_path = bamfiles,
                   file_size = file.info(bamfiles)$size,
                   file_format = "uBAM",
                   stringsAsFactors = F)

files = files %>%
  left_join(md5sums) %>%
  left_join(aliquots_master) %>%
  select(aliquot_barcode, file_name, file_size, file_md5sum, file_format, file_path) %>%
  distinct()

readgroups = data.frame(readgroup_id = "UNKWN.1",
                        readgroup_platform = "ILLUMINA",
                        readgroup_platform_unit = "Unknown",
                        readgroup_library = "Unknown",
                        aliquot_id_legacy = unique(aliquots$aliquot_id_legacy),
                        readgroup_center = "LU",
                        stringsAsFactors = F)

readgroups = readgroups %>% 
  left_join(select(aliquots_master, aliquot_barcode, aliquot_id_legacy)) %>%
  select(aliquot_barcode, readgroup_idtag=readgroup_id, readgroup_platform, readgroup_platform_unit, readgroup_library,
         readgroup_center, -aliquot_id_legacy) %>%
  mutate(readgroup_sample_id = aliquot_barcode)

# Select only those relevant fields. Cases will vary depending on which samples are included.
cases = aliquots_master %>% 
  mutate(case_barcode = substring(aliquot_barcode, 1, 12), 
         case_project = substring(aliquot_barcode, 1, 4),
         case_source = substr(aliquot_barcode, 6, 7)) %>% 
  select(case_project, case_barcode, case_source) %>% 
  distinct()

### Samples ####
# Grab last two characrters of barcode. These samples have multi-sector sampling. So the 
# number of sample_types will be smaller than all of the tumors.
samples = aliquots_master %>% 
  mutate(sample_type = substring(aliquot_barcode, 14, 15),
         case_barcode =  substring(aliquot_barcode, 1, 12),
         sample_barcode = substring(aliquot_barcode, 1, 15)) %>% 
  select(case_barcode, sample_barcode, sample_type) %>% 
  distinct()

## Filemap
tmp1 = files %>% select(file_name, aliquot_barcode)
tmp2 = readgroups %>% select(readgroup_idtag, readgroup_sample_id, aliquot_barcode)
files_readgroups = full_join(tmp1, tmp2) %>% select(file_name, readgroup_idtag, readgroup_sample_id)

## Write to database.
dbWriteTable(con, Id(schema="clinical",table="cases"), cases, append=T)
dbWriteTable(con, Id(schema="biospecimen",table="samples"), samples, append=T)
dbWriteTable(con, Id(schema="biospecimen",table="aliquots"), aliquots, append=T)
dbWriteTable(con, Id(schema="biospecimen",table="readgroups"), readgroups, append=T)
dbWriteTable(con, Id(schema="analysis",table="files"), files, append=T)
dbWriteTable(con, Id(schema="analysis",table="files_readgroups"), files_readgroups, append=T)
