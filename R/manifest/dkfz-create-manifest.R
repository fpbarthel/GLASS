##############################################
# Generate new barcodes and readgroups for DKFZ GLASS data (n=15).
# Updated: 2019.02.20
# Author: Kevin J.
##################################################

# Local directory for github repo.
mybasedir = "/Users/johnsk/Documents/Life-History/GLASS-WG/"
setwd(mybasedir)

#######################################################
# Necessary packages:
library(tidyverse)
library(DBI)
library(openxlsx)

######################################################## 
#Establish connection with the GLASS database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

# Load essential tables. This will be helpful so that you know which fields are required for each table.
cases = dbReadTable(con,  Id(schema="clinical", table="cases"))
case_sources = dbReadTable(con,  Id(schema="clinical", table="case_sources"))
aliquots = dbReadTable(con,  Id(schema="biospecimen", table="aliquots"))
samples = dbReadTable(con,  Id(schema="biospecimen", table="samples"))
readgroups = dbReadTable(con,  Id(schema="biospecimen", table="readgroups"))
files = dbReadTable(con,  Id(schema="analysis", table="files"))
files_readgroups = dbReadTable(con,  Id(schema="analysis", table="files_readgroups"))

# Read in the barcodes that Anzhela prepared:
batch2_barcodes = read.xlsx("/Users/johnsk/Documents/Life-History/glass-batch2/Barcodes_20190219.xlsx")
batch2_barcodes = data.frame(lapply(batch2_barcodes, trimws), stringsAsFactors = FALSE)
# We are missing surgery information for these -DK- samples at this time.
batch2_surgery = read.xlsx("/Users/johnsk/Documents/Life-History/glass-batch2/clinic_surgery20190219.xlsx")
batch2_surgery = data.frame(lapply(batch2_surgery, trimws), stringsAsFactors = FALSE)
# There should be 15 new cases.
batch2_cases = read.xlsx("/Users/johnsk/Documents/Life-History/glass-batch2/clinic_caseAM20190219.xlsx")
batch2_cases = data.frame(lapply(batch2_cases, trimws), stringsAsFactors = FALSE)

# Restrict the barcodes to just the DKFZ samples.
dkfz_barcodes = batch2_barcodes %>% 
  filter(grepl("-DK-", sample_barcode)) %>% 
  select(patient_barcode:uuid)

# DKFZ calculated md5 checksums. NOTE: This includes the updated checksums after the globus transfer of one file (see "dkfz_checksums.R").
dkfz_md5_jax = read.delim("/Volumes/johnskhelix/DKFZ_20190214/GLASS_transfer/file_checksums/jax_file_checksums.txt", sep = "", header = FALSE, stringsAsFactors = FALSE)
dkfz_md5_jax_fqs = dkfz_md5_jax %>% 
  filter(grepl(".fastq.gz", V2))

# Retrieve the information about the files, including path and size.
fqfiles = list.files('/Volumes/johnskhelix/DKFZ_20190214/GLASS_transfer/fastqs/', pattern = ".fastq.gz$", full.names = T)
files = data.frame(aliquot_id_legacy = gsub("_complete_filtered.fastq.gz", "", basename(fqfiles)),
                   file_name = basename(fqfiles),
                   file_path = gsub("/Volumes/johnskhelix/", "/projects/johnsk/", fqfiles),
                   file_size = file.info(fqfiles)$size,
                   file_format = "FASTQ",
                   stringsAsFactors = F)

#### Files ######
dkfz_files = files %>% 
  # For some reason an extra blackslash was added. Remove it.
  inner_join(dkfz_md5_jax_fqs, by=c("file_name"="V2")) %>% 
  mutate(aliquot_barcode = substr(file_name, 1, 30),
         file_path = gsub("fastqs//", "fastqs/", file_path)) %>% 
  select(aliquot_barcode, file_path, file_name, file_size, file_md5sum = V1, file_format) %>%
  distinct()

#### Case Sources ####
dkfz_case_sources <-data.frame("DK","Deutsches Krebsforschungszentrum (German Cancer Research Center, Heidelberg, Germany) ")
colnames(dkfz_case_sources) <- colnames(case_sources)

### Aliquots ####
dkfz_aliquots = dkfz_barcodes %>% 
  mutate(aliquot_barcode = sample_barcode,
         # Not including legacy IDs because it was never explicitly shared with us.
         aliquot_id_legacy = NA,
         sample_barcode = substr(aliquot_barcode, 1, 15), 
         aliquot_uuid_short = substr(aliquot_barcode, 25, 30), 
         aliquot_analyte_type = substr(aliquot_barcode, 19, 19), 
         aliquot_analysis_type = substr(aliquot_barcode, 21, 23),
         aliquot_portion = as.integer(substr(aliquot_barcode, 17, 18)),
         aliquot_batch = paste(substr(aliquot_barcode, 1, 7), substr(aliquot_barcode, 21, 23), sep="-")) %>% 
  select(aliquot_barcode, aliquot_id_legacy, sample_barcode, aliquot_uuid_short, 
         aliquot_analyte_type, aliquot_analysis_type, aliquot_portion, aliquot_batch)

### Cases ####
dkfz_cases = batch2_cases %>% 
  filter(grepl("DK", case_source)) 

### Samples ####
# Grab last two characrters of barcode. 
dkfz_samples = dkfz_barcodes %>% 
  mutate(sample_type = substring(sample_barcode, 14, 15),
         case_barcode =  substring(sample_barcode, 1, 12),
         sample_barcode = substring(sample_barcode, 1, 15)) %>% 
  select(case_barcode, sample_barcode, sample_type) %>% 
  distinct()

##### Surgery #######
dkfz_surgeries = batch2_surgery %>% 
  filter(grepl("-DK-", case_barcode)) 

# Determine the readgroup information from the fastq file (use header first line information):
dkfz_fastq_info = read.delim("/Volumes/johnskhelix/DKFZ_20190214/dfkz_fastq_information.txt", sep = "", header = FALSE, stringsAsFactors = FALSE)
colnames(dkfz_fastq_info) <- c("file_name", "seq_info", "mate")
dkfz_fastq_info = dkfz_fastq_info %>% 
  mutate(aliquot_barcode = substr(file_name, 1, 30), 
         flowcell_id = sapply(strsplit(dkfz_fastq_info$seq_info, ":"), "[[", 3),
         lane_id = sapply(strsplit(dkfz_fastq_info$seq_info, ":"), "[[", 4))

# Restrict to a specific set of variables.
dkfz_readgroups_all = dkfz_fastq_info %>% 
  mutate(readgroup_platform = "ILLUMINA",
         readgroup_platform_unit = paste(flowcell_id, lane_id, sep="."),
         readgroup_library = "Unknown",
         readgroup_timestamp = strftime(as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%dT%H:%M:%S%z"), 
         readgroup_sample_id = aliquot_barcode,
         readgroup_center = "DK",
         readgroup_idtag = paste(substr(readgroup_platform_unit, 1, 5), lane_id, sep="."),
         readgroup_idtag_legacy = NA, 
         readgroup_sample_id = aliquot_barcode)  
# Restrict to the readgroup columns that we need.
dkfz_readgroups = dkfz_readgroups_all %>%   
select(aliquot_barcode, readgroup_idtag, readgroup_idtag_legacy, readgroup_platform, readgroup_platform_unit, readgroup_library, 
         readgroup_center, readgroup_sample_id, readgroup_timestamp) %>% 
  distinct()

#### files_readgroups #####
dkfz_files_readgroups = dkfz_readgroups_all %>% select(file_name, readgroup_idtag, readgroup_sample_id) %>% distinct()

### OUTPUT ####
# Only upload NEW cases. You can query the database and perform an antijoin.
extant_cases <- dbReadTable(con,  Id(schema="clinical",table="cases"))
new_cases = dkfz_cases %>% 
  anti_join(extant_cases, by="case_barcode")
# New surgeries
extant_surgeries <- dbReadTable(con,  Id(schema="clinical",table="surgeries"))
new_surgeries = dkfz_surgeries %>% 
  anti_join(extant_surgeries, by="case_barcode")
# New samples.
extant_samples <- dbReadTable(con,  Id(schema="biospecimen",table="samples"))
new_samples = dkfz_samples %>% 
  anti_join(extant_samples, by="sample_barcode")
# New samples.
extant_aliquots <- dbReadTable(con,  Id(schema="biospecimen",table="aliquots"))
new_aliquots = dkfz_aliquots %>% 
  anti_join(extant_aliquots, by="aliquot_barcode")

## Write to database.
dbWriteTable(con, Id(schema="clinical", table="case_sources"), dkfz_case_sources, append=T)
dbWriteTable(con, Id(schema="clinical", table="cases"), new_cases, append=T)
dbWriteTable(con, Id(schema="biospecimen", table="samples"), new_samples, append=T)
dbWriteTable(con, Id(schema="clinical", table="surgeries"), new_surgeries, append=T)
dbWriteTable(con, Id(schema="biospecimen", table="aliquots"), new_aliquots, append=T)
dbWriteTable(con, Id(schema="biospecimen", table="readgroups"), dkfz_readgroups, append=T)
dbWriteTable(con, Id(schema="analysis", table="files"), dkfz_files, append=T)
dbWriteTable(con, Id(schema="analysis", table="files_readgroups"), dkfz_files_readgroups, append=T)

# From Anzhela's table there were three aliquots with incorrect UUIDs that needed to be updated.
rs = dbSendStatement(con, "UPDATE biospecimen.aliquots SET aliquot_barcode = 'GLSS-DK-0001-R1-01D-WXS-D6A176' , aliquot_uuid_short = 'D6A176' WHERE sample_barcode = 'GLSS-DK-0001-R1'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE biospecimen.aliquots SET aliquot_barcode = 'GLSS-DK-0003-R1-01D-WXS-89A5CD' , aliquot_uuid_short = '89A5CD' WHERE sample_barcode = 'GLSS-DK-0003-R1'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE biospecimen.aliquots SET aliquot_barcode = 'GLSS-DK-0009-TP-01D-WXS-291A6E' , aliquot_uuid_short = '291A6E' WHERE sample_barcode = 'GLSS-DK-0009-TP'")
dbClearResult(rs)

