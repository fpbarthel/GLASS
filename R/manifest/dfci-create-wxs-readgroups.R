##############################################
# Generate new barcodes and readgroups for Dana Farber Cancer Institute GLASS data.
# Updated: 2019.02.19
# Author: Kevin J.
##################################################

# Local working directory for this data preparation. 
mybasedir = "/Users/johnsk/Documents/Life-History/glass-batch2/"
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

# Anzhela's created GLASS barcodes for these newly received samples. She did not include samples that either had
# only multisector samples or were included in another dataset.
batch2_barcodes = read.xlsx("/Users/johnsk/Documents/Life-History/glass-batch2/Barcodes_20190219.xlsx")
batch2_barcodes = data.frame(lapply(batch2_barcodes, trimws), stringsAsFactors = FALSE)
batch2_surgery = read.xlsx("/Users/johnsk/Documents/Life-History/glass-batch2/clinic_surgery20190219.xlsx")
batch2_surgery = data.frame(lapply(batch2_surgery, trimws), stringsAsFactors = FALSE)
batch2_cases = read.xlsx("/Users/johnsk/Documents/Life-History/glass-batch2/clinic_caseAM20190219.xlsx")
batch2_cases = data.frame(lapply(batch2_cases, trimws), stringsAsFactors = FALSE)

# Read in file paths and calculated md5sums.
dfci_manifest <- read_tsv("/Users/johnsk/Documents/Life-History/glass-batch2/dfci/dfci-manifest.txt", col_names = TRUE)

# Read in calculated file size.
dfci_filesize <- read_tsv("/Users/johnsk/Documents/Life-History/glass-batch2/dfci/dfci_filesize.txt", col_names = FALSE)
colnames(dfci_filesize) <- c("file_size", "sample_name")
# Trim off the "./" in each sample name.
dfci_filesize$sample_name <- gsub("\\./", "", dfci_filesize$sample_name )

# It is necessary to manually review the readgroup file before inputting it here.
dfci_readgroups <- read_tsv("/Users/johnsk/Documents/Life-History/glass-batch2/dfci/dfci_readgroups.txt", col_names = FALSE)
colnames(dfci_readgroups) <- c("RG", "bam_rgid", "bam_rgsm", "bam_rglb", "bam_rgpl", "bam_rgpu", "bam_rgcn", "bam_rgdt")

# Remove RG identifier before colon.
drop_prefix = function(x) {gsub(".*:","", x) }
dfci_readgroups[2:7] <- lapply(dfci_readgroups[2:7], drop_prefix)
dfci_readgroups$bam_rgdt = gsub("DT:", "", dfci_readgroups$bam_rgdt)

# Determine the number of read groups per sample. Only TCGA has more than 5 readgroups.
read_group_totals = dfci_readgroups %>% 
  group_by(bam_rgsm) %>% 
  summarise(n = n())

# The readgroup ID was already provided for the DFCI samples.
# dfci_readgroups = dfci_readgroups %>% 
# mutate(glass_rg_id = paste0(substr(sapply(strsplit(bam_rgpu, "\\."), function(x) paste(x[1:2], collapse = ".")), 1 , 5), 
#                            substring(sapply(strsplit(bam_rgpu, "\\."), function(x) paste(x[1:2], collapse = ".")), nchar(sapply(strsplit(bam_rgpu, "\\."), function(x) paste(x[1:2], collapse = ".")))-1, nchar(sapply(strsplit(bam_rgpu, "\\."), function(x) paste(x[1:2], collapse = ".")))), "")) 

dfci_barcodes = batch2_barcodes %>% 
  filter(grepl("-DF-", sample_barcode)) %>% 
  inner_join(dfci_filesize, by = c("bam_name"="sample_name")) %>% 
  inner_join(dfci_manifest, by = "bam_name") %>%  
  mutate(file_format = "uBAM",
         merge_id = gsub(".bam", "", bam_name)) %>% 
  select(patient_barcode:sample_name, file_size:merge_id)
  

#### Files ######
dfci_files = dfci_barcodes %>% 
  select(aliquot_barcode = sample_barcode, file_path = bam_path, file_name = bam_name, file_size, file_md5sum = bam_md5, file_format) %>%
  distinct()

#### Case Sources ####
dfci_case_sources <-data.frame("DF","Dana Farber Cancer Institute")
colnames(dfci_case_sources)<- colnames(case_sources)

### Cases ####
dfci_cases = batch2_cases %>% 
  filter(grepl("DF", case_source)) 

### Samples ####
# Grab last two characrters of barcode. 
dfci_samples = dfci_barcodes %>% 
  mutate(sample_type = substring(sample_barcode, 14, 15),
         case_barcode =  substring(sample_barcode, 1, 12),
         sample_barcode = substring(sample_barcode, 1, 15)) %>% 
  select(case_barcode, sample_barcode, sample_type) %>% 
  distinct()

### Aliquots ####
dfci_aliquots = dfci_barcodes %>% 
  mutate(aliquot_barcode = sample_barcode,
         aliquot_id_legacy = sample_name,
         sample_barcode = substr(aliquot_barcode, 1, 15), 
         aliquot_uuid_short = substr(aliquot_barcode, 25, 30), 
         aliquot_analyte_type = substr(aliquot_barcode, 19, 19), 
         aliquot_analysis_type = substr(aliquot_barcode, 21, 23),
         aliquot_portion = as.integer(substr(aliquot_barcode, 17, 18)),
         aliquot_batch = paste(substr(aliquot_barcode, 1, 7), substr(aliquot_barcode, 21, 23), sep="-")) %>% 
  select(aliquot_barcode, aliquot_id_legacy, sample_barcode, aliquot_uuid_short, 
         aliquot_analyte_type, aliquot_analysis_type, aliquot_portion, aliquot_batch)
         
### Surgeries ####  
dfci_surgeries = batch2_surgery %>% 
  filter(grepl("-DF-", case_barcode)) 
  
### Readgroups ####
# Revised: Can't rename RGID, but pipeline is looking for old RGID.
dfci_readgroups_sub = dfci_readgroups %>% 
  inner_join(dfci_barcodes, by=c("bam_rgsm"="merge_id"))

readgroup_df = dfci_readgroups_sub %>% 
  mutate(readgroup_platform = bam_rgpl,
         legacy_readgroup_id = bam_rgid,
         readgroup_platform_unit = bam_rgpu,
         readgroup_date = bam_rgdt,
         readgroup_library = bam_rglb,
         readgroup_center = bam_rgcn,
         readgroup_sample_id = sample_barcode,
         readgroup_id = bam_rgid)

# Finalize readgroup information in predefined order.
dfci_readgroups = readgroup_df %>% 
  mutate(readgroup_platform = "ILLUMINA") %>% 
  select(aliquot_barcode = sample_barcode, readgroup_idtag=readgroup_id, readgroup_platform, readgroup_platform_unit, readgroup_library,
         readgroup_center, readgroup_idtag_legacy = legacy_readgroup_id) %>% 
  mutate(readgroup_sample_id = aliquot_barcode) %>% distinct()

## Filemap
tmp1 = files %>% select(file_name, aliquot_barcode)
tmp2 = dfci_readgroups %>% select(readgroup_idtag, readgroup_sample_id, aliquot_barcode)
dfci_files_readgroups = full_join(tmp1, tmp2) %>% select(file_name, readgroup_idtag, readgroup_sample_id)


### OUTPUT ####
# Only upload NEW cases. You can query the database and perform an antijoin.
extant_cases <- dbReadTable(con,  Id(schema="clinical",table="cases"))
new_cases = dfci_cases %>% 
  anti_join(extant_cases, by="case_barcode")
# New surgeries
extant_surgeries <- dbReadTable(con,  Id(schema="clinical",table="surgeries"))
new_surgeries = dfci_surgeries %>% 
  anti_join(extant_surgeries, by="case_barcode")
# New samples.
extant_samples <- dbReadTable(con,  Id(schema="biospecimen",table="samples"))
new_samples = dfci_samples %>% 
  anti_join(extant_samples, by="sample_barcode")
# New samples.
extant_aliquots <- dbReadTable(con,  Id(schema="biospecimen",table="aliquots"))
new_aliquots = dfci_aliquots %>% 
  anti_join(extant_aliquots, by="aliquot_barcode")
extant_aliquots <- dbReadTable(con,  Id(schema="biospecimen",table="aliquots"))
new_aliquots = dfci_aliquots %>% 
  anti_join(extant_aliquots, by="aliquot_barcode")


## Write to database.
dbWriteTable(con, Id(schema="clinical", table="case_sources"), dfci_case_sources, append=T)
dbWriteTable(con, Id(schema="clinical", table="cases"), new_cases, append=T)
dbWriteTable(con, Id(schema="clinical", table="surgeries"), new_surgeries, append=T)
dbWriteTable(con, Id(schema="biospecimen", table="samples"), new_samples, append=T)
dbWriteTable(con, Id(schema="biospecimen", table="aliquots"), new_aliquots, append=T)
dbWriteTable(con, Id(schema="biospecimen", table="readgroups"), dfci_readgroups, append=T)
dbWriteTable(con, Id(schema="analysis", table="files"), dfci_files, append=T)
dbWriteTable(con, Id(schema="analysis", table="files_readgroups"), dfci_files_readgroups, append=T)






