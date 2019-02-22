##############################################
# Generate new barcodes and readgroups for CW (Barnholtz et al.) GLASS data.
# Updated: 2019.02.22
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
#### NOTE:
### There were two patients GLSS-19-0275 (aka TEMP-159) and GLSS-19-0276 (aka TEMP-163) that do not have any normal bams.
### We already had GLSS-19-0265 (TCGA-19-0957) and GLSS-19-0281 (TCGA-FG-5965) in the dataset.
####

# Establish connection with the GLASS database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

# Load essential tables. This will be helpful so that you know which fields are required for each table.
cases = dbReadTable(con,  Id(schema="clinical", table="cases"))
case_sources = dbReadTable(con,  Id(schema="clinical", table="case_sources"))
sample_types = dbReadTable(con,  Id(schema="biospecimen", table="sample_types"))
aliquots = dbReadTable(con,  Id(schema="biospecimen", table="aliquots"))
samples = dbReadTable(con,  Id(schema="biospecimen", table="samples"))
readgroups = dbReadTable(con,  Id(schema="biospecimen", table="readgroups"))
files = dbReadTable(con,  Id(schema="analysis", table="files"))
files_readgroups = dbReadTable(con,  Id(schema="analysis", table="files_readgroups"))

# Read in the barcodes that Anzhela prepared:
batch2_barcodes = read.xlsx("/Users/johnsk/Documents/Life-History/glass-batch2/Barcodes_20190222.xlsx")
batch2_barcodes = data.frame(lapply(batch2_barcodes, trimws), stringsAsFactors = FALSE)
# We are missing surgery information for these -19- samples at this time.
batch2_surgery = read.xlsx("/Users/johnsk/Documents/Life-History/glass-batch2/clinic_surgery20190219.xlsx")
batch2_surgery = data.frame(lapply(batch2_surgery, trimws), stringsAsFactors = FALSE)
# There should be 13 new cases.
batch2_cases = read.xlsx("/Users/johnsk/Documents/Life-History/glass-batch2/clinic_caseAM20190219.xlsx")
batch2_cases = data.frame(lapply(batch2_cases, trimws), stringsAsFactors = FALSE)

# Restrict the barcodes to just the CW samples.
cwru_barcodes = batch2_barcodes %>% 
  filter(grepl("-19-", sample_barcode), !grepl("gbm*", sample_name), grepl("sorted", bam_name)) %>% 
  select(patient_barcode:Column1) %>% 
  mutate(sample_type = substr(surgery_barcode, 14, 15),
         rgID = ifelse(sample_type=="NB", "NORMAL", "TUMOR"))

# Read in file paths and md5sums that were calculated on helix.
cwru_bam_md5 = read.delim("/Users/johnsk/Documents/Life-History/glass-batch2/cwru/case_western_bam_md5.txt", sep = "", header = FALSE, stringsAsFactors = FALSE)
colnames(cwru_bam_md5) <- c("bam_md5", "file_path")

# Determine file size for the files on helix.
cwru_files = read.delim("/Users/johnsk/Documents/Life-History/glass-batch2/cwru/case_western_bams.txt", sep = "", header = FALSE, stringsAsFactors = FALSE)
colnames(cwru_files) = "file_path"
last <- function(x) { return( x[length(x)] ) }
cwru_files = cwru_files %>% 
  mutate(local_path = gsub("/projects/verhaak-lab/", "/Volumes/verhaak-lab/", file_path), 
         file_size = file.info(local_path)$size,
         file_name = sapply(strsplit(file_path, "/"), last), 
         file_format = "uBAM")

# This particular analysis is performed without the TCGA normal samples.
cwru_df = cwru_files %>% 
  inner_join(cwru_bam_md5, by="file_path") %>% 
  inner_join(cwru_barcodes, by = c("file_name" = "bam_name"))

#### Files ######
files = cwru_df %>% 
  select(aliquot_barcode = sample_barcode, file_path, file_name, file_size, file_md5sum = bam_md5, file_format) %>%
  distinct()

### Cases ####
cwru_cases = batch2_cases %>% 
  filter(grepl("GLSS-19", case_barcode)) 

### Samples ####
cwru_samples = cwru_df %>% 
  mutate(sample_type = substring(sample_barcode, 14, 15),
         case_barcode =  substring(sample_barcode, 1, 12),
         sample_barcode = substring(sample_barcode, 1, 15)) %>% 
  select(case_barcode, sample_barcode, sample_type) %>% 
  distinct()

### Surgeries ####  
cwru_surgeries = batch2_surgery %>% 
  filter(grepl("GLSS-19", case_barcode)) 

### Aliquots ####
cwru_aliquots = cwru_df %>% 
  mutate(aliquot_barcode = sample_barcode,
         aliquot_id_legacy = sample_name,
         sample_barcode = substr(aliquot_barcode, 1, 15), 
         aliquot_uuid_short = substr(aliquot_barcode, 25, 30), 
         aliquot_analyte_type = substr(aliquot_barcode, 19, 19), 
         aliquot_analysis_type = substr(aliquot_barcode, 21, 23),
         aliquot_portion = as.integer(substr(aliquot_barcode, 17, 18)),
         aliquot_batch = paste(substr(aliquot_barcode, 1, 7), substr(aliquot_barcode, 21, 23), sep="-")) %>% 
  select(aliquot_barcode, aliquot_id_legacy, sample_barcode, aliquot_uuid_short, 
         aliquot_analyte_type, aliquot_analysis_type, aliquot_portion, aliquot_batch) %>% 
  distinct()


# It is necessary to manually review the readgroup file before inputting it here.
cwru_readgroups <- read_tsv("/Users/johnsk/Documents/Life-History/glass-batch2/cwru/case_western_readgroups.txt", col_names = FALSE)
colnames(cwru_readgroups) <- c("RG", "bam_rgid", "bam_rgsm", "bam_rgpl","bam_rglb")

# Remove RG identifier before colon.
drop_prefix = function(x) {gsub(".*:","", x) }
cwru_readgroups[2:5] <- lapply(cwru_readgroups[2:5], drop_prefix)
# Restrict to the only two options (TUMOR and NORMAL).
cwru_readgroups = cwru_readgroups %>% distinct()

# Define the readgroups.
cwru_only_readgroups = cwru_barcodes %>% 
  left_join(cwru_readgroups, by=c("rgID" = "bam_rgid"))

# Separate so that it's easier to follow the renaming scheme.
readgroup_df = cwru_only_readgroups %>% 
  mutate(readgroup_platform = bam_rgpl,
         legacy_readgroup_id = rgID,
         readgroup_platform_unit = "UNKNOWN",
         readgroup_date = NA,
         readgroup_library = bam_rglb,
         readgroup_center = "CWRU",
         readgroup_sample_id = sample_barcode,
         readgroup_id = "UNKWN.1")

# Finalize readgroup information in predefined order.
cwru_only_readgroups = readgroup_df %>% 
  mutate(readgroup_platform = "ILLUMINA") %>% 
  select(aliquot_barcode = sample_barcode, readgroup_idtag=readgroup_id, readgroup_idtag_legacy = legacy_readgroup_id, readgroup_platform, readgroup_platform_unit, 
         readgroup_library, readgroup_center) %>% 
  mutate(readgroup_sample_id = aliquot_barcode) %>% distinct()

# Create a files_readgroup file.
tmp1 = files %>% select(file_name, aliquot_barcode) %>% distinct()
tmp2 = cwru_only_readgroups %>% select(readgroup_idtag, readgroup_sample_id, aliquot_barcode) %>% distinct()
cwru_files_readgroups = full_join(tmp1, tmp2) %>% select(file_name, readgroup_idtag, readgroup_sample_id)

### OUTPUT ####
# Only upload NEW cases. You can query the database and perform an antijoin.
extant_cases <- dbReadTable(con,  Id(schema="clinical",table="cases"))
new_cases = cwru_cases %>% 
  anti_join(extant_cases, by="case_barcode")
# New surgeries
extant_surgeries <- dbReadTable(con,  Id(schema="clinical",table="surgeries"))
new_surgeries = cwru_surgeries %>% 
  anti_join(extant_surgeries, by="case_barcode")
# New samples.
extant_samples <- dbReadTable(con,  Id(schema="biospecimen",table="samples"))
new_samples = cwru_samples %>% 
  anti_join(extant_samples, by="sample_barcode")
# New samples.
extant_aliquots <- dbReadTable(con,  Id(schema="biospecimen",table="aliquots"))
new_aliquots = cwru_aliquots %>% 
  anti_join(extant_aliquots, by="aliquot_barcode")
# New readgroups.
extant_readgroups <- dbReadTable(con,  Id(schema="biospecimen", table="readgroups"))
new_readgroups = cwru_only_readgroups %>% 
  anti_join(extant_readgroups, by="aliquot_barcode")
# New files
extant_files <- dbReadTable(con,  Id(schema="analysis", table="files"))
new_files = files %>% 
  anti_join(extant_files, by="aliquot_barcode")
# New files_readgroups.
extant_files_readgroups <- dbReadTable(con,  Id(schema="analysis", table="files_readgroups"))
new_files_readgroups = cwru_files_readgroups %>% 
  anti_join(extant_files_readgroups, by="readgroup_sample_id")

########################
####### OUTPUT #########
########################
dbWriteTable(con, Id(schema="clinical", table="cases"), new_cases, append=T)
dbWriteTable(con, Id(schema="biospecimen", table="samples"), new_samples, append=T)
dbWriteTable(con, Id(schema="clinical", table="surgeries"), new_surgeries, append=T)
dbWriteTable(con, Id(schema="biospecimen", table="aliquots"), new_aliquots, append=T)
dbWriteTable(con, Id(schema="biospecimen", table="readgroups"), new_readgroups, append=T)
dbWriteTable(con, Id(schema="analysis", table="files"), new_files, append=T)
dbWriteTable(con, Id(schema="analysis", table="files_readgroups"), new_files_readgroups, append=T)



######## TCGA normals in CW data #########
# For those samples with TCGA normals, keep the GLSS ID Anzhela has designated, BUT
# remove the aliquot_batch.
# Restrict the barcodes to just the TCGA samples also found in the CWRU dataset.
cw_tcga_barcodes = batch2_barcodes %>% 
  filter(grepl("-19-", sample_barcode), grepl("TCGA", bam_name)) %>% 
  select(patient_barcode:Column1) %>% 
  mutate(sample_type = substr(surgery_barcode, 14, 15),
         rgID = ifelse(sample_type=="NB", "NORMAL", "TUMOR"))

# Read in file paths and calculated md5sums.
cwru_tcga_md5 = read.delim("/Users/johnsk/Documents/Life-History/glass-batch2/cwru/case_western_tcga_md5.txt", sep = "", header = FALSE, stringsAsFactors = FALSE)
colnames(cwru_tcga_md5) <- c("bam_md5", "file_path")

# Determine file size.
cwru_tcga_files = read.delim("/Users/johnsk/Documents/Life-History/glass-batch2/cwru/case_western_tcga_bams.txt", sep = "", header = FALSE, stringsAsFactors = FALSE)
colnames(cwru_tcga_files) = "file_path"
last <- function(x) { return( x[length(x)] ) }
cwru_tcga_files = cwru_tcga_files %>% 
  mutate(local_path = gsub("/projects/verhaak-lab/", "/Volumes/verhaak-lab/", file_path), 
         file_size = file.info(local_path)$size,
         file_name = sapply(strsplit(file_path, "/"), last), 
         file_format = "uBAM")

# This analysis is performed without the TCGA samples.
cwru_tcga_df = cwru_tcga_files %>% 
  inner_join(cwru_tcga_md5, by="file_path") %>% 
  inner_join(cw_tcga_barcodes, by = c("file_name" = "bam_name"))

# Files
cw_tcga_files = cwru_tcga_df %>% 
  select(aliquot_barcode = sample_barcode, file_path, file_name, file_size, file_md5sum = bam_md5, file_format) %>%
  distinct()

### Samples ####
# Grab last two characrters of barcode. 
cw_tcga_samples = cwru_tcga_df %>% 
  mutate(sample_type = substring(sample_barcode, 14, 15),
         case_barcode =  substring(sample_barcode, 1, 12),
         sample_barcode = substring(sample_barcode, 1, 15)) %>% 
  select(case_barcode, sample_barcode, sample_type) %>% 
  distinct()

### Aliquots ####
cw_tcga_aliquots = cwru_tcga_df %>% 
  mutate(aliquot_barcode = sample_barcode,
         aliquot_id_legacy = sample_name,
         sample_barcode = substr(aliquot_barcode, 1, 15), 
         aliquot_uuid_short = substr(aliquot_barcode, 25, 30), 
         aliquot_analyte_type = substr(aliquot_barcode, 19, 19), 
         aliquot_analysis_type = substr(aliquot_barcode, 21, 23),
         aliquot_portion = as.integer(substr(aliquot_barcode, 17, 18)),
         aliquot_batch = NA) %>% 
  select(aliquot_barcode, aliquot_id_legacy, sample_barcode, aliquot_uuid_short, 
         aliquot_analyte_type, aliquot_analysis_type, aliquot_portion, aliquot_batch)

# It is necessary to manually review the readgroup file before inputting it here.
cw_tcga_readgroups <- read_tsv("/Users/johnsk/Documents/Life-History/glass-batch2/cwru/case_western_tcga_readgroups.txt", col_names = FALSE)
colnames(cw_tcga_readgroups) <- c("RG", "bam_rgid", "bam_rgpl", "bam_rgpu", "bam_rglb", "bam_rgdt", "bam_rgsm", "bam_rgcn")
cw_tcga_readgroups = cw_tcga_readgroups %>% 
  select(RG, bam_rgid, bam_rgpl, bam_rgpu, bam_rglb, bam_rgsm, bam_rgcn, bam_rgdt)

# Remove RG identifier before colon.
drop_prefix = function(x) {gsub(".*:","", x) }
cw_tcga_readgroups[2:7] <- lapply(cw_tcga_readgroups[2:7], drop_prefix)
cw_tcga_readgroups$bam_rgdt = gsub("DT:", "", cw_tcga_readgroups$bam_rgdt)

# Restrict to the only two options (TUMOR and NORMAL).
cw_tcga_readgroups = cw_tcga_readgroups %>% distinct()
cw_tcga_readgroups$mergeID = substr(cw_tcga_readgroups$bam_rgsm, 1, 25)

# Define the readgroups.
tcga_only_readgroups = cwru_tcga_df %>% 
  mutate(mergeID = substr(file_name, 1, 25)) %>% 
  left_join(cw_tcga_readgroups, by= "mergeID")

# Recreate the appropriate fields.
readgroup_df = tcga_only_readgroups %>% 
  mutate(readgroup_platform = bam_rgpl,
         legacy_readgroup_id = bam_rgid,
         readgroup_platform_unit = bam_rgpu,
         readgroup_date = bam_rgdt,
         readgroup_library = bam_rglb,
         readgroup_center = bam_rgcn,
         readgroup_sample_id = sample_barcode,
         readgroup_id = c("UNKWN.1","UNKWN.1", "UNKWN.2"))

# Finalize readgroup information in predefined order.
tcga_only_readgroups = readgroup_df %>% 
  mutate(readgroup_platform = "SOLiD") %>% 
  select(aliquot_barcode = sample_barcode, readgroup_idtag=readgroup_id, readgroup_platform, readgroup_platform_unit, readgroup_library,
         readgroup_center, readgroup_idtag_legacy = legacy_readgroup_id) %>% 
  mutate(readgroup_sample_id = aliquot_barcode) %>% distinct()

## file_readgroup order.
tmp1 = cw_tcga_files %>% select(file_name, aliquot_barcode) %>% distinct()
tmp2 = tcga_only_readgroups %>% select(readgroup_idtag, readgroup_sample_id, aliquot_barcode) %>% distinct()
tcga_files_readgroups = full_join(tmp1, tmp2) %>% select(file_name, readgroup_idtag, readgroup_sample_id)

########################
####### OUTPUT #########
########################
## Write to database.
dbWriteTable(con, Id(schema="biospecimen", table="samples"), cw_tcga_samples, append=T)
dbWriteTable(con, Id(schema="biospecimen", table="aliquots"), cw_tcga_aliquots, append=T)
dbWriteTable(con, Id(schema="biospecimen", table="readgroups"), tcga_only_readgroups, append=T)
dbWriteTable(con, Id(schema="analysis", table="files"), cw_tcga_files, append=T)
dbWriteTable(con, Id(schema="analysis", table="files_readgroups"), tcga_files_readgroups, append=T)

# I neglected to include the aliquot_legacy_id for the two TCGA normal samples.
rs = dbSendStatement(con, "UPDATE biospecimen.aliquots SET aliquot_id_legacy = 'TCGA-19-0955-10A-01W-0609-10' WHERE aliquot_barcode = 'GLSS-19-0266-NB-01D-WXS-GGPQS4'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE biospecimen.aliquots SET aliquot_id_legacy = 'TCGA-19-0963-10A-01W-0609-10' WHERE aliquot_barcode = 'GLSS-19-0267-NB-01D-WXS-RXL9E0'")
dbClearResult(rs)


####### DFCI in CW data ##########
## *** Replace ALIQUOT_BATCH FOR DFCI *** ##
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

# Determine those patients that part of the Dana Farber data set that also listed in the Case Western data set.
dfci_cw_barcodes = batch2_barcodes %>% 
  filter(grepl("-19-", sample_barcode), grepl("gbm*", sample_name)) %>% 
  inner_join(dfci_filesize, by = c("bam_name"="sample_name")) %>% 
  inner_join(dfci_manifest, by = "bam_name") %>%  
  mutate(file_format = "uBAM",
         merge_id = gsub(".bam", "", bam_name)) %>% 
  select(patient_barcode:sample_name, file_size:merge_id)

# Define the files to be uploaded for Case Western samples sequenced at Dana Farber.
dfci_cw_files = dfci_cw_barcodes %>% 
  select(aliquot_barcode = sample_barcode, file_path = bam_path, file_name = bam_name, file_size, file_md5sum = bam_md5, file_format) %>%
  distinct()

# Grab last two characrters of barcode. 
dfci_cw_samples = dfci_cw_barcodes %>% 
  mutate(sample_type = substring(sample_barcode, 14, 15),
         case_barcode =  substring(sample_barcode, 1, 12),
         sample_barcode = substring(sample_barcode, 1, 15)) %>% 
  select(case_barcode, sample_barcode, sample_type) %>% 
  distinct()

### Aliquots ####
dfci_cw_aliquots = dfci_cw_barcodes %>% 
  mutate(aliquot_barcode = sample_barcode,
         aliquot_id_legacy = sample_name,
         sample_barcode = substr(aliquot_barcode, 1, 15), 
         aliquot_uuid_short = substr(aliquot_barcode, 25, 30), 
         aliquot_analyte_type = substr(aliquot_barcode, 19, 19), 
         aliquot_analysis_type = substr(aliquot_barcode, 21, 23),
         aliquot_portion = as.integer(substr(aliquot_barcode, 17, 18)),
         aliquot_batch = "GLSS-DF-WXS") %>% 
  select(aliquot_barcode, aliquot_id_legacy, sample_barcode, aliquot_uuid_short, 
         aliquot_analyte_type, aliquot_analysis_type, aliquot_portion, aliquot_batch)


# Define the readgroups.
dfci_cw_readgroups = dfci_readgroups %>% 
  inner_join(dfci_cw_barcodes, by=c("bam_rgsm"="merge_id"))
# Prepare the data with the appropriate fields.
readgroup_df = dfci_cw_readgroups %>% 
  mutate(readgroup_platform = bam_rgpl,
         legacy_readgroup_id = bam_rgid,
         readgroup_platform_unit = bam_rgpu,
         readgroup_date = bam_rgdt,
         readgroup_library = bam_rglb,
         readgroup_center = bam_rgcn,
         readgroup_sample_id = sample_barcode,
         readgroup_id = bam_rgid)

# Finalize readgroup information in predefined order.
dfci_cw_readgroups = readgroup_df %>% 
  mutate(readgroup_platform = "ILLUMINA") %>% 
  select(aliquot_barcode = sample_barcode, readgroup_idtag=readgroup_id, readgroup_platform, readgroup_platform_unit, readgroup_library,
         readgroup_center, readgroup_idtag_legacy = legacy_readgroup_id) %>% 
  mutate(readgroup_sample_id = aliquot_barcode) %>% distinct()

## Filemap
tmp1 = dfci_cw_files %>% select(file_name, aliquot_barcode) %>% distinct()
tmp2 = dfci_cw_readgroups %>% select(readgroup_idtag, readgroup_sample_id, aliquot_barcode) %>% distinct()
dfci_cw_files_readgroups = full_join(tmp1, tmp2) %>% select(file_name, readgroup_idtag, readgroup_sample_id)

########################
####### OUTPUT #########
########################
dbWriteTable(con, Id(schema="biospecimen", table="aliquots"), dfci_cw_aliquots, append=T)
dbWriteTable(con, Id(schema="biospecimen", table="readgroups"), dfci_cw_readgroups, append=T)
dbWriteTable(con, Id(schema="analysis", table="files"), dfci_cw_files, append=T)
dbWriteTable(con, Id(schema="analysis", table="files_readgroups"), dfci_cw_files_readgroups, append=T)

#### Double-check that all UUIDs are correct ####
extant_aliquots <- dbReadTable(con,  Id(schema="biospecimen",table="aliquots"))
dim(extant_aliquots) # 1042.
n_distinct(extant_aliquots$aliquot_uuid_short) # 1042.

extant_samples <- dbReadTable(con,  Id(schema="biospecimen",table="samples"))
tmp1 = extant_samples %>% 
  filter(grepl('GLSS-19-', sample_barcode))
n_distinct(tmp1$case_barcode)  
