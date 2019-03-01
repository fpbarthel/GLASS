##############################################
# Generate new barcodes and readgroups for MGH (Brastianos et al.) GLASS data.
# Updated: 2019.02.21
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
sample_types = dbReadTable(con,  Id(schema="biospecimen", table="sample_types"))
aliquots = dbReadTable(con,  Id(schema="biospecimen", table="aliquots"))
samples = dbReadTable(con,  Id(schema="biospecimen", table="samples"))
readgroups = dbReadTable(con,  Id(schema="biospecimen", table="readgroups"))
files = dbReadTable(con,  Id(schema="analysis", table="files"))
files_readgroups = dbReadTable(con,  Id(schema="analysis", table="files_readgroups"))

# Read in the barcodes that Anzhela prepared:
batch2_barcodes = read.xlsx("/Users/johnsk/Documents/Life-History/glass-batch2/Barcodes_20190219.xlsx")
batch2_barcodes = data.frame(lapply(batch2_barcodes, trimws), stringsAsFactors = FALSE)
# We are missing surgery information for these -MG- samples at this time.
batch2_surgery = read.xlsx("/Users/johnsk/Documents/Life-History/glass-batch2/clinic_surgery20190219.xlsx")
batch2_surgery = data.frame(lapply(batch2_surgery, trimws), stringsAsFactors = FALSE)
# There should be 12 new cases.
batch2_cases = read.xlsx("/Users/johnsk/Documents/Life-History/glass-batch2/clinic_caseAM20190219.xlsx")
batch2_cases = data.frame(lapply(batch2_cases, trimws), stringsAsFactors = FALSE)

# Restrict the barcodes to just the MGH samples.
mgh_barcodes = batch2_barcodes %>% 
  filter(grepl("-MG-", sample_barcode)) %>% 
  select(patient_barcode:Column1)

# Load in the manifest file from SRA.
mgh_manifest = read.delim("/Volumes/verhaak-lab/kimh_BRASTIANOS_NPJ/wes/dbgap/manifest/SraRunTable.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# Confirm the samples with lung tissue listed: c("GLSS-MG-0008", "GLSS-MG-0014").
mgh_manifest[mgh_manifest$body_site=="Lung", ]

# Access information about readgroups and md5sums.
mgh_md5_checksums = read.delim("/Users/johnsk/Documents/Life-History/glass-batch2/mgh/mgh_brastianos_checksums.txt", sep = "", header = FALSE, stringsAsFactors = FALSE)
colnames(mgh_md5_checksums) <- c("bam_md5", "file_path")

# Curate the readgroup information to be passed to the alignment pipeline.
mgh_readgroups = read.delim("/Users/johnsk/Documents/Life-History/glass-batch2/mgh/mgh_brastianos_readgroups.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(mgh_readgroups) <- c("RG", "bam_rgid","bam_rgpl", "bam_rgpu", "bam_rglb", "bam_pi", "bam_rgdt", "bam_rgsm", "bam_rgcn")
mgh_readgroups = mgh_readgroups %>% 
  select(RG, bam_rgid, bam_rgpl, bam_rgpu, bam_rglb, bam_pi, bam_rgsm, bam_rgcn, bam_rgdt)

# Remove RG identifier before colon.
drop_prefix = function(x) {gsub(".*:","", x) }
mgh_readgroups[2:8] <- lapply(mgh_readgroups[2:8], drop_prefix)
mgh_readgroups$bam_rgdt = gsub("DT:", "", mgh_readgroups$bam_rgdt)

# Determine the number of read groups per sample. Only TCGA has more than 5 readgroups.
read_group_totals = mgh_readgroups %>% 
  group_by(bam_rgsm) %>% 
  summarise(n = n())

# Retrieve the information about the files, including path and size.
bam_files = list.files('/Volumes/verhaak-lab/kimh_BRASTIANOS_NPJ/wes/dbgap/bams', pattern = ".bam$", full.names = T)
files = data.frame(aliquot_id_legacy = gsub(".sra.bam", "", basename(bam_files)),
                   file_name = basename(bam_files),
                   file_path = gsub("/Volumes/", "/projects/", bam_files),
                   file_size = file.info(bam_files)$size,
                   file_format = "uBAM",
                   stringsAsFactors = F)

# Combine all metadata information in a single dataframe,
mgh_df = mgh_barcodes %>% 
  inner_join(mgh_manifest, by=c("sample_name" = "Sample_Name")) %>% 
  inner_join(files, by=c("Run"= "aliquot_id_legacy")) %>% 
  inner_join(mgh_md5_checksums, by="file_path")

#### Files ######
mgh_files = mgh_df %>% 
  select(aliquot_barcode = sample_barcode, file_path, file_name, file_size, file_md5sum = bam_md5, file_format) %>%
  distinct()

#### Case Sources ####
mgh_case_sources <-data.frame("MG","Massachusetts General Hospital")
colnames(mgh_case_sources)<- colnames(case_sources)

### Cases ####
mgh_cases = batch2_cases %>% 
  filter(grepl("MG", case_source)) 

### Samples ####
# Grab last two characrters of barcode. 
mgh_samples = mgh_barcodes %>% 
  mutate(sample_type = substring(sample_barcode, 14, 15),
         case_barcode =  substring(sample_barcode, 1, 12),
         sample_barcode = substring(sample_barcode, 1, 15)) %>% 
  select(case_barcode, sample_barcode, sample_type) %>% 
  distinct()

### Aliquots ####
mgh_aliquots = mgh_barcodes %>% 
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
mgh_surgeries = batch2_surgery %>% 
  filter(grepl("-MG-", case_barcode)) 

#### Readgroups #####
readgroup_df = mgh_readgroups %>% 
  inner_join(mgh_df, by=c("bam_rgsm" = "biospecimen_repository_sample_id")) %>% 
  mutate(readgroup_platform = bam_rgpl,
         legacy_readgroup_id = bam_rgid,
         readgroup_platform_unit = bam_rgpu,
         readgroup_date = bam_rgdt,
         readgroup_library = bam_rglb,
         readgroup_center = bam_rgcn,
         readgroup_sample_id = sample_barcode,
         readgroup_id = bam_rgid)

# Finalize readgroup information in predefined order.
mgh_readgroups = readgroup_df %>% 
  mutate(readgroup_platform = "ILLUMINA") %>% 
  select(aliquot_barcode = sample_barcode, readgroup_idtag=readgroup_id, readgroup_platform, readgroup_platform_unit, readgroup_library,
         readgroup_center, readgroup_idtag_legacy = legacy_readgroup_id) %>% 
  mutate(readgroup_sample_id = aliquot_barcode) %>% distinct()

## Filemap ##
tmp1 = mgh_files %>% select(file_name, aliquot_barcode) %>% distinct()
tmp2 = mgh_readgroups %>% select(readgroup_idtag, readgroup_sample_id, aliquot_barcode) %>% distinct()
mgh_files_readgroups = full_join(tmp1, tmp2) %>% select(file_name, readgroup_idtag, readgroup_sample_id)

# Only upload NEW cases. You can query the database and perform an antijoin.
extant_cases <- dbReadTable(con,  Id(schema="clinical",table="cases"))
new_cases = mgh_cases %>% 
  anti_join(extant_cases, by="case_barcode")
# New surgeries
extant_surgeries <- dbReadTable(con,  Id(schema="clinical",table="surgeries"))
new_surgeries = mgh_surgeries %>% 
  anti_join(extant_surgeries, by="case_barcode")
# New samples.
extant_samples <- dbReadTable(con,  Id(schema="biospecimen",table="samples"))
new_samples = mgh_samples %>% 
  anti_join(extant_samples, by="sample_barcode")
# New samples.
extant_aliquots <- dbReadTable(con,  Id(schema="biospecimen",table="aliquots"))
new_aliquots = mgh_aliquots %>% 
  anti_join(extant_aliquots, by="aliquot_barcode")
# New readgroups.
extant_readgroups <- dbReadTable(con,  Id(schema="biospecimen", table="readgroups"))
new_readgroups = mgh_readgroups %>% 
  anti_join(extant_readgroups, by="aliquot_barcode")
# New files
extant_files <- dbReadTable(con,  Id(schema="analysis", table="files"))
new_files = mgh_files %>% 
  anti_join(extant_files, by="aliquot_barcode")
# New files_readgroups.
extant_files_readgroups <- dbReadTable(con,  Id(schema="analysis", table="files_readgroups"))
new_files_readgroups = mgh_files_readgroups %>% 
  anti_join(extant_files_readgroups, by="readgroup_sample_id")

## Write to database.
dbWriteTable(con, Id(schema="clinical", table="case_sources"), mgh_case_sources, append=T)
dbWriteTable(con, Id(schema="clinical", table="cases"), new_cases, append=T)
dbWriteTable(con, Id(schema="biospecimen", table="samples"), new_samples, append=T)
dbWriteTable(con, Id(schema="clinical", table="surgeries"), new_surgeries, append=T)
dbWriteTable(con, Id(schema="biospecimen", table="aliquots"), new_aliquots, append=T)
dbWriteTable(con, Id(schema="biospecimen", table="readgroups"), new_readgroups, append=T)
dbWriteTable(con, Id(schema="analysis", table="files"), new_files, append=T)
dbWriteTable(con, Id(schema="analysis", table="files_readgroups"), new_files_readgroups, append=T)

## After initial processing of the data, a mistake was noted.
## Need to update -NT- to be -NM-; remove -NT- from sample types table; edit sample_types table.
rs = dbSendStatement(con, "DELETE FROM biospecimen.samples WHERE sample_barcode = 'GLSS-MG-0001-NT'")
dbClearResult(rs)
rs = dbSendStatement(con, "DELETE FROM biospecimen.samples WHERE sample_barcode = 'GLSS-MG-0002-NT'")
dbClearResult(rs)
rs = dbSendStatement(con, "DELETE FROM biospecimen.samples WHERE sample_barcode = 'GLSS-MG-0005-NT'")
dbClearResult(rs)
rs = dbSendStatement(con, "DELETE FROM biospecimen.samples WHERE sample_barcode = 'GLSS-MG-0008-NT'")
dbClearResult(rs)
rs = dbSendStatement(con, "DELETE FROM biospecimen.samples WHERE sample_barcode = 'GLSS-MG-0009-NT'")
dbClearResult(rs)
rs = dbSendStatement(con, "DELETE FROM biospecimen.samples WHERE sample_barcode = 'GLSS-MG-0010-NT'")
dbClearResult(rs)
rs = dbSendStatement(con, "DELETE FROM biospecimen.samples WHERE sample_barcode = 'GLSS-MG-0011-NT'")
dbClearResult(rs)
rs = dbSendStatement(con, "DELETE FROM biospecimen.samples WHERE sample_barcode = 'GLSS-MG-0014-NT'")
dbClearResult(rs)

# Update Normal Tissue classification.
rs = dbSendStatement(con, "DELETE FROM biospecimen.sample_types WHERE sample_type = 'NT'")
rs = dbSendStatement(con, "UPDATE biospecimen.sample_types SET sample_type_description = 'Normal tissue' WHERE sample_type = 'NM'")
dbClearResult(rs)



