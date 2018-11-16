#######################################################
# Examine some of the clinical data prepared by Anzhela
# Date: 2018.11.15 
# Author: Kevin J.
#######################################################

# Load the most up-to-date clinical table produce by Anzhela on GLASS github.
latest_clinical_data = '/Users/johnsk/Documents/Life-History/SurgeriesTableAM-201810241443-20181114.xlsx'

# Here is the codel status from the file I created from GATK copy number output (see `determine-codel-status.R`).
glass_codel_subject = "/Users/johnsk/Documents/Life-History/glass-subject-level-codel-status.txt"

#######################################################

# Necessary packages:
library(tidyverse)
library(GenomicRanges)
library(openxlsx)
library(DBI)

#######################################################
# Establish connection with the database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

## Retrieve the colnames() for the database to inform categories:
surgeries = dbReadTable(con,  Id(schema="clinical",table="surgeries"))
samples =  dbReadTable(con,  Id(schema="biospecimen",table="samples"))

## Load in clinical data, do not remove any samples with missing information.
glass_clinical = readWorkbook(latest_clinical_data, sheet = 1, startRow = 1, colNames = TRUE)
# Renaming the codel and idh columns to merge with database and copy number calls.
glass_clin = glass_clinical %>% 
  mutate(idh_status_old = idh_status,
         codel_status_old = codel_status) %>% 
  select(-idh_status, -codel_status)

## Make sure there is no trailing whitespace introduced by excel clinical sheet.
glass_clin$sample_barcode = trimws(glass_clin$sample_barcode, "r")

## Query the database to retrieve IDH1/IDH2 calls. Request specific variants.
glass_idh_muts <- dbGetQuery(con, "SELECT gt.aliquot_barcode, gene_symbol, hgvs_p, gt.called, gt.read_depth, gt.alt_count,
                             ROUND(CAST(gt.alt_count AS DECIMAL)/CAST(gt.read_depth AS DECIMAL),2) AS vaf,
                             (ROUND(CAST(gt.alt_count AS DECIMAL)/CAST(gt.read_depth AS DECIMAL),2) > 0.1 and read_depth > 3) AS manual_call
                             FROM analysis.snvs v, analysis.snv_genotypes gt
                             WHERE gt.chrom = v.chrom AND gt.start = v.start AND gt.end = v.end AND gt.alt = v.alt AND gene_symbol IN ('IDH1','IDH2') AND hgvs_p IN ('p.R132H','p.R132C','p.R132G','p.R132S','p.R172K')")

### Supervised approach for clinical IDH status.
# Use "manual_call" to increase sensitivity of mutation detection. Sample-level mutation status.
# This includes reverses the "manual_call" mentioned in GLASS-WE issue #18. Paired tumor shares same mutation.
glass_idh_muts$manual_call[glass_idh_muts$aliquot_barcode=="GLSS-MD-0093-R1-01D-WXS-PR89I0" & glass_idh_muts$hgvs_p=="p.R132H"] = "1"
glass_idh_muts$manual_call[glass_idh_muts$aliquot_barcode=="GLSS-MD-0084-R1-01D-WXS-O09B9L" & glass_idh_muts$hgvs_p=="p.R132H"] = "1"
glass_idh_muts$manual_call[glass_idh_muts$aliquot_barcode=="GLSS-SU-0239-R1-01D-WXS-7U2QLD" & glass_idh_muts$hgvs_p=="p.R132H"] = "1"
glass_idh_muts$manual_call[glass_idh_muts$aliquot_barcode=="GLSS-MD-0137-TP-01D-WXS-9LXGRU" & glass_idh_muts$hgvs_p=="p.R132H"] = "1"


## Include both Mutect2 and vaf threshold calls.
idh_mutants = glass_idh_muts %>% 
  mutate(consistent = called==manual_call) %>% 
  filter(called=="1"| consistent=="FALSE") %>% 
  mutate(sample_id = substr(aliquot_barcode, 1, 15)) %>% 
  distinct(sample_id, gene_symbol, hgvs_p) 

## These copy number calls were made using 1 kb bins.
glass_codel = read.delim(glass_codel_subject, as.is=T, header=T)

## This will generate a list of tumor samples that are 1p19q codeleted.
codel_tumor_samples = samples %>% 
  filter(!sample_type%in%c("NB", "NM")) %>% 
  left_join(glass_codel, by=c("case_barcode"="subject_id")) %>% 
  # Poor quality samples initially called co-deleted.
  filter(!case_barcode%in%c("GLSS-RV-LP07","GLSS-SF-0006")) 

#### Combine data to link with database. 
clin_idh_codel = left_join(glass_clin, idh_mutants, by=c("sample_barcode"="sample_id")) %>% 
  filter(!is.na(case_barcode)) %>% 
  left_join(codel_tumor_samples, by="sample_barcode") %>% 
  mutate(idh_status = ifelse(is.na(hgvs_p), "IDHwt", "IDHmut"),
         codel_status = ifelse(is.na(status_1p19q), "noncodel", status_1p19q),
         idh_codel_subtype = paste(idh_status, codel_status, sep="_")) %>%
  # Provide all fields for the database.
  select(case_barcode = case_barcode.x, surgery_number, sample_barcode, surgical_interval_mo, histology, grade, idh_status, codel_status, who_classification,
       mgmt_methylation, surgery_type, surgery_indication, surgery_extent_of_resection, surgery_laterality, surgery_location, treatment_tmz, treatment_tmz_cycles, treatment_tmz_cycles_6 = treatment_tmz_cycles_.6, treatment_tmz_cycles_notes, treatment_concurrent_tmz, 
         treatment_radiotherapy, treatment_radiation_dose_gy, treatment_radiation_fractions, treatment_radiation_other, treatment_chemotherapy_other,
         treatment_chemotherapy_other_cycles, comments, idh_codel_subtype) %>% 
  mutate(surgical_interval_mo = ceiling(as.numeric(surgical_interval_mo))) %>% 
  # Remove one particular sample.
  filter(case_barcode!="GLSS-MD-0134")

# Remove any molecular designation for samples that have not been sequenced.
# We do not have clinical data to confirm subtype in absence of molecular.
clin_idh_codel$idh_status[is.na(clin_idh_codel$sample_barcode)] = NA
clin_idh_codel$codel_status[is.na(clin_idh_codel$sample_barcode)] = NA
clin_idh_codel$idh_codel_subtype[is.na(clin_idh_codel$sample_barcode)] = NA
clin_idh_codel$idh_codel_subtype[is.na(clin_idh_codel$sample_barcode)] = NA

# For samples where evidence for calls set to NA.
clin_idh_codel$idh_status[clin_idh_codel$sample_barcode=="GLSS-MD-LP03-R1"] = NA
clin_idh_codel$codel_status[clin_idh_codel$sample_barcode=="GLSS-MD-LP03-R1"] = NA
clin_idh_codel$idh_codel_subtype[clin_idh_codel$sample_barcode=="GLSS-MD-LP03-R1"] = NA
clin_idh_codel$idh_status[clin_idh_codel$sample_barcode=="GLSS-MD-LP08-R1"] = NA
clin_idh_codel$codel_status[clin_idh_codel$sample_barcode=="GLSS-MD-LP08-R1"] = NA
clin_idh_codel$idh_codel_subtype[clin_idh_codel$sample_barcode=="GLSS-MD-LP08-R1"] = NA

# Trim all whitespace that may have been introduced with excel formatting.
clin_idh_codel$sample_barcode = trimws(clin_idh_codel$sample_barcode, "r")
clin_idh_codel$grade = trimws(clin_idh_codel$grade, "r")
clin_idh_codel$mgmt_methylation = trimws(clin_idh_codel$mgmt_methylation, "r")
clin_idh_codel$treatment_tmz = trimws(clin_idh_codel$treatment_tmz, "r")
clin_idh_codel$treatment_chemotherapy_other = trimws(clin_idh_codel$treatment_chemotherapy_other, "r")
clin_idh_codel$case_barcode = trimws(clin_idh_codel$case_barcode, "r")


## Check to see whether any samples are not found in the surgeries sheet:
tumor_samples = samples %>% 
  filter(!grepl("-NB",sample_barcode)) %>% 
  filter(!grepl("-NM",sample_barcode))

# Determine potential missing samples:
missing_clinical_indices = which(tumor_samples$sample_barcode%in%clin_idh_codel$sample_barcode==FALSE)
tumor_samples[missing_clinical_indices, ]
# Missing in sample sheet:
missing_sample_indices = which(clin_idh_codel$sample_barcod%in%tumor_samples$sample_barcode==FALSE)
clin_idh_codel[missing_sample_indices, ] %>% 
  filter(!is.na(sample_barcode))

## Write to database.
dbWriteTable(con, Id(schema="clinical",table="surgeries"), clin_idh_codel, append=T)


##############################
# Scripts in this section are only to update clinical.surgeries tabl
##############################
## Update clinical data in the database where there is currently an NA.
# treatment_tmz_cycles, treatment_tmz_cycles_6, treatment_tmz_cycles_notes, treatment_concurrent_tmz

# Floris manually changed this file in the database. Changing here for consistency.
clin_idh_codel$sample_barcode[clin_idh_codel$sample_barcode=="GLSS-MD-0112-R2"] = "GLSS-MD-0112-TP"

### UPDATE: treatment_tmz_cycles
for (i in 1:dim(clin_idh_codel)[1]){
  if(!is.na(clin_idh_codel$treatment_tmz_cycles[i])){
  rs = dbSendStatement(con, sprintf("UPDATE clinical.surgeries SET treatment_tmz_cycles = '%s' WHERE case_barcode = '%s' AND surgery_number= '%s'", clin_idh_codel$treatment_tmz_cycles[i], clin_idh_codel$case_barcode[i], clin_idh_codel$surgery_number[i]))
  dbClearResult(rs)
  } 
}

### UPDATE: treatment_tmz_cycles_6
for (i in 1:dim(clin_idh_codel)[1]){
  if(!is.na(clin_idh_codel$treatment_tmz_cycles_6[i])){
    rs = dbSendStatement(con, sprintf("UPDATE clinical.surgeries SET treatment_tmz_cycles_6 = '%s' WHERE case_barcode = '%s' AND surgery_number= '%s'", clin_idh_codel$treatment_tmz_cycles_6[i], clin_idh_codel$case_barcode[i], clin_idh_codel$surgery_number[i]))
    dbClearResult(rs)
    } 
}

### UPDATE: treatment_tmz_cycles_notes
for (i in 1:dim(clin_idh_codel)[1]){
  if(!is.na(clin_idh_codel$treatment_tmz_cycles_notes[i])){
    rs = dbSendStatement(con, sprintf("UPDATE clinical.surgeries SET treatment_tmz_cycles_notes = '%s' WHERE case_barcode = '%s' AND surgery_number= '%s'", clin_idh_codel$treatment_tmz_cycles_notes[i], clin_idh_codel$case_barcode[i], clin_idh_codel$surgery_number[i]))
    dbClearResult(rs)
    } 
}

### UPDATE: treatment_concurrent_tmz
for (i in 1:dim(clin_idh_codel)[1]){
  if(!is.na(clin_idh_codel$treatment_concurrent_tmz[i])){
   rs = dbSendStatement(con, sprintf("UPDATE clinical.surgeries SET treatment_concurrent_tmz = '%s' WHERE case_barcode = '%s' AND surgery_number= '%s'", clin_idh_codel$treatment_concurrent_tmz[i], clin_idh_codel$case_barcode[i], clin_idh_codel$surgery_number[i]))
   dbClearResult(rs)
  } 
}

## Create new "merge" field to merge case_barcode with surgery number as not all surgeries had sample_barcode.
clin_idh_codel$merge_id = paste0(clin_idh_codel$case_barcode, "_", clin_idh_codel$surgery_number)
surgeries$merge_id = paste0(surgeries$case_barcode, "_", surgeries$surgery_number)
surgeries = surgeries %>% 
  # GLSS-NS-0001
  filter(case_barcode!="GLSS-NS-0001")

# Reorder both tables.
clin_idh_codel = clin_idh_codel[order(clin_idh_codel$merge_id), ]
surgeries = surgeries[order(surgeries$merge_id), ]

# Do all of the IDs match between the surgeries table and clinical.input table?
# Note: this may change as I am directly loading the latest database surgeries table.
clin_idh_codel$merge_id==surgeries$merge_id

# The database recodes "Y' as "1" so set them equal to one another to check whether the new values line up.
clin_idh_codel$treatment_concurrent_tmz[clin_idh_codel$treatment_concurrent_tmz=="Y"] <- "1"
clin_idh_codel$treatment_concurrent_tmz[clin_idh_codel$treatment_concurrent_tmz=="N"] <- "0"
clin_idh_codel$treatment_tmz_cycles_6[clin_idh_codel$treatment_tmz_cycles_6=="Y"] <- "1"
clin_idh_codel$treatment_tmz_cycles_6[clin_idh_codel$treatment_tmz_cycles_6=="N"] <- "0"
# All should be TRUE.
all(clin_idh_codel$treatment_concurrent_tmz==surgeries$treatment_concurrent_tmz, na.rm = T)
all(clin_idh_codel$treatment_tmz_cycles==surgeries$treatment_tmz_cycles, na.rm = T)
all(clin_idh_codel$treatment_tmz_cycles_6==surgeries$treatment_tmz_cycles_6, na.rm = T)
