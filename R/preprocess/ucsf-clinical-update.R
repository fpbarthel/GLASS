#######################################################
# Update newly received clinical data for UCSF samples.
# Date: 2018.11.23
# Author: Kevin J.
#######################################################

# Matthew Grimmer provided additional clinical data on November 16th 2018.
ucsf_clinical_sheet = '/Users/johnsk/Documents/Life-History/ClinicalData/UCSF/2018-1116_glass_wes_clinic_table-costello_Roel.xlsx'

#######################################################
# Necessary packages:
library(tidyverse)
library(openxlsx)
library(DBI)
library(stringr)

#######################################################
# Establish connection with the database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

## Load in clinical data, it may require some processing before use.
ucsf_clinical = readWorkbook(ucsf_clinical_sheet, sheet = 1, startRow = 1, colNames = TRUE)

# Retrieve the case_sources and biospecimen_aliquots from the Database.
cases = dbReadTable(con,  Id(schema="clinical",table="cases"))
surgeries = dbReadTable(con,  Id(schema="clinical",table="surgeries"))

# Subset the cases and surgeries tables to just the patients from UCSF.
ucsf_cases = cases %>% 
  filter(grepl("GLSS-SF-", case_barcode))
 ucsf_surgeries = surgeries %>% 
  filter(grepl("GLSS-SF-", case_barcode))

 # Gather the format of the variables to be uploaded to the cases table.
 str(ucsf_cases)
 
 # Revise these variables to be uploaded to the database.
 ucsf_clinical_db = ucsf_clinical %>% 
   filter(tm_sampletype2 == "TP") %>% 
   filter(!(is.na(age.at.diagnosis))) %>% 
   mutate(patient_number = gsub("patient", "", patientid)) %>% 
   mutate_at("patient_number", str_pad, width = 4, side='left', pad = 0) %>% 
     mutate(case_barcode = paste("GLSS-SF", patient_number, sep="-")) %>% 
   left_join(ucsf_cases, by="case_barcode") %>% 
 mutate(revise_case_vital_status = recode(vital.status, "A" = "alive", "D"="dead"),
        revise_case_overall_survival_mo = round(as.numeric(overall.survival.mo)), 
        revise_case_age_diagnosis_years = floor(age.at.diagnosis)) 
 
# First update the `case_age_diagnosis_years` variable.
for (i in 1:dim(ucsf_clinical_db)[1]){
if(is.na(ucsf_clinical_db$case_age_diagnosis_years[i])){
    rs = dbSendStatement(con, sprintf("UPDATE clinical.cases SET case_age_diagnosis_years = '%s' WHERE case_barcode = '%s'", ucsf_clinical_db$revise_case_age_diagnosis_years[i], ucsf_clinical_db$case_barcode[i]))
    dbClearResult(rs)
  print(ucsf_clinical_db$case_barcode[i])
  print(ucsf_clinical_db$revise_case_age_diagnosis_years[i])
  } 
}

 # Next update the `case_vital_status` variable.
 for (i in 1:dim(ucsf_clinical_db)[1]){
   if(is.na(ucsf_clinical_db$case_vital_status[i])){
     rs = dbSendStatement(con, sprintf("UPDATE clinical.cases SET case_vital_status = '%s' WHERE case_barcode = '%s'", ucsf_clinical_db$revise_case_vital_status[i], ucsf_clinical_db$case_barcode[i]))
     dbClearResult(rs)
     print(ucsf_clinical_db$case_barcode[i])
     print(ucsf_clinical_db$revise_case_vital_status[i])
   } 
 }
 
 # Finally, update the `case_overall_survival_mo` variable.
 for (i in 1:dim(ucsf_clinical_db)[1]){
   if(is.na(ucsf_clinical_db$case_overall_survival_mo[i])){
     rs = dbSendStatement(con, sprintf("UPDATE clinical.cases SET case_overall_survival_mo = '%s' WHERE case_barcode = '%s'", ucsf_clinical_db$revise_case_overall_survival_mo[i], ucsf_clinical_db$case_barcode[i]))
     dbClearResult(rs)
     print(ucsf_clinical_db$case_barcode[i])
     print(ucsf_clinical_db$revise_case_overall_survival_mo[i])
   } 
 }
 # NOTE: GLSS-SF-0081 is still missing the `case_overall_survival_mo` variable.

# It's more difficult to amend the surgeries table because of the clinical variables' format.
# Instead of using a loop, the objective is to manually enter each field.
 ucsf_surgery_db = ucsf_clinical %>% 
   mutate(patient_number = gsub("patient", "", patientid)) %>% 
   mutate_at("patient_number", str_pad, width = 4, side='left', pad = 0) %>% 
   mutate(sample_barcode = paste("GLSS-SF", patient_number, tm_sampletype2, sep="-"))
 
#########################
# Manually update values. Use dbClearResult to prevent error.
#########################
### surgical_interval_mo
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET surgical_interval_mo = '4' WHERE sample_barcode = 'GLSS-SF-0131-R1'")
dbClearResult(rs) 
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET surgical_interval_mo = '26' WHERE sample_barcode = 'GLSS-SF-0157-R1'")
dbClearResult(rs) 
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET surgical_interval_mo = '63' WHERE sample_barcode = 'GLSS-SF-0334-R1'")
dbClearResult(rs) 
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET surgical_interval_mo = '18' WHERE sample_barcode = 'GLSS-SF-0339-R1'")
dbClearResult(rs) 
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET surgical_interval_mo = '43' WHERE sample_barcode = 'GLSS-SF-0060-R1'")
dbClearResult(rs) 
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET surgical_interval_mo = '39' WHERE sample_barcode = 'GLSS-SF-0081-R1'")
dbClearResult(rs) 

### histology 
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET histology = 'Glioblastoma' WHERE sample_barcode = 'GLSS-SF-0131-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET histology = 'Glioblastoma' WHERE sample_barcode = 'GLSS-SF-0131-R1'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET histology = 'Oligoastrocytoma' WHERE sample_barcode = 'GLSS-SF-0157-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET histology = 'Oligoastrocytoma' WHERE sample_barcode = 'GLSS-SF-0157-R1'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET histology = 'Oligodendroglioma' WHERE sample_barcode = 'GLSS-SF-0334-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET histology = 'Oligodendroglioma' WHERE sample_barcode = 'GLSS-SF-0334-R1'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET histology = 'Astrocytoma' WHERE sample_barcode = 'GLSS-SF-0338-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET histology = 'Astrocytoma' WHERE sample_barcode = 'GLSS-SF-0338-R1'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET histology = 'Astrocytoma' WHERE sample_barcode = 'GLSS-SF-0339-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET histology = 'Astrocytoma' WHERE sample_barcode = 'GLSS-SF-0339-R1'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET histology = 'Astrocytoma' WHERE sample_barcode = 'GLSS-SF-0081-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET histology = 'Astrocytoma' WHERE sample_barcode = 'GLSS-SF-0081-R1'")
dbClearResult(rs)

### who_classification
 
### grade
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET grade = 'IV' WHERE sample_barcode = 'GLSS-SF-0131-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET grade = 'IV' WHERE sample_barcode = 'GLSS-SF-0131-R1'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET grade = 'II' WHERE sample_barcode = 'GLSS-SF-0157-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET grade = 'IV' WHERE sample_barcode = 'GLSS-SF-0157-R1'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET grade = 'II' WHERE sample_barcode = 'GLSS-SF-0159-R1'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET grade = 'II' WHERE sample_barcode = 'GLSS-SF-0334-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET grade = 'II' WHERE sample_barcode = 'GLSS-SF-0334-R1'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET grade = 'II' WHERE sample_barcode = 'GLSS-SF-0338-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET grade = 'III' WHERE sample_barcode = 'GLSS-SF-0338-R1'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET grade = 'II' WHERE sample_barcode = 'GLSS-SF-0339-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET grade = 'II' WHERE sample_barcode = 'GLSS-SF-0339-R1'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET grade = 'II' WHERE sample_barcode = 'GLSS-SF-0065-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET grade = 'II' WHERE sample_barcode = 'GLSS-SF-0081-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET grade = 'II' WHERE sample_barcode = 'GLSS-SF-0081-R1'")
dbClearResult(rs)

### treatment_tmz 
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_tmz = '1' WHERE sample_barcode = 'GLSS-SF-0131-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_tmz = '1' WHERE sample_barcode = 'GLSS-SF-0157-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_tmz = '1' WHERE sample_barcode = 'GLSS-SF-0334-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_tmz = '1' WHERE sample_barcode = 'GLSS-SF-0338-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_tmz = '0' WHERE sample_barcode = 'GLSS-SF-0339-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_tmz = '1' WHERE sample_barcode = 'GLSS-SF-0081-TP'")
dbClearResult(rs)


### treatment_chemotherapy_other") treatment_chemotherapy_other_cycles
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_chemotherapy_other = 'Irinotecan, Optune, CBD' WHERE sample_barcode = 'GLSS-SF-0131-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_chemotherapy_other_cycles = '2' WHERE sample_barcode = 'GLSS-SF-0131-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_chemotherapy_other = 'Lomustine' WHERE sample_barcode = 'GLSS-SF-0157-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_chemotherapy_other_cycles = '12' WHERE sample_barcode = 'GLSS-SF-0157-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_chemotherapy_other = 'Lomustine' WHERE sample_barcode = 'GLSS-SF-0170-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_chemotherapy_other_cycles = '12, 10' WHERE sample_barcode = 'GLSS-SF-0170-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_chemotherapy_other = 'Steroids' WHERE sample_barcode = 'GLSS-SF-0032-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_chemotherapy_other_cycles = '4' WHERE sample_barcode = 'GLSS-SF-0032-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_chemotherapy_other = 'Everolimus' WHERE sample_barcode = 'GLSS-SF-0039-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_chemotherapy_other_cycles = '12' WHERE sample_barcode = 'GLSS-SF-0039-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_chemotherapy_other = 'Lomustine' WHERE sample_barcode = 'GLSS-SF-0081-TP'")
dbClearResult(rs)

### treatment_radiotherapy") treatment_radiation_other
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_radiotherapy = '1' WHERE sample_barcode = 'GLSS-SF-0131-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_radiotherapy = '1' WHERE sample_barcode = 'GLSS-SF-0157-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_radiation_other = 'IMRT' WHERE sample_barcode = 'GLSS-SF-0157-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_radiotherapy = '1' WHERE sample_barcode = 'GLSS-SF-0159-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_radiation_other = 'IMRT' WHERE sample_barcode = 'GLSS-SF-0159-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_radiotherapy = '1' WHERE sample_barcode = 'GLSS-SF-0032-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_radiation_other = 'IMRT' WHERE sample_barcode = 'GLSS-SF-0032-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_radiotherapy = '1' WHERE sample_barcode = 'GLSS-SF-0338-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_radiotherapy = '1' WHERE sample_barcode = 'GLSS-SF-0039-TP'")
dbClearResult(rs)
dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_radiation_other = 'IMRT' WHERE sample_barcode = 'GLSS-SF-0039-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_radiotherapy = '1' WHERE sample_barcode = 'GLSS-SF-0060-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_radiation_other = 'Proton beam' WHERE sample_barcode = 'GLSS-SF-0060-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_radiotherapy = '1' WHERE sample_barcode = 'GLSS-SF-0069-TP'")
dbClearResult(rs)
dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_radiation_other = 'IMRT' WHERE sample_barcode = 'GLSS-SF-0069-TP'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE clinical.surgeries SET treatment_radiotherapy = '1' WHERE sample_barcode = 'GLSS-SF-0081-TP'")
dbClearResult(rs)

## Update who_classification
rs =dbSendStatement(con, "UPDATE clinical.surgeries SET who_classification = 'Diffuse Astrocytoma, IDH-mutant' WHERE histology = 'Astrocytoma' AND grade = 'II' AND idh_status = 'IDHmut'")
dbClearResult(rs)
rs =dbSendStatement(con, "UPDATE clinical.surgeries SET who_classification = 'Diffuse Astrocytoma, IDH-wildtype' WHERE histology = 'Astrocytoma' AND grade = 'II' AND idh_status = 'IDHwt'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE  clinical.surgeries SET who_classification = 'Anaplastic Astrocytoma, IDH-mutant' WHERE histology = 'Astrocytoma' AND grade = 'III' AND idh_status = 'IDHmut'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE  clinical.surgeries SET who_classification = 'Anaplastic Astrocytoma, IDH-wildtype' WHERE histology = 'Astrocytoma' AND grade = 'III' AND idh_status = 'IDHwt'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE  clinical.surgeries SET who_classification = 'Glioblastoma, IDH-wildtype' WHERE histology = 'Glioblastoma' AND grade = 'IV' AND idh_status = 'IDHwt'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE  clinical.surgeries SET who_classification = 'Glioblastoma, IDH-mutant' WHERE histology = 'Glioblastoma' AND grade = 'IV' AND idh_status = 'IDHmut'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE  clinical.surgeries SET who_classification = 'Oligodendroglioma, IDH-mutant and 1p/19q-codeleted' WHERE histology = 'Oligodendroglioma' AND grade = 'II' AND idh_status = 'IDHmut' AND codel_status = 'codel'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE  clinical.surgeries SET who_classification = 'Anaplastic Oligodendroglioma, IDH-mutant and 1p/19q-codeleted' WHERE histology = 'Oligodendroglioma' AND grade = 'III' AND idh_status = 'IDHmut' AND codel_status = 'codel'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE  clinical.surgeries SET who_classification = 'Oligodendroglioma, IDH-mutant and 1p/19q-codeleted' WHERE histology IS NOT NULL AND grade = 'II' AND idh_status = 'IDHmut' AND codel_status = 'codel'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE  clinical.surgeries SET who_classification = 'Anaplastic Oligodendroglioma, IDH-mutant and 1p/19q-codeleted' WHERE histology IS NOT NULL AND (grade = 'III' OR grade = 'IV') AND idh_status = 'IDHmut' AND codel_status = 'codel'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE  clinical.surgeries SET who_classification = 'Diffuse Astrocytoma, IDH-mutant' WHERE histology = 'Oligoastrocytoma' AND grade = 'II' AND idh_status = 'IDHmut' AND codel_status = 'noncodel'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE  clinical.surgeries SET who_classification = 'Diffuse Astrocytoma, IDH-wildtype' WHERE histology = 'Oligoastrocytoma' AND grade = 'II' AND idh_status = 'IDHwt' AND codel_status = 'noncodel'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE  clinical.surgeries SET who_classification = 'Anaplastic Astrocytoma, IDH-mutant' WHERE histology = 'Oligoastrocytoma' AND grade = 'III' AND idh_status = 'IDHmut' AND codel_status = 'noncodel'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE  clinical.surgeries SET who_classification = 'Anaplastic Astrocytoma, IDH-wildtype' WHERE histology = 'Oligoastrocytoma' AND grade = 'III' AND idh_status = 'IDHwt' AND codel_status = 'noncodel'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE  clinical.surgeries SET who_classification = 'Diffuse Astrocytoma, IDH-mutant' WHERE histology = 'Oligodendroglioma' AND grade = 'II' AND idh_status = 'IDHmut' AND codel_status = 'noncodel'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE  clinical.surgeries SET who_classification = 'Diffuse Astrocytoma, IDH-wildtype' WHERE histology = 'Oligodendroglioma' AND grade = 'II' AND idh_status = 'IDHwt' AND codel_status = 'noncodel'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE  clinical.surgeries SET who_classification = 'Anaplastic Astrocytoma, IDH-mutant' WHERE histology = 'Oligodendroglioma' AND grade = 'III' AND idh_status = 'IDHmut' AND codel_status = 'noncodel'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE  clinical.surgeries SET who_classification = 'Anaplastic Astrocytoma, IDH-wildtype' WHERE histology = 'Oligodendroglioma' AND grade = 'III' AND idh_status = 'IDHwt' AND codel_status = 'noncodel'")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE  clinical.surgeries SET who_classification = 'Diffuse Astrocytoma, NOS' WHERE histology = 'Astrocytoma' AND grade = 'II' AND idh_status IS NULL AND codel_status IS NULL")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE  clinical.surgeries SET who_classification = 'Anaplastic Astrocytoma, NOS' WHERE histology = 'Astrocytoma' AND grade = 'III' AND idh_status IS NULL AND codel_status IS NULL")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE  clinical.surgeries SET who_classification = 'Oligodendroglioma, NOS' WHERE histology = 'Oligodendroglioma' AND grade = 'II' AND idh_status IS NULL AND codel_status IS NULL")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE  clinical.surgeries SET who_classification = 'Anaplastic Oligodendroglioma, NOS' WHERE histology = 'Oligodendroglioma' AND (grade = 'III' OR grade = 'IV') AND idh_status IS NULL AND codel_status IS NULL")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE  clinical.surgeries SET who_classification = 'Glioblastoma, NOS' WHERE histology = 'Glioblastoma' AND grade = 'IV' AND idh_status IS NULL AND codel_status IS NULL")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE  clinical.surgeries SET who_classification = 'Oligoastrocytoma, NOS' WHERE histology = 'Oligoastrocytoma' AND grade = 'II' AND idh_status IS NULL AND codel_status IS NULL")
dbClearResult(rs)
rs = dbSendStatement(con, "UPDATE  clinical.surgeries SET who_classification = 'Anaplastic Oligoastrocytoma, NOS' WHERE histology = 'Oligoastrocytoma' AND grade = 'III' AND idh_status IS NULL AND codel_status IS NULL")
dbClearResult(rs)

