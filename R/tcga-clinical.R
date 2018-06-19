##################################################################################
# Retrieve the TCGA clinical data for the Primary-Recurrent Tumor Pairs (LGG, GBM)
# Date: 2018.06.19
# Author: Kevin Johnson, Floris Barthel (1st part)
##################################################################################

## Create a clinical manifest of TCGA whole genome (WGS) files from LGG and GBM cohorts
## Limit to primary-recurrent triplets and 2nd recurrences
library(GenomicDataCommons)
library(listviewer)
library(tidyverse)
library(TCGAbiolinks) # To extract information on radiation, drug, follow-up.
library(data.table)
library(dplyr)

# The following code was provided by Floris. This may assume all WGS were cataloged in Legacy.
unpaired_json = "data/ref/TCGA_WGS_GDC_legacy_UUIDs.json"

## Make sure to include case ids
# "cases.project.project_id" = project (eg. TCGA-LGG)
# "cases.samples.sample_type" = sample type (eg. Primary Tumor)
# "cases.samples.portions.analytes.aliquots.submitter_id" = TCGA barcode
# "cases.case_id" = case.id. It turns out regular "id" is also present in the dataset.
# However, it was more reassuring to specify "case_ids".
add_fields = c("cases.project.project_id",
               "cases.samples.sample_type", 
               "cases.samples.portions.analytes.aliquots.submitter_id",
               "cases.case_id")

# Inspect data types.
files(legacy = TRUE) %>% facet(c('data_type')) %>% aggregations()

# Get a list of all WGS aligned BAM files from LGG and GBM tumors.
fq = files(legacy = TRUE) %>% 
  GenomicDataCommons::filter( ~ cases.project.project_id %in% c("TCGA-LGG", "TCGA-GBM") & 
                                experimental_strategy == "WGS" & 
                                data_type == "Aligned reads") %>% 
  GenomicDataCommons::select(c(default_fields(files()), add_fields))

# Provide the total number of whole genomes for LGG and GBM.
message(sprintf("Found %s hits", GenomicDataCommons::count(fq)))

# Extract results.
fres = results_all(fq)

# Inspect list using listviewer::jsonedit.
jsonedit(fres)

# Flatten nested variables.
fres$case_id = unlist(map(fres$cases, "case_id")) # To explicitly grab "case_id" instead of just "id" in case of discrepancies.
fres$project = map(fres$cases, "project") %>% map_chr("project_id")
fres$sample_type = map(fres$cases, "samples") %>% map(unlist) %>% map_chr("sample_type")
fres$aliquot_id = map(fres$cases, "samples") %>% map(unlist) %>% map_chr("portions.analytes.aliquots.submitter_id")

# Convert to dataframe.
df = as.data.frame(fres[-which(names(fres) %in% c("cases", "acl", "analysis"))], stringsAsFactors = F) %>%
  select(id, aliquot_id, project, sample_type, case_id, experimental_strategy, file_size, md5sum, file_name, created_datetime, updated_datetime) %>%
  filter(grepl("TCGA", project)) %>%
  mutate(sample_id = substr(aliquot_id,1,16),
         ali_id = substr(aliquot_id,1,12))

## Pair primary-recurrent-2ndrecurrence samples.
filtered_files = df %>% 
  group_by(sample_id) %>%
  mutate(p = order(file_size, decreasing = T)) %>%
  ungroup() %>%
  group_by(ali_id) %>%
  mutate(hasRec = any(sample_type == "Recurrent Tumor")) %>%
  ungroup() %>%
  filter(hasRec, p == 1) %>%
  select(-hasRec, -p, -created_datetime, -updated_datetime, -experimental_strategy)

# Create a simplified sheet for labelling the sample status.
paired_files = filtered_files %>%
  mutate(sample_type_numeric = recode_factor(substr(aliquot_id, 14,16), "01A" = "P", "01B" = "P", "02A" = "R1", "02B" = "R2", "10A" = "N", "10B" = "N", "10D" = "N")) %>%
  select(ali_id, project, sample_type_numeric, id) %>%
  spread(sample_type_numeric, id)

# Output files in JSON format. Some clinical data present.
write(jsonlite::toJSON(filtered_files, pretty = T), file = unpaired_json)



## Clinical Data (Treatment, Diagnosis)
# We need to recover additional data about these particular samples.
# Get case_ids for each sample to feed into gdc_clinical() from GenomicDataCommons package.
qcases = cases() %>% GenomicDataCommons::select(available_fields('cases'))
# files.cases.diagnoses.treatment
# "files.analysis.metadata.read_groups.read_length"
# "diagnoses.treatments"
clinResults = cases() %>%
  GenomicDataCommons::select(NULL) %>%
  #  GenomicDataCommons::expand(expands) %>%
  results(size=10)








# gdc_clinical grabs "diagnoses", "demographic", and "exposures".
clin_res = gdc_clinical(filtered_files$case_id)

# This generates four separate files: demographic, Dx, exposures, and main. Subject-level.
sapply(clin_res, dim) %>%
  t() %>%
  data.frame() %>%
  set_names(c('rows','columns'))

# Print the variables that are measured in these clincal sheets.
sapply(clin_res, colnames)

# Use lef_join to merge all of the covariate files to the "main" sample sheet.
basic_clin = with(clin_res,
                 main %>%
                   left_join(demographic, by = "case_id") %>%
                   left_join(exposures, by = "case_id") %>%
                   left_join(diagnoses, by = "case_id"))


###########################
# There is STILL some missing clinical information with this tool. We are also interested in 
# drugs, radiation, and follow-up.
############################
# Use TCGABioLinks to retrieve treatment information. 
# Note that the "Clinical" sheet may not have updated information about
# follow-up. Must specifcy clinical.info = "follow_up" for latest update.
LGG_clin = basic_clin %>% filter(disease_type=="Brain Lower Grade Glioma")
LGG_IDs = LGG_clin$submitter_id.x

# Using this tool I had to hardcode the specific sample names. 
# I ran into issues when I fed barcode = LGG_IDs. The software also failed
# when I tried combining the LGG and GBM queries. 
print(LGG_IDs)
query_LGG <- GDCquery(project = "TCGA-LGG", 
                  data.category = "Clinical", 
                  barcode = c("TCGA-DU-6404","TCGA-DU-6397","TCGA-DU-5870",
                              "TCGA-TQ-A7RK", "TCGA-DH-A669", "TCGA-FG-A4MT", "TCGA-FG-5965",
                              "TCGA-TM-A7CF", "TCGA-TQ-A8XE", "TCGA-DU-5872", "TCGA-DU-6407", 
                              "TCGA-TQ-A7RV", "TCGA-DU-7304"))

# The files are downloaded into your home directory. Here it's /Users/johnsk/
GDCdownload(query_LGG)
# Each clinical sheet requires its own preparation.
LGG_clinical = GDCprepare_clinic(query_LGG, clinical.info = "patient")
LGG_clinical_drug = GDCprepare_clinic(query_LGG, clinical.info = "drug")
# Rename to avoid conflict between drug and radiation.
LGG_clinical_drug = LGG_clinical_drug %>% rename(regimen_indication_drug = regimen_indication, regimen_indication_drug_notes = regimen_indication_notes, measure_of_response_drug = measure_of_response) 
LGG_clinical_followup = GDCprepare_clinic(query_LGG, clinical.info = "follow_up")
LGG_clinical_radiation = GDCprepare_clinic(query_LGG, clinical.info = "radiation")
LGG_clinical_radiation = LGG_clinical_radiation %>% rename(regimen_indication_radiation = regimen_indication, 
           regimen_indication_radiation_notes = regimen_indication_notes,
            measure_of_response_radiation = measure_of_response)
# Avoid confusion by removing repeated columns generated by left_join.
# Combine all LGG files into a single sheet.
full_LGG_clinical = LGG_clinical %>% select(-vital_status, -days_to_death, -days_to_last_followup,
                                            -person_neoplasm_cancer_status, -day_of_form_completion, -month_of_form_completion,
                                            -year_of_form_completion, -performance_status_scale_timing, -radiation_therapy,
                                            -karnofsky_performance_score, -eastern_cancer_oncology_group, -days_to_performance_status_assessment,
                                            -targeted_molecular_therapy, -primary_therapy_outcome_success,
                                            -day_of_form_completion, -month_of_form_completion) %>% 
                   left_join(LGG_clinical_drug, by = "bcr_patient_barcode") %>% 
                   left_join(LGG_clinical_followup, by = "bcr_patient_barcode") %>%
                   left_join(LGG_clinical_radiation, by = "bcr_patient_barcode") 
# Double check to see whether any variable was duplicated. Remove for clarity.
colnames(full_LGG_clinical)[grepl("*\\.x", colnames(full_LGG_clinical))]
colnames(full_LGG_clinical)[grepl("*\\.y", colnames(full_LGG_clinical))]

# Spot check to make sure essential variables remain the same.
full_LGG_clinical$day_of_form_completion.x==full_LGG_clinical$day_of_form_completion.y
LGG_clinical_combined = full_LGG_clinical %>% select(-day_of_form_completion.y, -month_of_form_completion.y, -year_of_form_completion.y)

####
# Perform same analyses on GBM (n = 10).
####
GBM_clin = full_clin %>% filter(disease_type=="Glioblastoma Multiforme")
GBM_IDs = GBM_clin$submitter_id.x
query_GBM <- GDCquery(project = "TCGA-GBM", 
                      data.category = "Clinical", 
                      barcode = c("TCGA-06-0221", "TCGA-06-0210", "TCGA-19-1389", "TCGA-14-1402",
                                  "TCGA-14-1034", "TCGA-06-0190", "TCGA-06-0152",
                                  "TCGA-06-0125", "TCGA-06-0211", "TCGA-06-0171"))
GDCdownload(query_GBM)
GBM_clinical <- GDCprepare_clinic(query_GBM, clinical.info = "patient")
GBM_clinical_drug <- GDCprepare_clinic(query_GBM, clinical.info = "drug")
# Rename variables that will be used in other data files
GBM_clinical_drug = GBM.clinical.drug %>% rename(regimen_indication_drug = regimen_indication, regimen_indication_drug_notes = regimen_indication_notes, measure_of_response_drug = measure_of_response) 
GBM_clinical_followup <- GDCprepare_clinic(query_GBM, clinical.info = "follow_up")
GBM_clinical_radiation <- GDCprepare_clinic(query_GBM, clinical.info = "radiation")
# Rename variables not unique to radiation clinical data to improve clarity.
GBM_clinical_radiation = GBM_clinical_radiation %>% rename(regimen_indication_radiation = regimen_indication, 
                                                           regimen_indication_radiation_notes = regimen_indication_notes,
                                                           measure_of_response_radiation = measure_of_response)
# Generate a final cleaned version of the GBM covariates.
full_GBM_clinical = GBM_clinical %>% select(-vital_status, -days_to_death, -days_to_last_followup,
       -person_neoplasm_cancer_status, -day_of_form_completion, -month_of_form_completion,
       -year_of_form_completion, -radiation_therapy,-postoperative_rx_tx, -primary_therapy_outcome_success,
       -karnofsky_performance_score, -eastern_cancer_oncology_group, -performance_status_scale_timing,
       -day_of_form_completion, -month_of_form_completion) %>% 
    left_join(GBM_clinical_drug, by = "bcr_patient_barcode") %>% 
   left_join(GBM_clinical_followup, by = "bcr_patient_barcode") %>%
   left_join(GBM_clinical_radiation, by = "bcr_patient_barcode") 

# Determine whether any variables were present in multiple sheets.
colnames(full_GBM_clinical)[grepl("*\\.x", colnames(full_GBM_clinical))]
colnames(full_GBM_clinical)[grepl("*\\.y", colnames(full_GBM_clinical))]

# Spot check to make sure essential variables remain the same.
full_GBM_clinical$day_of_form_completion.x==full_GBM_clinical$day_of_form_completion.y # This variable actually happens to change in multiple datasets. The rest remain the same.
# Discard repeat variables.
GBM_clinical_combined = full_GBM_clinical %>% select(-day_of_form_completion.y, -month_of_form_completion.y, -year_of_form_completion.y)

#####
# Determine which variables are measured in LGG, but not GBM. What is shared between two sets.
#####
dplyr::setdiff(colnames(LGG_clinical_combined), (colnames(GBM_clinical_combined)))
dplyr::intersect(colnames(LGG_clinical_combined), (colnames(GBM_clinical_combined)))

# Combine both the LGG and GBM clinical data.
TCGA_LGG_GBM_clinical_covars = merge(LGG_clinical_combined, GBM_clinical_combined, all=TRUE)


