##################################################
# Determine association between mut_freq and aneuploidy with survival
# Updated: 2019.05.28
# Author: Kevin J.
##################################################

## OBJECTIVE: Test whether mutational frequency or chromosome instability at 
# initial tumor or recurrent tumor are associated with poor survival. 

# Working directory for this analysis in the GLASS-analysis project. 
mybasedir = "/Volumes/verhaak-lab/GLASS-analysis/"
setwd(mybasedir)

##################################################

# Necessary packages:
library(tidyverse)
library(DBI)
library(survival)

##################################################
# Establish connection with database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB2")

# Essential tables:
gold_set = dbReadTable(con,  Id(schema="analysis", table="gold_set"))
subtypes = dbReadTable(con,  Id(schema="clinical", table="subtypes"))
cases = dbReadTable(con,  Id(schema="clinical", table="cases"))
tumor_mut_comparison = dbReadTable(con,  Id(schema="analysis", table="tumor_mut_comparison"))

# Construct a table that provides subject-level information about clinical variables between two timepoints.
clinical_tumor_pairs_query = read_file("sql/clinical/clinical-tumor-pairs-db2.sql")
clinical_tumor_pairs <- dbGetQuery(con, clinical_tumor_pairs_query)

# Define aneuploidy pairs from the gold set.
aneuploidy_pairs <- dbGetQuery(con, "SELECT tumor_pair_barcode, case_barcode, tumor_barcode_a, tumor_barcode_b, a1.prop_aneuploidy AS aneuploidy_a, a1.aneuploidy_score::integer AS aneuploidy_score_a, 
                               a1.aneuploidy_amp_score::integer AS aneuploidy_amp_score_a, a1.aneuploidy_del_score::integer AS aneuploidy_del_score_a, 
                               a2.prop_aneuploidy AS aneuploidy_b, a2.aneuploidy_score::integer AS aneuploidy_score_b, a2.aneuploidy_amp_score::integer AS aneuploidy_amp_score_b, a2.aneuploidy_del_score::integer AS aneuploidy_del_score_b
                               FROM analysis.gold_set gs
                               LEFT JOIN analysis.gatk_aneuploidy a1 ON a1.aliquot_barcode = gs.tumor_barcode_a
                               LEFT JOIN analysis.gatk_aneuploidy a2 ON a2.aliquot_barcode = gs.tumor_barcode_b")

# Incorporate tables containing age, survival, and subtype.
aneuploidy_pairs_clin = aneuploidy_pairs %>% 
  left_join(clinical_tumor_pairs, by=c("tumor_pair_barcode", "case_barcode", "tumor_barcode_a", "tumor_barcode_b")) %>% 
  filter(tumor_pair_barcode %in% gold_set$tumor_pair_barcode) %>% 
  left_join(subtypes, by="case_barcode") %>% 
  left_join(cases, by = "case_barcode")
aneuploidy_pairs_clin$patient_vital = ifelse(aneuploidy_pairs_clin$case_vital_status=="alive", 0, 1)


# How are the aneuploidy values distributed?
# aneuploidy_value is approx. normally distributed.
ggplot(aneuploidy_pairs_clin, aes(x=aneuploidy_a)) + geom_histogram() + theme_bw()
ggplot(aneuploidy_pairs_clin, aes(x=aneuploidy_b)) + geom_histogram() + theme_bw()
ggplot(aneuploidy_pairs_clin, aes(x=aneuploidy_score_a)) + geom_histogram() + theme_bw()
ggplot(aneuploidy_pairs_clin, aes(x=aneuploidy_score_b)) + geom_histogram() + theme_bw()

# Build a Cox proportion hazards model for each individual aliquot-level variable. General aneuploidy.
aneu_cox_model = coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + idh_codel_subtype + aneuploidy_a, data = aneuploidy_pairs_clin)
summary(aneu_cox_model)
aneu_cox_model = coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + idh_codel_subtype + aneuploidy_b, data = aneuploidy_pairs_clin)
summary(aneu_cox_model)
# Arm-level aneuploidy for the initial tumor.
aneu_cox_model = coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + idh_codel_subtype + aneuploidy_score_a, data = aneuploidy_pairs_clin)
summary(aneu_cox_model)
aneu_cox_model = coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + idh_codel_subtype + aneuploidy_amp_score_a, data = aneuploidy_pairs_clin)
summary(aneu_cox_model)
## The only singificant finding**
aneu_cox_model = coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + idh_codel_subtype + aneuploidy_del_score_a, data = aneuploidy_pairs_clin)
summary(aneu_cox_model)

# Arm-level aneuploidy for the recurrent tumor.
aneu_cox_model = coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + idh_codel_subtype + aneuploidy_score_b, data = aneuploidy_pairs_clin)
summary(aneu_cox_model)
aneu_cox_model = coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + idh_codel_subtype + aneuploidy_amp_score_b, data = aneuploidy_pairs_clin)
summary(aneu_cox_model)
aneu_cox_model = coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + idh_codel_subtype + aneuploidy_del_score_b, data = aneuploidy_pairs_clin)
summary(aneu_cox_model)

# Now, inspect the mutational frequencies.
mf_pairs <- dbGetQuery(con, "SELECT
	tmc.tumor_pair_barcode,
tmc.case_barcode,
tmc.tumor_barcode_a,
tmc.tumor_barcode_b,
idh_codel_subtype,
received_alk,
hypermutator_status,
0 AS time_birth,
ca.case_age_diagnosis_years AS time_initial,
ROUND(ca.case_age_diagnosis_years + (tmc.surgical_interval_mo / 12.0),2) AS time_recurrence,
0 AS mf_birth,
mf1.coverage_adj_mut_freq AS mf_initial,
mf2.coverage_adj_mut_freq AS mf_recurrence,
tmc.count_a,
tmc.count_b,
tmc.union_ab,
tmc.intersection_ab,
tmc.setdiff_a,
tmc.setdiff_b,
mf1.cumulative_coverage AS cov_a,
mf2.cumulative_coverage AS cov_b,
LEAST(mf1.cumulative_coverage, mf2.cumulative_coverage) AS min_cov,
ROUND(setdiff_a::decimal / mf1.cumulative_coverage * 1e6, 4) AS mf_private_a,
ROUND(setdiff_b::decimal / mf2.cumulative_coverage * 1e6, 4) AS mf_private_b,
ROUND(intersection_ab::decimal / LEAST(mf1.cumulative_coverage, mf2.cumulative_coverage) * 1e6, 4) AS mf_shared
FROM analysis.tumor_mut_comparison tmc
INNER JOIN analysis.gold_set stp ON tmc.tumor_pair_barcode = stp.tumor_pair_barcode
LEFT JOIN analysis.tumor_clinical_comparison ctp ON ctp.tumor_pair_barcode = stp.tumor_pair_barcode
LEFT JOIN analysis.mut_freq mf1 ON mf1.aliquot_barcode = tmc.tumor_barcode_a 
LEFT JOIN analysis.mut_freq mf2 ON mf2.aliquot_barcode = tmc.tumor_barcode_b 
LEFT JOIN clinical.subtypes su ON su.case_barcode = stp.case_barcode
LEFT JOIN clinical.cases ca ON ca.case_barcode = stp.case_barcode")

# Inspect the mutational frequency at initial tumor and recurrence.
mf_pairs_annot = mf_pairs %>% 
left_join(clinical_tumor_pairs, by=c("tumor_pair_barcode", "case_barcode", "tumor_barcode_a", "tumor_barcode_b")) %>% 
  filter(tumor_pair_barcode %in% gold_set$tumor_pair_barcode) %>% 
  left_join(cases, by = "case_barcode")
mf_pairs_annot$patient_vital = ifelse(mf_pairs_annot$case_vital_status=="alive", 0, 1)

# mf_a
mf1_cox_model = coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + idh_codel_subtype + log10(mf_initial), data = mf_pairs_annot)
summary(mf1_cox_model)
# mf_b
mf2_cox_model = coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + idh_codel_subtype + log10(mf_recurrence), data = mf_pairs_annot)
summary(mf2_cox_model)

