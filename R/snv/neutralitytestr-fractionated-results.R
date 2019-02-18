##############################################
# Analyze results from neutralitytestR applied to fractionated mutations
# Updated: 2019.01.26
# Author: Kevin J.
##################################################

# Working directory for this analysis in the GLASS-analysis project. 
mybasedir = "/Volumes/verhaak-lab/GLASS-analysis/"
setwd(mybasedir)

#######################################################
# Necessary packages:
library(tidyverse)
library(DBI)
library(neutralitytestr)
library(survminer)
library(survival)
library(ggExtra)
library(EnvStats)

#######################################################

# Establish connection with Floris' database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

# Load additional tables from the database.
pairs = dbReadTable(con,  Id(schema="analysis",table="pairs"))
tumor_pairs = dbReadTable(con,  Id(schema="analysis", table="tumor_pairs"))
aliquots = dbReadTable(con,  Id(schema="biospecimen", table="aliquots"))
surgeries = dbReadTable(con,  Id(schema="clinical", table="surgeries"))
cases = dbReadTable(con,  Id(schema="clinical", table="cases"))
titan_param = dbReadTable(con,  Id(schema="analysis",table="titan_params"))
mut_freq = dbReadTable(con,  Id(schema="analysis",table="mutation_freq"))
neutrality_aliquots = dbReadTable(con,  Id(schema="analysis",table="neutrality_aliquots"))
# These were the fractionated results calculated with primary, shared, and recurrence.
neutrality_tumor_pairs = dbReadTable(con,  Id(schema="analysis",table="neutrality_tumor_pairs"))

# These tables **MAY** change, especially the driver table.
clinal_tumor_pairs = dbGetQuery(con,"SELECT * FROM analysis.clinical_by_tumor_pair")
drivers = dbGetQuery(con, "SELECT * FROM analysis.driver_status")
silver_set = dbGetQuery(con, "SELECT * FROM analysis.silver_set")

# Combine with pairs file to access "tumor_barcode".
titan_info = titan_param %>% 
  left_join(pairs, by="pair_barcode")

################################################
# **Fractionated** neutrality results
################################################
# Prepare fractionated neutrality results and merge with SILVER set.
fract_neut_annot = neutrality_tumor_pairs %>% 
  left_join(clinal_tumor_pairs, by="tumor_pair_barcode") %>% 
  inner_join(silver_set, by="tumor_pair_barcode") %>% 
  left_join(cases, by=c("case_barcode.x"="case_barcode")) %>% 
  left_join(drivers, by="tumor_pair_barcode")

# Create the appropriate vital status for these events.
fract_neut_annot$patient_vital = ifelse(fract_neut_annot$case_vital_status=="alive", 0, 1)

# There is a good amount of missingness for neutralitytestR in the shared variants.
# data.frame(134 (P), 80 (S_a), 69 (S_b), 106 (R)).
sum(!is.na(fract_neut_annot$primary_evolution))
sum(!is.na(fract_neut_annot$shared_a_evolution))
sum(!is.na(fract_neut_annot$shared_b_evolution))
sum(!is.na(fract_neut_annot$recurrence_evolution))

# Retrieve summary statistics for median surgical interval:
fract_neut_annot %>% 
  group_by(recurrence_evolution) %>% 
  summarise(med_int = median(surgical_interval),
            sample_size = n())

# Perform statistical tests for all associated **FRACTIONATED** evolutionary patterns.
# No variables were significantly associated with the mode of evolution for the recurrent fraction.
wilcox.test(fract_neut_annot$surgical_interval~fract_neut_annot$recurrence_evolution)
fisher.test(table(fract_neut_annot$recurrence_location,fract_neut_annot$recurrence_evolution))
fisher.test(table(fract_neut_annot$driver_stability, fract_neut_annot$recurrence_evolution))
fisher.test(table(fract_neut_annot$received_rt,fract_neut_annot$recurrence_evolution))
fisher.test(table(fract_neut_annot$received_tmz,fract_neut_annot$recurrence_evolution))
fisher.test(table(fract_neut_annot$grade_change,fract_neut_annot$recurrence_evolution))
fisher.test(table(fract_neut_annot$hypermutator_status, fract_neut_annot$recurrence_evolution))
fisher.test(table(fract_neut_annot$case_sex, fract_neut_annot$recurrence_evolution))
fisher.test(table(fract_neut_annot$idh_codel_subtype, fract_neut_annot$recurrence_evolution))

## Survival curve for silver set pairs with evolutionary status of the mutations PRIVATE to the recurrence.
frac_recur <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ recurrence_evolution,
                      data = fract_neut_annot)
ggsurvplot(frac_recur, data = fract_neut_annot, risk.table = TRUE, pval= TRUE, pval.method = TRUE)

# Similar analysis with the primary data. Note: Shared evolution ONLY had ONE neutrality event (mostly selected). Did not conduct test.
primary_recur <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ primary_evolution,
                         data = fract_neut_annot)
ggsurvplot(primary_recur, data = fract_neut_annot, risk.table = TRUE, pval= TRUE, pval.method = TRUE)

## Surgical interval for recurrence-evolution. Log-rank test vs. kruskal wallis or wilcoxon.
frac_int_recur <- survfit(Surv(surgical_interval, rep(1, length(fract_neut_annot$surgical_interval))) ~ recurrence_evolution,
                          data = fract_neut_annot)
ggsurvplot(frac_int_recur, data = fract_neut_annot, risk.table = TRUE, pval= TRUE, ylab = "Proportion Progression Free")

frac_int_primary <- survfit(Surv(surgical_interval, rep(1, length(fract_neut_annot$surgical_interval))) ~ primary_evolution,
                            data = fract_neut_annot)
ggsurvplot(frac_int_primary, data = fract_neut_annot, risk.table = TRUE, pval= TRUE, ylab = "Proportion Progression Free")




