##############################################
# Analyze results from neutralitytestR applied to each aliquot
# Updated: 2019.01.24
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
require(alluvial)
library(ggExtra)
library(EnvStats)

#######################################################
# Establish connection with Floris' database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

# Load additional tables from the database.
pairs = dbReadTable(con,  Id(schema="analysis", table="pairs"))
tumor_pairs = dbReadTable(con,  Id(schema="analysis", table="tumor_pairs"))
aliquots = dbReadTable(con,  Id(schema="biospecimen", table="aliquots"))
surgeries = dbReadTable(con,  Id(schema="clinical", table="surgeries"))
cases = dbReadTable(con,  Id(schema="clinical", table="cases"))
titan_param = dbReadTable(con,  Id(schema="analysis", table="titan_params"))
mut_freq = dbReadTable(con,  Id(schema="analysis", table="mutation_freq"))
aliquot_neutrality = dbReadTable(con,  Id(schema="analysis", table="neutrality_aliquots"))
clinal_tumor_pairs = dbReadTable(con,  Id(schema="clinical", table="clinical_by_tumor_pair"))  

# These tables **MAY** change, especially the driver table.
silver_set = dbGetQuery(con, "SELECT * FROM analysis.silver_set")
all_drivers = dbGetQuery(con, "SELECT * FROM analysis.driver_status")

# Define set of variables to use:
silver_drivers = silver_set %>% 
  inner_join(all_drivers, by="tumor_pair_barcode") %>% 
  inner_join(clinal_tumor_pairs, by="tumor_pair_barcode")

tmp = silver_drivers %>% 
  group_by(idh_codel_subtype, arm_driver_stability, snv_driver_stability, cnv_driver_stability) %>% 
  summarise(samle_size = n())

## Create revised survival variables for Kaplan-Meier curves and log-rank tests.
# The status indicator, normally 0=alive, 1=dead.
silver_drivers$patient_vital = ifelse(silver_drivers$case_vital_status=="alive", 0, 1)

# Censor all subjects still alive at 36 months (not appropriate, but RV requested it).
silver_drivers$case_survival_cutoff = ifelse(silver_drivers$case_overall_survival_mo>=36, 36, silver_drivers$case_overall_survival_mo)
silver_drivers$patient_vital_cutoff = ifelse(silver_drivers$case_overall_survival_mo>=36, 0, silver_drivers$patient_vital)

# Create a variable for post-recurrence survival.
silver_drivers$patient_post_recur_surv = silver_drivers$case_overall_survival_mo-silver_drivers$surgical_interval
# Treat "neutrality" and "selection" as group by itself.
silver_drivers$binary_mode = ifelse(silver_drivers$evolution_mode=="N-N", "Neutral-Neutral", "Some-Selection")
silver_drivers$s_binary_mode = ifelse(silver_drivers$evolution_mode=="S-S", "Selection-Selection", "Some-Neutrality")
silver_drivers$any_driver_stability = ifelse(silver_drivers$snv_driver_stability=="Driver unstable" | silver_drivers$cnv_driver_stability=="Driver unstable", "Driver unstable", "Driver stable")

# Subset to IDHwt, IDHmut noncodel, and IDHmut codel groups for silver set.
silver_drivers_IDHwt = silver_drivers %>% 
  filter(idh_codel_subtype %in% c("IDHwt_noncodel"))
silver_drivers_IDHmut_noncodel = silver_drivers %>% 
  filter(idh_codel_subtype %in% c("IDHmut_noncodel"))
silver_drivers_IDHmut_codel = silver_drivers %>% 
  filter(idh_codel_subtype %in% c("IDHmut_codel"))

# Perform subtype-specific analyses for associations with different driver changes.
kruskal.test(silver_drivers$surgical_interval, as.factor(silver_drivers$snv_driver_stability))
kruskal.test(silver_drivers_IDHwt$surgical_interval, as.factor(silver_drivers_IDHwt$snv_driver_stability))
kruskal.test(silver_drivers_IDHmut_noncodel$surgical_interval, as.factor(silver_drivers_IDHmut_noncodel$snv_driver_stability))
kruskal.test(silver_drivers_IDHmut_codel$surgical_interval, as.factor(silver_drivers_IDHmut_codel$snv_driver_stability))
kruskal.test(silver_drivers_IDHwt$surgical_interval, as.factor(silver_drivers_IDHwt$cnv_driver_stability))
kruskal.test(silver_drivers_IDHmut_noncodel$surgical_interval, as.factor(silver_drivers_IDHmut_noncodel$cnv_driver_stability))
kruskal.test(silver_drivers_IDHmut_codel$surgical_interval, as.factor(silver_drivers_IDHmut_codel$cnv_driver_stability))
kruskal.test(silver_drivers_IDHwt$surgical_interval, as.factor(silver_drivers_IDHwt$any_driver_stability))
kruskal.test(silver_drivers_IDHmut_noncodel$surgical_interval, as.factor(silver_drivers_IDHmut_noncodel$any_driver_stability))
kruskal.test(silver_drivers_IDHmut_codel$surgical_interval, as.factor(silver_drivers_IDHmut_codel$any_driver_stability))

fisher.test(table(silver_drivers$snv_driver_stability, silver_drivers$idh_codel_subtype))
fisher.test(table(silver_drivers$snv_driver_stability, silver_drivers$recurrence_location))
fisher.test(table(silver_drivers$snv_driver_stability, silver_drivers$received_tmz))
fisher.test(table(silver_drivers$snv_driver_stability, silver_drivers$received_rt))
fisher.test(table(silver_drivers$snv_driver_stability, silver_drivers$hypermutator_status))
fisher.test(table(silver_drivers$snv_driver_stability, silver_drivers$grade_change))
fisher.test(table(silver_drivers$snv_driver_stability, silver_drivers$case_sex))

fisher.test(table(silver_drivers$cnv_driver_stability, silver_drivers$idh_codel_subtype))
fisher.test(table(silver_drivers$cnv_driver_stability, silver_drivers$recurrence_location))
fisher.test(table(silver_drivers$cnv_driver_stability, silver_drivers$received_tmz))
fisher.test(table(silver_drivers$cnv_driver_stability, silver_drivers$received_rt))
fisher.test(table(silver_drivers$cnv_driver_stability, silver_drivers$hypermutator_status))
fisher.test(table(silver_drivers$cnv_driver_stability, silver_drivers$grade_change))
fisher.test(table(silver_drivers$cnv_driver_stability, silver_drivers$case_sex))

### IDHmut codels
fisher.test(table(silver_drivers_IDHmut_codel$snv_driver_stability, silver_drivers_IDHmut_codel$recurrence_location))
fisher.test(table(silver_drivers_IDHmut_codel$snv_driver_stability, silver_drivers_IDHmut_codel$received_tmz))
fisher.test(table(silver_drivers_IDHmut_codel$snv_driver_stability, silver_drivers_IDHmut_codel$received_rt))
fisher.test(table(silver_drivers_IDHmut_codel$snv_driver_stability, silver_drivers_IDHmut_codel$hypermutator_status))
fisher.test(table(silver_drivers_IDHmut_codel$snv_driver_stability, silver_drivers_IDHmut_codel$grade_change))

fisher.test(table(silver_drivers_IDHmut_codel$cnv_driver_stability, silver_drivers_IDHmut_codel$recurrence_location))
fisher.test(table(silver_drivers_IDHmut_codel$cnv_driver_stability, silver_drivers_IDHmut_codel$received_tmz))
fisher.test(table(silver_drivers_IDHmut_codel$cnv_driver_stability, silver_drivers_IDHmut_codel$received_rt))
fisher.test(table(silver_drivers_IDHmut_codel$cnv_driver_stability, silver_drivers_IDHmut_codel$hypermutator_status))
fisher.test(table(silver_drivers_IDHmut_codel$cnv_driver_stability, silver_drivers_IDHmut_codel$grade_change))

fisher.test(table(silver_drivers_IDHmut_codel$any_driver_stability, silver_drivers_IDHmut_codel$recurrence_location))
fisher.test(table(silver_drivers_IDHmut_codel$any_driver_stability, silver_drivers_IDHmut_codel$received_tmz))
fisher.test(table(silver_drivers_IDHmut_codel$any_driver_stability, silver_drivers_IDHmut_codel$received_rt))
fisher.test(table(silver_drivers_IDHmut_codel$any_driver_stability, silver_drivers_IDHmut_codel$hypermutator_status))
fisher.test(table(silver_drivers_IDHmut_codel$any_driver_stability, silver_drivers_IDHmut_codel$grade_change))

### IDHmut noncodel
fisher.test(table(silver_drivers_IDHmut_noncodel$snv_driver_stability, silver_drivers_IDHmut_noncodel$recurrence_location))
# Driver stability associated with TMZ treatment
fisher.test(table(silver_drivers_IDHmut_noncodel$snv_driver_stability, silver_drivers_IDHmut_noncodel$received_tmz))
fisher.test(table(silver_drivers_IDHmut_noncodel$snv_driver_stability, silver_drivers_IDHmut_noncodel$received_rt))
fisher.test(table(silver_drivers_IDHmut_noncodel$snv_driver_stability, silver_drivers_IDHmut_noncodel$hypermutator_status))
fisher.test(table(silver_drivers_IDHmut_noncodel$snv_driver_stability, silver_drivers_IDHmut_noncodel$grade_change))

fisher.test(table(silver_drivers_IDHmut_noncodel$cnv_driver_stability, silver_drivers_IDHmut_noncodel$recurrence_location))
fisher.test(table(silver_drivers_IDHmut_noncodel$cnv_driver_stability, silver_drivers_IDHmut_noncodel$received_tmz))
# Relationship between CNV driver stability and RT.
fisher.test(table(silver_drivers_IDHmut_noncodel$cnv_driver_stability, silver_drivers_IDHmut_noncodel$received_rt))
fisher.test(table(silver_drivers_IDHmut_noncodel$cnv_driver_stability, silver_drivers_IDHmut_noncodel$hypermutator_status))
fisher.test(table(silver_drivers_IDHmut_noncodel$cnv_driver_stability, silver_drivers_IDHmut_noncodel$grade_change))

fisher.test(table(silver_drivers_IDHmut_noncodel$any_driver_stability, silver_drivers_IDHmut_noncodel$recurrence_location))
fisher.test(table(silver_drivers_IDHmut_noncodel$any_driver_stability, silver_drivers_IDHmut_noncodel$received_tmz))
fisher.test(table(silver_drivers_IDHmut_noncodel$any_driver_stability, silver_drivers_IDHmut_noncodel$received_rt))
fisher.test(table(silver_drivers_IDHmut_noncodel$any_driver_stability, silver_drivers_IDHmut_noncodel$hypermutator_status))
# Association between grade change and any driver instability
fisher.test(table(silver_drivers_IDHmut_noncodel$any_driver_stability, silver_drivers_IDHmut_noncodel$grade_change))
# Association between grade change and hypermutation status. Grade change not associated with treatment (P >> 0.05)
fisher.test(table(silver_drivers_IDHmut_noncodel$hypermutator_status, silver_drivers_IDHmut_noncodel$grade_change))

### IDHwt tumors


# Driver stability:
fit_all_stability <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ driver_stability,
                             data = silver_drivers)
ggsurvplot(fit_all_stability, data = silver_drivers, risk.table = TRUE, pval= TRUE)
fit_codel_stability <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ driver_stability,
                               data = silver_drivers_IDHmut_codel)
ggsurvplot(fit_codel_stability, data = silver_drivers_IDHmut_codel, risk.table = TRUE, pval= TRUE)
fit_noncodel_stability <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ driver_stability,
                                  data = silver_drivers_IDHmut_noncodel)
ggsurvplot(fit_noncodel_stability, data = silver_drivers_IDHmut_noncodel, risk.table = TRUE, pval= TRUE)
fit_wt_stability <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ driver_stability,
                            data = silver_drivers_IDHwt)
ggsurvplot(fit_wt_stability, data = silver_drivers_IDHwt, risk.table = TRUE, pval= TRUE)


res_cox_drivers <- coxph(Surv(case_overall_survival_mo, patient_vital)~ case_age_diagnosis_years + idh_codel_subtype + any_driver_stability, data = silver_drivers)
summary(res_cox_drivers)

