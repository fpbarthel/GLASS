##############################################
# NeutralitytestR applied to each aliquot
# Updated: 2019.01.14
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
neutrality_tumor_pairs = dbReadTable(con,  Id(schema="analysis",table="neutrality_tumor_pairs"))
aneuploidy = dbReadTable(con,  Id(schema="analysis", table="aneuploidy"))



# These tables **MAY** change, especially the driver table.
clinal_tumor_pairs = dbGetQuery(con,"SELECT * FROM analysis.clinical_by_tumor_pair")
drivers = dbGetQuery(con, "SELECT * FROM analysis.driver_status")
silver_set = dbGetQuery(con, "SELECT * FROM analysis.silver_set")

# Combine with pairs file to access "tumor_barcode".
titan_info = titan_param %>% 
  left_join(pairs, by="pair_barcode")

#######################################################
# Traditional neutralitytestR analysis of all variants in each tumor
#######################################################
# To mirror what Georgette is doing with
# the SubclonalSelection script, we will include 
# all tumors
neutrality_input_aliquot_mutect2 = read_file("sql/neutralitytestr-input-aliquot-level.sql")
glass_single_vaf <- dbGetQuery(con, neutrality_input_aliquot_mutect2)

# 514 distinct aliquots.
n_distinct(glass_single_vaf$aliquot_barcode)

# Determine the number of mutations that are considered subclonal in each tumor.
aliquot_mutation_counts = glass_single_vaf %>% 
  filter(vaf >= 0.1 & vaf <= 0.25) %>% 
  group_by(aliquot_barcode) %>% 
  summarize(subclonal_mut = n())

# Store results from per-level neutralitytestR.
aliquot_results = matrix(NA, nrow = length(unique(glass_single_vaf$tumor_pair_barcode)), ncol = 12)
colnames(aliquot_results) =  c("aliquot_barcode", "ploidy", "cellularity", "mutation_rate", "model_rsq","model_pval", "area_value", "area_pval","meanDist_value", "meanDist_pval", "Dk_value", "Dk_pval")
aliquot_results = as.data.frame(aliquot_results)

# Set directory for output.
setwd("/Volumes/verhaak-lab/GLASS-analysis/results/neutral_evolution/neutral-evolution-mutect2/aliquot-level/")

## Generate results for "tumor_barcode_b" using "vaf_b".
for ( i in 1:length(unique(glass_single_vaf$aliquot_barcode))){ 
  # Create subsetted data for this case.
  Yi = glass_single_vaf[glass_single_vaf$aliquot_barcode==unique(glass_single_vaf$aliquot_barcode)[i], ]
  xi = unique(Yi$aliquot_barcode)
  
  print(sprintf("Analyzing all mutations for: %s",  xi))
  # Store sample information.
  aliquot_results[i, "aliquot_barcode"] = xi
  aliquot_results[i, "ploidy"] = titan_info[titan_info$tumor_barcode==xi, "ploidy"]
  aliquot_results[i, "cellularity"] = titan_info[titan_info$tumor_barcode==xi, "purity"]
  
  # The model can sometimes fail. Skip over samples that error out.
  fit_a =  try(neutralitytest(Yi$vaf, ploidy = titan_info[titan_info$tumor_barcode==xi, "ploidy"], cellularity = titan_info[titan_info$tumor_barcode==xi, "purity"]), silent = T) 
  if(!inherits(fit_a, "try-error")) {
    aliquot_results[i, "mutation_rate"] = fit_a$mutation.rate
    aliquot_results[i, "model_rsq"] = fit_a$model$rsq
    aliquot_results[i, "model_pval"] = fit_a$model$pval
    aliquot_results[i, "area_value"] = fit_a$area$metric
    aliquot_results[i, "area_pval"] = fit_a$area$pval
    aliquot_results[i, "meanDist_value"] = fit_a$meanDist$metric
    aliquot_results[i, "meanDist_pval"] = fit_a$meanDist$pval
    aliquot_results[i, "Dk_value"] = fit_a$Dk$metric
    aliquot_results[i, "Dk_pval"] = fit_a$Dk$pval
    
    # Print out model plots to directory of interest.
    #pdf(paste(unique(glass_single_vaf$aliquot_barcode)[i],"-_all",".pdf", sep=""), width=9, height=6)
    print(plot_all(fit_a))
    dev.off()
  }
}

# Combine data and examine the patterns.
aliquot_neutrality = aliquot_results %>% 
  left_join(aliquot_mutation_counts, by="aliquot_barcode") 

# Export to database:
# dbWriteTable(con, Id(schema="analysis",table="neutrality_aliquots"), as.data.frame(aliquot_neutrality))

################################################
# Analyses with fractionated neutrality results
################################################
# Prepare fractionated neutrality results and merge with silver set.
fract_neut_annot = neutrality_tumor_pairs %>% 
  left_join(clinal_tumor_pairs, by="tumor_pair_barcode") %>% 
  inner_join(silver_set, by="tumor_pair_barcode") %>% 
  left_join(cases, by=c("case_barcode.x"="case_barcode")) %>% 
  left_join(drivers, by="tumor_pair_barcode")

# Create the appropriate vital status for these events.
fract_neut_annot$patient_vital = ifelse(fract_neut_annot$case_vital_status=="alive", 0, 1)

# Retrieve summary statistics:
fract_neut_annot %>% 
  group_by(recurrence_evolution) %>% 
  summarise(med_int = median(surgical_interval),
            sample_size = n())

# Perform statistical tests for all associated evolutionary patterns.
wilcox.test(fract_neut_annot$surgical_interval~fract_neut_annot$recurrence_evolution)
fisher.test(table(fract_neut_annot$recurrence_location,fract_neut_annot$recurrence_evolution))
fisher.test(table(fract_neut_annot$driver_stability, fract_neut_annot$recurrence_evolution))
fisher.test(table(fract_neut_annot$received_rt,fract_neut_annot$recurrence_evolution))
fisher.test(table(fract_neut_annot$received_tmz,fract_neut_annot$recurrence_evolution))
fisher.test(table(fract_neut_annot$grade_change,fract_neut_annot$recurrence_evolution))
fisher.test(table(fract_neut_annot$hypermutator_status, fract_neut_annot$recurrence_evolution))
fisher.test(table(fract_neut_annot$case_sex, fract_neut_annot$recurrence_evolution))
fisher.test(table(fract_neut_annot$idh_codel_subtype, fract_neut_annot$recurrence_evolution))

## Survival curve for silver set pairs with evolutionary status of the mutations private to the recurrence.
frac_recur <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ recurrence_evolution,
                         data = fract_neut_annot)
ggsurvplot(frac_recur, data = fract_neut_annot, risk.table = TRUE, pval= TRUE)
# Similar analysis with the primary data. Note: Shared evolution only had one neutrality event.
primary_recur <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ primary_evolution,
                      data = fract_neut_annot)
ggsurvplot(primary_recur, data = fract_neut_annot, risk.table = TRUE, pval= TRUE)

## Surgical interval for recurrence-evolution.
frac_int_primary <- survfit(Surv(surgical_interval, rep(1, length(fract_neut_annot$surgical_interval))) ~ primary_evolution,
                            data = fract_neut_annot)
ggsurvplot(frac_int_primary, data = fract_neut_annot, risk.table = TRUE, pval= TRUE, ylab = "Proportion Progression Free")
frac_int_recur <- survfit(Surv(surgical_interval, rep(1, length(fract_neut_annot$surgical_interval))) ~ recurrence_evolution,
                      data = fract_neut_annot)
ggsurvplot(frac_int_recur, data = fract_neut_annot, risk.table = TRUE, pval= TRUE, ylab = "Proportion Progression Free")


# Although all aliquot results are present here, BUT need to filter low subclonal mutation count AND ploidy > 3.
# Build tumor pairs from individual aliquots. 
case_level_subtype = surgeries %>% 
  select(case_barcode, idh_codel_subtype) %>% 
  distinct() %>% 
  filter(!is.na(idh_codel_subtype)) 

# Combine the neutrality at the aliquot with specific tumor_pairs.
aliquot_neutrality = neutrality_aliquots
per_sample_neutrality = tumor_pairs %>% 
  left_join(aliquot_neutrality, by=c("tumor_barcode_a"="aliquot_barcode")) %>% 
  left_join(aliquot_neutrality, by=c("tumor_barcode_b"="aliquot_barcode")) %>% 
  filter(subclonal_mut.x >= 12, subclonal_mut.y >= 12) %>% 
  mutate(primary_evolution = ifelse(model_rsq.x > 0.98, "N", "S"), 
         recurrence_evolution = ifelse(model_rsq.y > 0.98, "N", "S"),
         evolution_mode = paste(primary_evolution, recurrence_evolution, sep="-")) %>% 
  left_join(case_level_subtype, by="case_barcode") %>% 
  left_join(drivers, by = "tumor_pair_barcode") %>% 
  left_join(clinal_tumor_pairs) %>% 
  mutate(driver_status_revised = recode(driver_status, "Driver gain" = "Driver unstable",
                                        "Driver loss" = "Driver unstable", "Driver switch" = "Driver unstable", "Driver null" = "Driver stable"))
per_sample_neutrality$tmz_std_care <- ifelse(per_sample_neutrality$received_tmz_sum_cycles >= 6, "1", "0")

# Subset to silver set so that a case is only represented at a single time point.
silver_neutrality = per_sample_neutrality %>% 
  inner_join(silver_set, by="tumor_pair_barcode") %>% 
  left_join(cases, by=c("case_barcode.x"="case_barcode"))

# Evolution is marginally associated with idh_subtype.
fisher.test(table(silver_neutrality$idh_codel_subtype, silver_neutrality$recurrence_evolution))
fisher.test(table(silver_neutrality$idh_codel_subtype, silver_neutrality$evolution_mode))

## Survival
# The status indicator, normally 0=alive, 1=dead.
silver_neutrality$patient_vital = ifelse(silver_neutrality$case_vital_status=="alive", 0, 1)

# Censor all subjects still alive at 36 months.
ifelse(silver_neutrality$case_overall_survival_mo>=36, 36, silver_neutrality$case_overall_survival_mo)




# Subset to IDHwt,
silver_neutrality_IDHwt = silver_neutrality %>% 
  filter(idh_codel_subtype %in% c("IDHwt_noncodel"))
silver_neutrality_IDHmut_noncodel = silver_neutrality %>% 
  filter(idh_codel_subtype %in% c("IDHmut_noncodel"))
silver_neutrality_IDHmut_codel = silver_neutrality %>% 
  filter(idh_codel_subtype %in% c("IDHmut_codel"))

# All
silver_neutrality %>% 
  group_by(recurrence_evolution, idh_codel_subtype) %>% 
  summarise(med_int = median(surgical_interval, na.rm = T),
            sample_size = n())
wilcox.test(silver_neutrality$surgical_interval~silver_neutrality$recurrence_evolution)
fisher.test(table(silver_neutrality$recurrence_evolution, silver_neutrality$idh_codel_subtype))
fisher.test(table(silver_neutrality$recurrence_evolution, silver_neutrality$driver_status))
fisher.test(table(silver_neutrality$recurrence_evolution, silver_neutrality$driver_stability))
fisher.test(table(silver_neutrality$recurrence_evolution, silver_neutrality$recurrence_location))
fisher.test(table(silver_neutrality$recurrence_evolution, silver_neutrality$case_sex))
fisher.test(table(silver_neutrality$recurrence_evolution, silver_neutrality$received_tmz))
fisher.test(table(silver_neutrality$recurrence_evolution, silver_neutrality$tmz_std_care))
fisher.test(table(silver_neutrality$recurrence_evolution, silver_neutrality$received_rt))
fisher.test(table(silver_neutrality$recurrence_evolution, silver_neutrality$hypermutator_status))
fisher.test(table(silver_neutrality$recurrence_evolution, silver_neutrality$grade_change))

# IDHmut codel
wilcox.test(silver_neutrality_IDHmut_codel$surgical_interval~silver_neutrality_IDHmut_codel$recurrence_evolution)
fisher.test(table(silver_neutrality_IDHmut_codel$recurrence_evolution, silver_neutrality_IDHmut_codel$driver_status))
fisher.test(table(silver_neutrality_IDHmut_codel$recurrence_evolution, silver_neutrality_IDHmut_codel$driver_stability))
fisher.test(table(silver_neutrality_IDHmut_codel$recurrence_evolution, silver_neutrality_IDHmut_codel$recurrence_location))
fisher.test(table(silver_neutrality_IDHmut_codel$recurrence_evolution, silver_neutrality_IDHmut_codel$received_tmz))
fisher.test(table(silver_neutrality_IDHmut_codel$recurrence_evolution, silver_neutrality_IDHmut_codel$tmz_std_care))
fisher.test(table(silver_neutrality_IDHmut_codel$recurrence_evolution, silver_neutrality_IDHmut_codel$received_rt))
fisher.test(table(silver_neutrality_IDHmut_codel$recurrence_evolution, silver_neutrality_IDHmut_codel$hypermutator_status))
fisher.test(table(silver_neutrality_IDHmut_codel$recurrence_evolution, silver_neutrality_IDHmut_codel$grade_change))
# IDHmut noncodel
wilcox.test(silver_neutrality_IDHmut_noncodel$surgical_interval~silver_neutrality_IDHmut_noncodel$recurrence_evolution)
fisher.test(table(silver_neutrality_IDHmut_noncodel$recurrence_evolution, silver_neutrality_IDHmut_noncodel$driver_status))
fisher.test(table(silver_neutrality_IDHmut_noncodel$recurrence_evolution, silver_neutrality_IDHmut_noncodel$driver_stability))
fisher.test(table(silver_neutrality_IDHmut_noncodel$recurrence_evolution, silver_neutrality_IDHmut_noncodel$recurrence_location))
fisher.test(table(silver_neutrality_IDHmut_noncodel$recurrence_evolution, silver_neutrality_IDHmut_noncodel$received_tmz))
fisher.test(table(silver_neutrality_IDHmut_noncodel$recurrence_evolution, silver_neutrality_IDHmut_noncodel$tmz_std_care))
fisher.test(table(silver_neutrality_IDHmut_noncodel$recurrence_evolution, silver_neutrality_IDHmut_noncodel$received_rt))
fisher.test(table(silver_neutrality_IDHmut_noncodel$recurrence_evolution, silver_neutrality_IDHmut_noncodel$hypermutator_status))
fisher.test(table(silver_neutrality_IDHmut_noncodel$recurrence_evolution, silver_neutrality_IDHmut_noncodel$grade_change))

# IDHwt
wilcox.test(silver_neutrality_IDHwt$surgical_interval~silver_neutrality_IDHwt$recurrence_evolution)
fisher.test(table(silver_neutrality_IDHwt$recurrence_evolution, silver_neutrality_IDHwt$driver_status))
fisher.test(table(silver_neutrality_IDHwt$recurrence_evolution, silver_neutrality_IDHwt$driver_stability))
fisher.test(table(silver_neutrality_IDHwt$recurrence_evolution, silver_neutrality_IDHwt$recurrence_location))
fisher.test(table(silver_neutrality_IDHwt$recurrence_evolution, silver_neutrality_IDHwt$received_tmz))
fisher.test(table(silver_neutrality_IDHwt$recurrence_evolution, silver_neutrality_IDHwt$tmz_std_care))
fisher.test(table(silver_neutrality_IDHwt$recurrence_evolution, silver_neutrality_IDHwt$received_rt))
fisher.test(table(silver_neutrality_IDHwt$recurrence_evolution, silver_neutrality_IDHwt$hypermutator_status))
fisher.test(table(silver_neutrality_IDHwt$recurrence_evolution, silver_neutrality_IDHwt$grade_change))

## Evolution Mode
# All
silver_neutrality %>% 
  group_by(idh_codel_subtype, evolution_mode) %>% 
  summarise(med_int = median(surgical_interval, na.rm = T),
            sample_size = n())
kruskal.test(silver_neutrality$surgical_interval, as.factor(silver_neutrality$evolution_mode))
kruskal.test(silver_neutrality_IDHmut_codel$surgical_interval, as.factor(silver_neutrality_IDHmut_codel$evolution_mode))
kruskal.test(silver_neutrality_IDHmut_noncodel$surgical_interval, as.factor(silver_neutrality_IDHmut_noncodel$evolution_mode))
kruskal.test(silver_neutrality_IDHwt$surgical_interval, as.factor(silver_neutrality_IDHwt$evolution_mode))

fisher.test(table(silver_neutrality$evolution_mode, silver_neutrality$idh_codel_subtype))
fisher.test(table(silver_neutrality$evolution_mode, silver_neutrality$driver_stability))
fisher.test(table(silver_neutrality$evolution_mode, silver_neutrality$recurrence_location))
fisher.test(table(silver_neutrality$evolution_mode, silver_neutrality$received_tmz))
fisher.test(table(silver_neutrality$evolution_mode, silver_neutrality$tmz_std_care))
fisher.test(table(silver_neutrality$evolution_mode, silver_neutrality$received_rt))
fisher.test(table(silver_neutrality$evolution_mode, silver_neutrality$hypermutator_status))
fisher.test(table(silver_neutrality$evolution_mode, silver_neutrality$grade_change))
fisher.test(table(silver_neutrality$evolution_mode, silver_neutrality$case_sex))

# IDHwt
fisher.test(table(silver_neutrality_IDHwt$evolution_mode, silver_neutrality_IDHwt$driver_status))
fisher.test(table(silver_neutrality_IDHwt$evolution_mode, silver_neutrality_IDHwt$driver_stability))
fisher.test(table(silver_neutrality_IDHwt$evolution_mode, silver_neutrality_IDHwt$recurrence_location))
fisher.test(table(silver_neutrality_IDHwt$evolution_mode, silver_neutrality_IDHwt$received_tmz))
fisher.test(table(silver_neutrality_IDHwt$evolution_mode, silver_neutrality_IDHwt$tmz_std_care))
fisher.test(table(silver_neutrality_IDHwt$evolution_mode, silver_neutrality_IDHwt$received_rt))
fisher.test(table(silver_neutrality_IDHwt$evolution_mode, silver_neutrality_IDHwt$hypermutator_status))
fisher.test(table(silver_neutrality_IDHwt$evolution_mode, silver_neutrality_IDHwt$grade_change))

# Model of evolution at recurrence:
fit_all_recur <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ recurrence_evolution,
               data = silver_neutrality)
ggsurvplot(fit_all_recur, data = silver_neutrality, risk.table = TRUE, pval= TRUE)
fit_codel_recur <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ recurrence_evolution,
                         data = silver_neutrality_IDHmut_codel)
ggsurvplot(fit_codel_recur, data = silver_neutrality_IDHmut_codel, risk.table = TRUE, pval= TRUE)
fit_noncodel_recur <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ recurrence_evolution,
                           data = silver_neutrality_IDHmut_noncodel)
ggsurvplot(fit_noncodel_recur, data = silver_neutrality_IDHmut_noncodel, risk.table = TRUE, pval= TRUE)
fit_wt_recur <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ recurrence_evolution,
                              data = silver_neutrality_IDHwt)
ggsurvplot(fit_wt_recur, data = silver_neutrality_IDHwt, risk.table = TRUE, pval= TRUE)


# Mode of evolution:
fit_all_mode <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ evolution_mode,
                data = silver_neutrality)
ggsurvplot(fit_all_mode, data = silver_neutrality, risk.table = TRUE, pval= TRUE) 
fit_codel_mode <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ evolution_mode,
                           data = silver_neutrality_IDHmut_codel)
ggsurvplot(fit_codel_mode, data = silver_neutrality_IDHmut_codel, risk.table = TRUE, pval= TRUE)
fit_noncodel_mode <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ evolution_mode,
                              data = silver_neutrality_IDHmut_noncodel)
ggsurvplot(fit_noncodel_mode, data = silver_neutrality_IDHmut_noncodel, risk.table = TRUE, pval= TRUE)
fit_wt_mode <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ evolution_mode,
                        data = silver_neutrality_IDHwt)
ggsurvplot(fit_wt_mode, data = silver_neutrality_IDHwt, risk.table = TRUE, pval= TRUE)

# RECURRENCE WITH SURGICAL INTERVAL.
fit_all_interval <- survfit(Surv(surgical_interval, rep(1, length(silver_neutrality$surgical_interval))) ~ recurrence_evolution,
                            data = silver_neutrality)
ggsurvplot(fit_all_interval, data = silver_neutrality, risk.table = TRUE, pval= TRUE, pval.method = TRUE, ylab = "Proportion Progression Free")
fit_codel_interval <- survfit(Surv(surgical_interval, rep(1, length(silver_neutrality_IDHmut_codel$surgical_interval))) ~ recurrence_evolution,
                              data = silver_neutrality_IDHmut_codel)
ggsurvplot(fit_codel_interval, data = silver_neutrality_IDHmut_codel, risk.table = TRUE, pval= TRUE, pval.method = TRUE, ylab = "Proportion Progression Free")
fit_noncodel_interval <- survfit(Surv(surgical_interval, rep(1, length(silver_neutrality_IDHmut_noncodel$surgical_interval))) ~ recurrence_evolution,
                                 data = silver_neutrality_IDHmut_noncodel)
ggsurvplot(fit_noncodel_interval, data = silver_neutrality_IDHmut_noncodel, risk.table = TRUE, pval= TRUE, pval.method = TRUE, ylab = "Proportion Progression Free")
fit_IDHwt_interval <- survfit(Surv(surgical_interval, rep(1, length(silver_neutrality_IDHwt$surgical_interval))) ~ recurrence_evolution,
                              data = silver_neutrality_IDHwt)
ggsurvplot(fit_IDHwt_interval, data = silver_neutrality_IDHwt, risk.table = TRUE, pval= TRUE, pval.method = TRUE, ylab = "Proportion Progression Free")

# EVOLUTION MODE WITH SURGICAL INTERVAL.
fit_all_interval <- survfit(Surv(surgical_interval, rep(1, length(silver_neutrality$surgical_interval))) ~ evolution_mode,
                            data = silver_neutrality)
ggsurvplot(fit_all_interval, data = silver_neutrality, risk.table = TRUE, pval= TRUE, pval.method = TRUE, ylab = "Proportion Progression Free")
fit_codel_interval <- survfit(Surv(surgical_interval, rep(1, length(silver_neutrality_IDHmut_codel$surgical_interval))) ~ evolution_mode,
                              data = silver_neutrality_IDHmut_codel)
ggsurvplot(fit_codel_interval, data = silver_neutrality_IDHmut_codel, risk.table = TRUE, pval= TRUE, pval.method = TRUE, ylab = "Proportion Progression Free")
fit_noncodel_interval <- survfit(Surv(surgical_interval, rep(1, length(silver_neutrality_IDHmut_noncodel$surgical_interval))) ~ evolution_mode,
                                 data = silver_neutrality_IDHmut_noncodel)
ggsurvplot(fit_noncodel_interval, data = silver_neutrality_IDHmut_noncodel, risk.table = TRUE, pval= TRUE, pval.method = TRUE, ylab = "Proportion Progression Free")
fit_IDHwt_interval <- survfit(Surv(surgical_interval, rep(1, length(silver_neutrality_IDHwt$surgical_interval))) ~ evolution_mode,
                              data = silver_neutrality_IDHwt)
ggsurvplot(fit_IDHwt_interval, data = silver_neutrality_IDHwt, risk.table = TRUE, pval= TRUE, pval.method = TRUE, ylab = "Proportion Progression Free")

surgeries_test = surgeries %>% 
  select(sample_barcode, surgery_number)

test = silver_neutrality %>% 
  mutate(sample_barcode_a = substr(tumor_barcode_a.x, 1, 15),
         sample_barcode_b = substr(tumor_barcode_b.x, 1, 15)) %>% 
left_join(surgeries_test, by=c("sample_barcode_a" = "sample_barcode")) %>% 
left_join(surgeries_test, by=c("sample_barcode_b" = "sample_barcode")) %>% 
  mutate(surgery_pair = paste(surgery_number.x, surgery_number.y, sep="-"))
table(test$surgery_pair)

full_test = silver_set %>% 
  mutate(sample_barcode_a = substr(tumor_barcode_a, 1, 15),
         sample_barcode_b = substr(tumor_barcode_b, 1, 15)) %>% 
  left_join(surgeries_test, by=c("sample_barcode_a" = "sample_barcode")) %>% 
  left_join(surgeries_test, by=c("sample_barcode_b" = "sample_barcode")) %>% 
  mutate(surgery_pair = paste(surgery_number.x, surgery_number.y, sep="-")) %>% 
  left_join(clinal_tumor_pairs, by="tumor_pair_barcode") %>% 
  mutate(case_barcode = substr(tumor_pair_barcode, 1, 12)) %>% 
  left_join(case_level_subtype, by="case_barcode") %>% 
  filter(received_tmz == 1, surgery_pair%in%c("1-2", "2-3","3-4"))
table(full_test$surgery_pair, full_test$hypermutator_status, full_test$idh_codel_subtype)

full_test_IDHwt = full_test %>% 
  filter(idh_codel_subtype == "IDHwt_noncodel")
full_test_IDH_noncodel = full_test %>% 
  filter(idh_codel_subtype == "IDHmut_noncodel")
full_test_IDH_codel = full_test %>% 
  filter(idh_codel_subtype == "IDHmut_codel")
wilcox.test(full_test_IDHwt$surgical_interval~full_test_IDHwt$hypermutator_status)
wilcox.test(full_test_IDH_noncodel$surgical_interval~full_test_IDH_noncodel$hypermutator_status)
wilcox.test(full_test_IDH_codel$surgical_interval~full_test_IDH_codel$hypermutator_status)


test_sub = test %>% 
  filter(surgery_pair == "1-2")
test_all_interval <- survfit(Surv(surgical_interval, rep(1, length(test_sub$surgical_interval))) ~ recurrence_evolution,
                            data = test_sub)
ggsurvplot(test_all_interval, data = test_sub, risk.table = TRUE, pval= TRUE, pval.method = TRUE, ylab = "Proportion Progression Free")
test_all_survival <- survfit(Surv(case_overall_survival_mo, rep(1, length(test_sub$case_overall_survival_mo))) ~ recurrence_evolution,
                             data = test_sub)
ggsurvplot(test_all_survival, data = test_sub, risk.table = TRUE, pval= TRUE, pval.method = TRUE)

 ### DRIVER SWITCH


silver_drivers_tmz = drivers %>% 
  inner_join(silver_set, by="tumor_pair_barcode") %>% 
  left_join(clinal_tumor_pairs, by = "tumor_pair_barcode") %>% 
  left_join(case_level_subtype, by="case_barcode") %>% 
  left_join(cases, by="case_barcode") %>% 
  filter(received_tmz==1) %>% 
  mutate(post_recurrence_surv = case_overall_survival_mo-surgical_interval)
  
silver_drivers_IDHwt = silver_drivers_tmz %>% 
  filter(idh_codel_subtype %in% c("IDHwt_noncodel"))
silver_drivers_IDHmut_noncodel = silver_drivers_tmz %>% 
  filter(idh_codel_subtype %in% c("IDHmut_noncodel"))
silver_drivers_IDHmut_codel = silver_drivers_tmz %>% 
  filter(idh_codel_subtype %in% c("IDHmut_codel"))

silver_drivers %>% 
  group_by(idh_codel_subtype, driver_stability) %>% 
  summarise(med_int = median(surgical_interval, na.rm = T),
            sample_size = n())

silver_drivers_tmz %>% 
  group_by(idh_codel_subtype, hypermutator_status) %>% 
  summarise(med_int = median(surgical_interval, na.rm = T),
            sample_size = n())
# Examine hypermutation:
wilcox.test(silver_drivers_IDHwt$surgical_interval~silver_drivers_IDHwt$hypermutator_status)
wilcox.test(silver_drivers_IDHmut_noncodel$surgical_interval~silver_drivers_IDHmut_noncodel$hypermutator_status)
wilcox.test(silver_drivers_IDHmut_codel$surgical_interval~silver_drivers_IDHmut_codel$hypermutator_status)

tmz_post_recur <- survfit(Surv(post_recurrence_surv, patient_vital) ~ silver_drivers_IDHwt,
                         data = silver_neutrality)
ggsurvplot(fit_all_recur, data = silver_neutrality, risk.table = TRUE, pval= TRUE)
fit_codel_recur <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ recurrence_evolution,
                           data = silver_neutrality_IDHmut_codel)



kruskal.test(silver_drivers$surgical_interval, as.factor(silver_drivers$driver_stability))
kruskal.test(silver_drivers_IDHwt$surgical_interval, as.factor(silver_drivers_IDHwt$driver_stability))
kruskal.test(silver_drivers_IDHmut_noncodel$surgical_interval, as.factor(silver_drivers_IDHmut_noncodel$driver_stability))
kruskal.test(silver_drivers_IDHmut_codel$surgical_interval, as.factor(silver_drivers_IDHmut_codel$driver_stability))

fisher.test(table(silver_drivers$driver_stability, silver_drivers$idh_codel_subtype))
fisher.test(table(silver_drivers$driver_stability, silver_drivers$recurrence_location))
fisher.test(table(silver_drivers$driver_stability, silver_drivers$received_tmz))
fisher.test(table(silver_drivers$driver_stability, silver_drivers$received_rt))
fisher.test(table(silver_drivers$driver_stability, silver_drivers$hypermutator_status))
fisher.test(table(silver_drivers$driver_stability, silver_drivers$grade_change))
fisher.test(table(silver_drivers$driver_stability, silver_drivers$case_sex))




# Driver stability:
fit_all_stability <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ driver_stability,
                         data = silver_neutrality)
ggsurvplot(fit_all_stability, data = silver_neutrality, risk.table = TRUE, pval= TRUE)
fit_codel_stability <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ driver_stability,
                           data = silver_neutrality_IDHmut_codel)
ggsurvplot(fit_codel_stability, data = silver_neutrality_IDHmut_codel, risk.table = TRUE, pval= TRUE)
fit_noncodel_stability <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ driver_stability,
                              data = silver_neutrality_IDHmut_noncodel)
ggsurvplot(fit_noncodel_stability, data = silver_neutrality_IDHmut_noncodel, risk.table = TRUE, pval= TRUE)
fit_wt_stability <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ driver_stability,
                        data = silver_neutrality_IDHwt)
ggsurvplot(fit_wt_stability, data = silver_neutrality_IDHwt, risk.table = TRUE, pval= TRUE)




kruskal.test(as.numeric(silver_neutrality_IDHmut_noncodel$surgical_interval), as.factor(silver_neutrality_IDHmut_noncodel$recurrence_evolution))
fisher.test(table(silver_neutrality_IDHmut_noncodel$driver_stability, silver_neutrality_IDHmut_noncodel$recurrence_evolution))
fisher.test(table(silver_neutrality_IDHmut_noncodel$driver_stability, silver_neutrality_IDHmut_noncodel$evolution_mode))







##############################
# Silver set visualizations  #
##############################
stacked_driver = silver_neutrality %>% 
  group_by(idh_codel_subtype, driver_stability) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n))
ggplot(stacked_driver, aes(x=idh_codel_subtype, y=freq, fill=driver_stability)) + geom_bar(stat="identity") +
  labs(fill="Driver change") + xlab("") + theme_bw()

stacked_neutrality = silver_neutrality %>% 
  group_by(idh_codel_subtype, evolution_mode) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n))
ggplot(stacked_neutrality, aes(x=idh_codel_subtype, y=freq, fill=evolution_mode)) + geom_bar(stat="identity") +
  labs(fill="Evolution route") + xlab("") + theme_bw() + ylab("Proportion of silver set pairs")

stacked_neutrality_recur = silver_neutrality %>% 
  group_by(idh_codel_subtype, recurrence_evolution) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n))
ggplot(stacked_neutrality_recur, aes(x=idh_codel_subtype, y=freq, fill=recurrence_evolution)) + geom_bar(stat="identity") +
  labs(fill="Evolution status at recurrence") + xlab("") + theme_bw() + ylab("Proportion of silver set pairs")


#################################
# Examine results at case-level
################################
# database query:
clinal_tumor_pairs = dbGetQuery(con, "WITH t1 AS
(
    SELECT cs.case_barcode, aliquot_analysis_type, string_agg(surgery_number::character(1), '-') OVER w AS surgeries, string_agg(evolution, '-') OVER w AS evolution, row_number() OVER w2 AS rnk
    FROM analysis.neutrality_aliquots an
    LEFT JOIN biospecimen.aliquots ba ON ba.aliquot_barcode = an.aliquot_barcode
    LEFT JOIN clinical.surgeries cs ON cs.sample_barcode = ba.sample_barcode
    WINDOW w AS (PARTITION BY case_barcode, aliquot_analysis_type ORDER BY surgery_number), w2 AS (PARTITION BY case_barcode, aliquot_analysis_type ORDER BY surgery_number DESC)
)
SELECT * FROM t1 WHERE rnk = 1")


# Separate out samples into IDHmut and IDHwt groups.
sample_neutrality_IDHmutcodel = per_sample_neutrality %>% 
  filter(idh_codel_subtype %in% c("IDHmut_codel"))
sample_neutrality_IDHmutnoncodel = per_sample_neutrality %>% 
  filter(idh_codel_subtype %in% c("IDHmut_noncodel"))
sample_neutrality_IDHwt = per_sample_neutrality %>% 
  filter(idh_codel_subtype %in% c("IDHwt_noncodel"))

# Note that this is a subset of the data.
drivers_collapsed = drivers %>% 
  left_join(clinal_tumor_pairs, by="tumor_pair_barcode") %>% 
mutate(driver_status_revised = recode(driver_status, "Driver gain" = "Driver unstable",
                                      "Driver loss" = "Driver unstable", "Driver switch" = "Driver unstable", "Driver null" = "Driver stable"))

# Analyses pertaining to all available driver information.
fisher.test(table(drivers_collapsed$driver_status, drivers_collapsed$recurrence_location))
fisher.test(table(drivers_collapsed$driver_status_revised, drivers_collapsed$recurrence_location))
fisher.test(table(drivers_collapsed$driver_status, drivers_collapsed$received_rt))
fisher.test(table(drivers_collapsed$driver_status_revised, drivers_collapsed$received_rt))
# Driver switch is strongly associated with hypermutation event (Obviously). 
fisher.test(table(drivers_collapsed$driver_status_revised, drivers_collapsed$hypermutator_status))


# Explore some associations with RECURRENCE EVOLUTION:
fisher.test(table(per_sample_neutrality$recurrence_evolution, per_sample_neutrality$idh_codel_subtype))
fisher.test(table(per_sample_neutrality$recurrence_evolution, per_sample_neutrality$driver_status))
fisher.test(table(per_sample_neutrality$recurrence_evolution, per_sample_neutrality$driver_status_revised))
fisher.test(table(per_sample_neutrality$recurrence_evolution, per_sample_neutrality$recurrence_location))
fisher.test(table(per_sample_neutrality$recurrence_evolution, per_sample_neutrality$received_tmz))
fisher.test(table(per_sample_neutrality$recurrence_evolution, per_sample_neutrality$tmz_std_care))
# RT use seems to be associated with selection.
fisher.test(table(per_sample_neutrality$recurrence_evolution, per_sample_neutrality$received_rt))
fisher.test(table(per_sample_neutrality$recurrence_evolution, per_sample_neutrality$hypermutator_status))
fisher.test(table(per_sample_neutrality$recurrence_evolution, per_sample_neutrality$grade_change))
# Neutral or Selected at recurrence was associated with surgical interval.
kruskal.test(as.numeric(per_sample_neutrality$surgical_interval), as.factor(per_sample_neutrality$recurrence_evolution))
## IDH wild-type neutrality differences still seem to have an impact on survival.
kruskal.test(as.numeric(sample_neutrality_IDHwt$surgical_interval), as.factor(sample_neutrality_IDHwt$recurrence_evolution))
kruskal.test(as.numeric(sample_neutrality_IDHmutnoncodel$surgical_interval), as.factor(sample_neutrality_IDHmutnoncodel$recurrence_evolution))
kruskal.test(as.numeric(sample_neutrality_IDHmutcodel$surgical_interval), as.factor(sample_neutrality_IDHmutcodel$recurrence_evolution))
# This might explain the overall relationship with surgical interval.
ggplot(per_sample_neutrality, aes(x=as.factor(recurrence_evolution), y=surgical_interval)) + geom_boxplot() + theme_bw() + ylab("Surgical interval") +
  xlab("Evolution (tumor_a - tumor_b)") + facet_grid(~idh_codel_subtype, scales="free")
ggplot(sample_neutrality_IDHwt, aes(x=as.factor(recurrence_evolution), y=surgical_interval)) + geom_boxplot() + theme_bw() + ylab("Surgical interval") +
  xlab("IDHwt Evolution (tumor_a - tumor_b)") + geom_text(x = 1.5, y = 65, label="Kruskal-Wallis P = 6.0E-03" )


# Explore some associations with MODE of EVOLUTION:
fisher.test(table(per_sample_neutrality$evolution_mode, per_sample_neutrality$driver_status_revised))
fisher.test(table(per_sample_neutrality$evolution_mode, per_sample_neutrality$recurrence_location))
fisher.test(table(per_sample_neutrality$evolution_mode, per_sample_neutrality$received_tmz))
fisher.test(table(per_sample_neutrality$evolution_mode, per_sample_neutrality$tmz_std_care))
# RT use seems to be associated with selection.
fisher.test(table(per_sample_neutrality$evolution_mode, per_sample_neutrality$received_rt))
fisher.test(table(per_sample_neutrality$evolution_mode, per_sample_neutrality$hypermutator_status))
fisher.test(table(per_sample_neutrality$evolution_mode, per_sample_neutrality$grade_change))


## IDH wild-type neutrality differences still seem to have an impact on survival.
kruskal.test(as.numeric(sample_neutrality_IDHwt$surgical_interval), as.factor(sample_neutrality_IDHwt$evolution_mode))
kruskal.test(as.numeric(sample_neutrality_IDHmutnoncodel$surgical_interval), as.factor(sample_neutrality_IDHmutnoncodel$evolution_mode))
kruskal.test(as.numeric(sample_neutrality_IDHmutcodel$surgical_interval), as.factor(sample_neutrality_IDHmutcodel$evolution_mode))

# This might explain the overall relationship with surgical interval.
ggplot(per_sample_neutrality, aes(x=as.factor(evolution_mode), y=surgical_interval)) + geom_boxplot() + theme_bw() + ylab("Surgical interval") +
  xlab("Evolution (tumor_a - tumor_b)") + facet_grid(~idh_codel_subtype, scales="free")
ggplot(sample_neutrality_IDHwt, aes(x=as.factor(evolution_mode), y=surgical_interval)) + geom_boxplot() + theme_bw() + ylab("Surgical interval") +
  xlab("IDHwt Evolution (tumor_a - tumor_b)") + geom_text(x = 1.5, y = 65, label="" )



# No substantial difference in paired R-squared values.
wilcox.test(silver_neutrality$model_rsq.x, silver_neutrality$model_rsq.y, paired = T)

# Quickly compare the R-squared values.
plot_neutrality_rsq = silver_neutrality %>% 
  gather("rsq_time", "rsq_value", c(model_rsq.x, model_rsq.y), -tumor_pair_barcode)
ggplot(plot_neutrality_rsq, aes(x = rsq_time, y = rsq_value, group = tumor_pair_barcode)) +
  geom_line(linetype="solid", size=1, alpha = 0.1) + 
  geom_point(color="black", size=2) + theme_bw() + ylab("R-squared values") + xlab("") + geom_hline(yintercept = 0.98, linetype="dotted") +
  geom_text(x = 1, y = 0.7, label="P=0.45") + ggtitle("Mutect2, primary_all, recurrence_all (n=167 selected pairs)")

## Compare the number of subclonal mutations in a tumor with the mutational frequency in the recurrence.
# Synchronous (subclonal mutation rate in primary) vs metachronous mutations (recurrence only mutation rate) 
# between three subtypes.
head(aliquot_mutation_counts)

# Step 1: Define query for extracting mutational frequency data by mutation type.
mut_freq_tumor_type <- dbGetQuery(con, read_file("sql/mutation_freq_private_shared.sql"))


# Step 2: Define mutation counts by silver set and merge tumor_barcode_b with mut_freq table.
mut_counts_silver = silver_set %>% 
  left_join(aliquot_mutation_counts, by=c("tumor_barcode_a" = "aliquot_barcode")) %>% 
  left_join(aliquot_mutation_counts, by=c("tumor_barcode_b" = "aliquot_barcode")) %>% 
  select(tumor_pair_barcode:tumor_barcode_b, subclonal_mut_a = subclonal_mut.x, subclonal_mut_b = subclonal_mut.y) %>% 
  left_join(mut_freq, by=c("tumor_barcode_a" = "aliquot_barcode")) %>% 
  left_join(mut_freq, by=c("tumor_barcode_b" = "aliquot_barcode")) %>% 
  mutate_if(bit64::is.integer64, as.double) %>% 
  left_join(case_level_subtype, by="case_barcode") %>% 
  left_join(clinal_tumor_pairs, by="tumor_pair_barcode") %>% 
  select(-tumor_barcode_a.y, -tumor_barcode_b.y) %>% 
  mutate(hypermutant_manual = ifelse(coverage_adj_mut_freq.y > 10, "hyper", "non-hyper")) %>% 
  mutate(subclonal_proportion = subclonal_mut_a / mutation_count.x) %>% 
  left_join(mut_freq_tumor_type, by="tumor_pair_barcode")
  
# Step 3: Test whether number of subclonal mutations is associated with recurrence mutation freq.
ggplot(mut_counts_silver, aes(x=subclonal_proportion*mf_a, mf_b)) + geom_point() + geom_smooth(method='lm',formula=y~x) +
  facet_grid(hypermutant_manual~idh_codel_subtype, scales = "free") + theme_bw() + xlab("Primary subclonal mut. rate") + ylab("Private to recurrence mut. rate")
 
