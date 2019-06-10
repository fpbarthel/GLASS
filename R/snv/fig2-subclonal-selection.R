##############################################
# Analyze SubClonalSelection (Georgette Tanner) results for GLASS samples
# Updated: 2019.06.15
# Author: Kevin J.
##################################################

# Working directory for this analysis in the GLASS-analysis project. 
mybasedir = "/Users/johnsk/Documents/Life-History/GLASS-WG/" 
setwd(mybasedir)

#######################################################
# Necessary packages:
library(tidyverse)
library(DBI)
library(survminer)
library(survival)
library(alluvial)

#######################################################
# Establish connection with the second version of the GLASS database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB2")

# Load additional tables from the database.
subtypes = dbReadTable(con,  Id(schema="clinical", table="subtypes"))  
cases = dbReadTable(con,  Id(schema="clinical", table="cases"))  
tumor_pairs = dbReadTable(con,  Id(schema="analysis", table="tumor_pairs"))
silver_set = dbReadTable(con,  Id(schema="analysis", table="silver_set"))
gold_set = dbReadTable(con,  Id(schema="analysis", table="gold_set"))
subclonalselection = dbReadTable(con,  Id(schema="analysis", table="subclonalselection"))
titan_params = dbReadTable(con,  Id(schema="variants", table="titan_params"))

# Get the clinical_tumor_pairs table.
clinical_tumor_pairs_query = read_file("sql/clinical/clinical-tumor-pairs-db2.sql")
clinical_tumor_pairs <- dbGetQuery(con, clinical_tumor_pairs_query)

# Merge with gold_set for mode of evolution at recurrence.
gold_selection = gold_set %>% 
  left_join(subclonalselection, by=c("tumor_barcode_a"="aliquot_barcode")) %>% 
  left_join(subclonalselection, by=c("tumor_barcode_b"="aliquot_barcode")) %>% 
  left_join(cases, by="case_barcode") %>% 
  left_join(subtypes, by="case_barcode") %>% 
  mutate(patient_vital = ifelse(case_vital_status=="alive", 0, 1),
         # Re-classify "most_probably_classification" to "recurrence_threshold_0.5".
         r_thres_5 = most_probable_classification.y,
         r_thres_6 = ifelse(probability_neutral.y > 0.6, "N", ifelse(probability_neutral.y < 0.4, "S", NA)),
         r_thres_7 = ifelse(probability_neutral.y > 0.7, "N", ifelse(probability_neutral.y < 0.3, "S", NA))) %>% 
  left_join(clinical_tumor_pairs, by=c("tumor_pair_barcode", "case_barcode", "tumor_barcode_a", "tumor_barcode_b")) 

# Create subsets of the data by separating out subtypes.
IDHmut_codel_gold_selection = gold_selection %>% 
  filter(idh_codel_subtype=="IDHmut-codel")
IDHmut_noncodel_gold_selections = gold_selection %>% 
  filter(idh_codel_subtype=="IDHmut-noncodel")
IDHwt_gold_selection = gold_selection %>% 
  filter(idh_codel_subtype=="IDHwt")


###################
# IDHwt survival - selection
###################
# Restrict only to the IDHwt tumors since this is the main comparison presented in the figure of the manuscript
# Survival analysis taking the most probable threshold.
fit_wt_recurrence <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ r_thres_5,
                             data = IDHwt_gold_selection)
pdf(file = "/Users/johnsk/Documents/selection-R-thres5.pdf", height = 5, width = 7, bg = "transparent", useDingbats = FALSE)
ggsurvplot(fit_wt_recurrence, data = IDHwt_gold_selection, risk.table = FALSE, pval= TRUE, pval.coord = c(100, 0.75),
           palette = c("#A3BD67", "#BD8167"),
           ylab = "Overall survival \n probability", xlab = "Time (months)")
dev.off()

# Increase the minimum threshold for a sample to called either "Neutral" or "Selected".
fit_wt_recurrence <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ r_thres_6,
                             data = IDHwt_gold_selection)
pdf(file = "/Users/johnsk/Documents/selection-R-thres6.pdf", height = 5, width = 7, bg = "transparent", useDingbats = FALSE)
ggsurvplot(fit_wt_recurrence, data = IDHwt_gold_selection, risk.table = FALSE, pval= TRUE, pval.coord = c(100, 0.75),
           palette = c("#A3BD67", "#BD8167"),
           ylab = "Overall survival \n probability", xlab = "Time (months)")
dev.off()
# Further increase the threshold for a sample to be classified to P > 0.7.
fit_wt_recurrence <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ r_thres_7,
                             data = IDHwt_gold_selection)
pdf(file = "/Users/johnsk/Documents/selection-R-thres7.pdf", height = 5, width = 7, bg = "transparent", useDingbats = FALSE)
ggsurvplot(fit_wt_recurrence, data = IDHwt_gold_selection, risk.table = FALSE, pval= TRUE, pval.coord = c(100, 0.75),
           palette = c("#A3BD67", "#BD8167"),
           ylab = "Overall survival \n probability", xlab = "Time (months)")
dev.off()

# Perform Cox analysis using categorical variable of selection across groups.
neutrality_cox <- coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + idh_codel_subtype + r_thres_5, data = gold_selection)
summary(neutrality_cox)
neutrality_cox <- coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + idh_codel_subtype + r_thres_6, data = gold_selection)
summary(neutrality_cox)
neutrality_cox <- coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + idh_codel_subtype + r_thres_7, data = gold_selection)
summary(neutrality_cox)

# Perform Cox analysis using continuous variable of selection across groups.
neutrality_cox <- coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + idh_codel_subtype + probability_neutral.y, data = gold_selection)
summary(neutrality_cox) # Neutrality is protective across subtypes while adjusting for age.

###############################
# Investigate Evolution modes #
###############################
gold_selection_mode = gold_selection %>% 
  filter(!is.na(most_probable_classification.x), !is.na(most_probable_classification.y)) %>% 
  mutate(evo_mode = paste(most_probable_classification.x, most_probable_classification.y, sep="-")) %>% 
  mutate(evo_binary = ifelse(evo_mode=="N-N", "Neutral-Neutral", "N-S|S-N|S-S")) 

table(gold_selection_mode$evo_mode)

# Define "evolution mode", selection status at primary and recurrence is necessary for this information.
IDHwt_gold_selection_mode = gold_selection_mode %>% 
  filter(idh_codel_subtype=="IDHwt")
# Perform survival KM survival analysis for evolution modes in IDHwt tumors.
fit_wt_recurrence <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ evo_mode,
                             data = IDHwt_gold_selection_mode)
ggsurvplot(fit_wt_recurrence, data = IDHwt_gold_selection_mode, risk.table = FALSE, pval= TRUE, pval.coord = c(100, 0.75),
           ylab = "Overall survival \n probability", xlab = "Time (months)")

# Across all subtypes is there a relationship between evolution mode and overall survival.
neutrality_cox <- coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + idh_codel_subtype + evo_mode, data = gold_selection_mode)
summary(neutrality_cox)

# Created a binary version of N-N versus all other categories.
fit_wt_recurrence <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ evo_binary,
                             data = IDHwt_gold_selection_mode)
ggsurvplot(fit_wt_recurrence, data = IDHwt_gold_selection_mode, risk.table = TRUE, pval= TRUE, pval.coord = c(100, 0.75),
           ylab = "Overall survival \n probability", xlab = "Time (months)")

# Perform survival analysis with the binarized evolution modes.
neutrality_cox <- coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + idh_codel_subtype + evo_binary, data = gold_selection_mode)
summary(neutrality_cox)

#######################################
### Create selection Sankey plots  ##
#######################################
## Color neutral-neutral or selection-selection.
# Use subtype, driver stability, and evolution pattern.
gold_select_dat = gold_selection %>% 
  filter(!is.na(most_probable_classification.x), !is.na(most_probable_classification.y)) %>% 
  mutate(evo_mode = paste(most_probable_classification.x, most_probable_classification.y, sep="-")) %>% 
  left_join(clinical_tumor_pairs, by=c("tumor_pair_barcode", "case_barcode", "tumor_barcode_a", "tumor_barcode_b"))

neutral_sankey = gold_select_dat %>% 
  select(tumor_pair_barcode, idh_codel_subtype, hypermutator_status = hypermutator_status.y, received_alk = received_alk.y, received_rt = received_rt.y, treatment = received_treatment.y, primary_evolution = most_probable_classification.x, recurrence_evolution = most_probable_classification.y) %>% 
  mutate(primary_evolution = recode(primary_evolution, "N" = "Neutral", "S" = "Selected"),
         recurrence_evolution = recode(recurrence_evolution, "N" = "Neutral", "S" = "Selection"),
         treatment = ifelse(is.na(treatment), "Unknown", treatment),
         treatment = recode(treatment, "1" = "Yes", "0" = "No"))  %>% 
  group_by(idh_codel_subtype, primary_evolution, treatment, recurrence_evolution) %>% 
  summarise(Freq = n())

# Adjust the alluvial plot so that different colors map to different selection strength.
# Color scheme for evolution trajectory.
#619CFF - IDHwt
#00BA38 - IDHmut noncodel
#F8766D - IDHmut codel
pal = c("#619CFF", "#00BA38", "#F8766D", "lightgray")
neutral_sankey$col_type = ifelse(neutral_sankey$idh_codel_subtype == "IDHwt", pal[1], ifelse(neutral_sankey$idh_codel_subtype == "IDHmut-noncodel", pal[2], pal[3]))
neutral_sankey$col_type = ifelse(neutral_sankey$treatment == "Unknown", pal[4], neutral_sankey$col_type)

pdf(file = "/Users/johnsk/Documents/f2b-kcj.pdf", height = 4, width = 8, bg = "transparent", useDingbats = FALSE)
alluvial(neutral_sankey[,1:4], freq=neutral_sankey$Freq, border=NA,
         col= neutral_sankey$col_type,
         cw = 0.25, axis_labels = c("Glioma subtype", "Primary tumor", "Treatment", "Recurrent tumor"))
dev.off()

# Create legend for alluvial plot.
pdf(file = "/Users/johnsk/Documents/f2b-legend-kcj.pdf", height = 4, width = 8, bg = "transparent", useDingbats = FALSE)
plot.new()
legend("center", title="Glioma subtype",
       c("IDHwt","IDHmut-noncodel","IDHmut-codel", "Treatment unknown"), fill=pal, cex=1.5, bty = "n")
dev.off()


####################################################
### Association between treatment and selection   ##
####################################################
# Test by method of **Selection at Recurrence**
## IDH MUT CODELs. The codels have a very small sample size and may error out.
## Only 3 selected samples. Not worth performing association tests.
fisher.test(table(IDHmut_codel_gold_selection$r_thres_5, IDHmut_codel_gold_selection$received_tmz))
fisher.test(table(IDHmut_codel_gold_selection$r_thres_5, IDHmut_codel_gold_selection$received_rt))

## IDH MUT NONCODELs. Neutral (n=32) and Selection (n=11).
wilcox.test(IDHmut_noncodel_gold_selections$surgical_interval~IDHmut_noncodel_gold_selections$r_thres_5)
fisher.test(table(IDHmut_noncodel_gold_selections$r_thres_5, IDHmut_noncodel_gold_selections$received_tmz))
fisher.test(table(IDHmut_noncodel_gold_selections$r_thres_5, IDHmut_noncodel_gold_selections$received_rt))
fisher.test(table(IDHmut_noncodel_gold_selections$r_thres_5, IDHmut_noncodel_gold_selections$received_treatment))
# Selection is found in 5 / 8 Hypermutated cases, but only 6 out of 35 non-hypermutators. 
# P = 0.01
# Create an image that may be useful moving forward.
fisher.test(table(IDHmut_noncodel_gold_selections$r_thres_5, IDHmut_noncodel_gold_selections$hypermutator_status))
stacked_selection = IDHmut_noncodel_gold_selections %>% 
  filter(!is.na(r_thres_5)) %>% 
  group_by(hypermutator_status, r_thres_5) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n))

pdf(file = "/Users/johnsk/Documents/f1c-kcj.pdf", height = 5, width = 7, bg = "transparent", useDingbats = FALSE)
ggplot(stacked_selection, aes(x=hypermutator_status, y=freq, fill=r_thres_5)) + geom_bar(stat="identity") + scale_fill_manual(values=c("royalblue4", "tomato3")) +
  labs(fill="Selection and hypermutation") + xlab("") + theme_bw() + ylab("Proportion of recurrent\nsamples selected")
dev.off()

fisher.test(table(IDHmut_noncodel_gold_selections$r_thres_5, IDHmut_noncodel_gold_selections$grade_change))

## IDH WT.
wilcox.test(IDHwt_gold_selection$surgical_interval~IDHwt_gold_selection$r_thres_5)
fisher.test(table(IDHwt_gold_selection$r_thres_5, IDHwt_gold_selection$recurrence_location))
fisher.test(table(IDHwt_gold_selection$r_thres_5, IDHwt_gold_selection$received_tmz))
fisher.test(table(IDHwt_gold_selection$r_thres_5, IDHwt_gold_selection$received_rt))
fisher.test(table(IDHwt_gold_selection$r_thres_5, IDHwt_gold_selection$received_treatment))
fisher.test(table(IDHwt_gold_selection$r_thres_5, IDHwt_gold_selection$hypermutator_status))
fisher.test(table(IDHwt_gold_selection$r_thres_5, IDHwt_gold_selection$grade_change))



