##############################################
# Examine the survival associations with hypermutation.
# Updated: 2019.05.14
# Author: Kevin J.
##################################################

# Working directory for this analysis in the GLASS-analysis project. 
mybasedir = "/Volumes/verhaak-lab/GLASS-analysis/"
setwd(mybasedir)

#######################################################
# Necessary packages:
library(tidyverse)
library(DBI)
library(survminer)
library(survival)
library(RColorBrewer)

######################################################## 
# Establish connection with the GLASS database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB2")

# Load additional tables from the database.
pairs = dbReadTable(con,  Id(schema="analysis", table="pairs"))
tumor_pairs = dbReadTable(con,  Id(schema="analysis", table="tumor_pairs"))
aliquots = dbReadTable(con,  Id(schema="biospecimen", table="aliquots"))
surgeries = dbReadTable(con,  Id(schema="clinical", table="surgeries"))
cases = dbReadTable(con,  Id(schema="clinical", table="cases"))
titan_params = dbReadTable(con,  Id(schema="variants", table="titan_params"))
mut_freq = dbReadTable(con,  Id(schema="analysis", table="mut_freq"))
silver_set = dbReadTable(con,  Id(schema="analysis", table="silver_set"))
gold_set = dbReadTable(con,  Id(schema="analysis", table="gold_set"))
subtypes = dbReadTable(con,  Id(schema="clinical", table="subtypes"))

# Get the new clinical_tumor_pairs table.
clinical_tumor_pairs_query = read_file("sql/clinical/clinical-tumor-pairs-db2.sql")
clinical_tumor_pairs <- dbGetQuery(con, clinical_tumor_pairs_query)

# Restrict the analyses to the gold set.
clinical_tumor_gold = clinical_tumor_pairs %>%  
  filter(tumor_pair_barcode %in% gold_set$tumor_pair_barcode) %>% 
  left_join(subtypes, by="case_barcode") %>% 
  left_join(cases, by="case_barcode") %>% 
  left_join(mut_freq, by=c("tumor_barcode_b" = "aliquot_barcode")) %>% 
  mutate(patient_vital =  ifelse(case_vital_status=="alive", 0, 1))

# Is there a difference between the hypermutants survival versus all others? Answer: No.
hyper_cox_model = coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + hypermutator_status + idh_codel_subtype, data = clinical_tumor_gold)
summary(hyper_cox_model)

# What about a difference between the alkylating agent associated hypermutants survival versus all others treated with alkylating agents? Answer: No.
hyper_cox_model = coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + alk_assoc_hypermutator_status + idh_codel_subtype, data = clinical_tumor_gold)
summary(hyper_cox_model)

# What if the hypermutation status changes? 2 standard deviations (12.4 muts/Mb) or Campbell et al. Cell 2017.
clinical_tumor_gold_2sd = clinical_tumor_gold %>% 
  filter(!is.na(alk_assoc_hypermutator_status)) %>% 
  mutate(two_sd_hyper = ifelse(alk_assoc_hypermutator_status == 1 & coverage_adj_mut_freq > 12.4, 1, 0))
hyper_cox_model = coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + two_sd_hyper + idh_codel_subtype, data = clinical_tumor_gold_2sd)
summary(hyper_cox_model)

# Ultra-hypermutation as defined by Campbell et al. Cell 2017.
clinical_tumor_gold_ultra = clinical_tumor_gold %>% 
  filter(!is.na(alk_assoc_hypermutator_status)) %>% 
  mutate(pancan_ultra = ifelse(alk_assoc_hypermutator_status == 1 & coverage_adj_mut_freq > 100, 1, 0))
hyper_cox_model = coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + pancan_ultra + idh_codel_subtype, data = clinical_tumor_gold_ultra)
summary(hyper_cox_model)

# Again, pay attention that we are restricting these analyses to TMZ-alkylating agent treated samples. 
stacked_tmz = clinical_tumor_gold %>% 
  filter(!is.na(alk_assoc_hypermutator_status)) %>% 
  group_by(idh_codel_subtype, alk_assoc_hypermutator_status) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n))

pdf(file = "/Users/johnsk/Documents/f1c-kcj.pdf", height = 5, width = 7, bg = "transparent", useDingbats = FALSE)
ggplot(stacked_tmz, aes(x=idh_codel_subtype, y=freq, fill=alk_assoc_hypermutator_status)) + geom_bar(stat="identity") + scale_fill_manual(values=c("royalblue4", "tomato3")) +
  labs(fill="Alkylating agent\nassociated hypermutation") + xlab("") + theme_bw() + ylab("")
dev.off()
######################
##### Survival #######
######################
### IDHwt
# Overall survival for IDHwt tumors treated with alkylating agent (n=99). 
clinical_tumor_gold_IDHwt = clinical_tumor_gold %>% 
  filter(idh_codel_subtype == "IDHwt") %>% 
  filter(received_alk == 1)
  
fit_IDHwt_surv <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ alk_assoc_hypermutator_status,
                          data = clinical_tumor_gold_IDHwt)
pdf(file = "/Users/johnsk/Documents/f1d-idhwt-kcj.pdf", height = 5, width = 7, bg = "transparent", useDingbats = FALSE)
ggsurvplot(fit_IDHwt_surv, data = clinical_tumor_gold_IDHwt, risk.table = FALSE, pval= TRUE, pval.coord = c(100, 0.50),
           palette = c("royalblue4", "tomato3"), 
           surv.median.line = "v", ylab = "Overall survival \n probability", xlab = "Time (months)") 
dev.off()
# Does this change when age is included as a variable? Answer: No.
hyper_cox_model = coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + alk_assoc_hypermutator_status, data = clinical_tumor_gold_IDHwt)
summary(hyper_cox_model) # P = 0.93


### IDHmut-noncodel
# Overall survival for IDHmut-noncodel tumors treated with alkylating agent (n=32). 
clinical_tumor_gold_IDHmut_noncodel = clinical_tumor_gold %>% 
  filter(idh_codel_subtype == "IDHmut-noncodel") %>% 
  filter(received_alk == 1)

# Hypermutation is tenatively associated with a grade change (P = 0.05).
fisher.test(table(clinical_tumor_gold_IDHmut_noncodel$alk_assoc_hypermutator_status, clinical_tumor_gold_IDHmut_noncodel$grade_change))

# Restricting to alkylating-agent associated hypermutation.
fit_IDHmutnon_surv <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ alk_assoc_hypermutator_status,
                          data = clinical_tumor_gold_IDHmut_noncodel)
pdf(file = "/Users/johnsk/Documents/f1d-idhmut-kcj.pdf", height = 5, width = 7, bg = "transparent", useDingbats = FALSE)
ggsurvplot(fit_IDHmutnon_surv, data = clinical_tumor_gold_IDHmut_noncodel, risk.table = FALSE, pval= TRUE, pval.coord = c(100, 0.50),
           palette = c("royalblue4", "tomato3"), 
           surv.median.line = "v", ylab = "Overall survival \n probability", xlab = "Time (months)") 
dev.off()

# Does this change when age is included as a variable? Answer: No.
hyper_cox_model = coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + alk_assoc_hypermutator_status, data = clinical_tumor_gold_IDHmut_noncodel)
summary(hyper_cox_model) # P = 0.536

############################
#### Surgical interval #####
############################
# Restrict to only those sequential samples to reflect progression-free survival.
# As you can see in the section below including the small percentage of non-sequential samples does have an impact.
fit_IDHwt_int <- survfit(Surv(surgical_interval, rep(1, length(clinical_tumor_gold_IDHwt$surgical_interval))) ~ alk_assoc_hypermutator_status,
                          data = clinical_tumor_gold_IDHwt)
ggsurvplot(fit_IDHwt_int, data = clinical_tumor_gold_IDHwt, risk.table = FALSE, pval= TRUE, pval.method = TRUE) + ylab("Progression-Free probability")
wilcox.test(clinical_tumor_gold_IDHwt$surgical_interval~clinical_tumor_gold_IDHwt$alk_assoc_hypermutator_status)

fit_IDHmut_int <- survfit(Surv(surgical_interval, rep(1, length(clinical_tumor_gold_IDHmut_noncodel$surgical_interval))) ~ alk_assoc_hypermutator_status,
                         data = clinical_tumor_gold_IDHmut_noncodel)
ggsurvplot(fit_IDHmut_int, data = clinical_tumor_gold_IDHmut_noncodel, risk.table = TRUE, pval= TRUE, pval.method = TRUE,
           surv.median.line = "v", palette = c("royalblue4", "tomato3"))
wilcox.test(clinical_tumor_gold_IDHmut_noncodel$surgical_interval~clinical_tumor_gold_IDHmut_noncodel$alk_assoc_hypermutator_status)


# True progression-free survival (sequentially sampled tumor samples)
clinical_surgery_seq = clinical_tumor_gold %>% 
mutate(sample_barcode_a = substr(tumor_barcode_a, 1, 15),
       sample_barcode_b = substr(tumor_barcode_b, 1, 15)) %>% 
  left_join(surgeries, by=c("sample_barcode_a" = "sample_barcode")) %>% 
  left_join(surgeries, by=c("sample_barcode_b" = "sample_barcode")) %>% 
  mutate(surgery_pair = paste(surgery_number.x, surgery_number.y, sep="-")) %>% 
  filter(received_alk== 1, surgery_pair%in%c("1-2", "2-3","3-4"))

# Divide into IDHwt and IDHmut noncodel.
clinical_gold_IDHwt_seq = clinical_surgery_seq %>% 
  filter(idh_codel_subtype == "IDHwt") 
clinical_gold_IDHmut_seq = clinical_surgery_seq %>% 
  filter(idh_codel_subtype == "IDHmut-noncodel") 

# IDHwt - PFS. n = 89 (non-hypermutants = 76, hypermutants = 13).
fit_IDHwt_seq_int <- survfit(Surv(surgical_interval, rep(1, length(clinical_gold_IDHwt_seq$surgical_interval))) ~ alk_assoc_hypermutator_status,
                         data = clinical_gold_IDHwt_seq)
ggsurvplot(fit_IDHwt_seq_int, data = clinical_gold_IDHwt_seq, risk.table = FALSE, pval= TRUE, pval.method = FALSE, pval.coord = c(40, 0.50), 
           surv.median.line = "v", palette = c("royalblue4", "tomato3"), 
           ylab = "Progression-free survival \n probability")


# IDHmut noncodels. n = 24 (n = 13 (non-hyper), n = 11 (hyper)).
fit_IDHmut_seq_int <- survfit(Surv(surgical_interval, rep(1, length(clinical_gold_IDHmut_seq$surgical_interval))) ~ alk_assoc_hypermutator_status,
                             data = clinical_gold_IDHmut_seq)
ggsurvplot(fit_IDHmut_seq_int, data = clinical_gold_IDHmut_seq, risk.table = FALSE, pval= TRUE, pval.method = FALSE, pval.coord = c(40, 0.75),
           surv.median.line = "v", palette = c("royalblue4", "tomato3"), 
           ylab = "Progression-free survival \n probability")

#######################################################
# Inclusion of MGMT status for hypermutation analyses #
#######################################################
# Join the clinical annotation with the surgery table to extract MGMT status.
surgery_mgmt = surgeries %>% 
  filter(!is.na(mgmt_methylation)) %>% 
  select(sample_barcode, surgery_number, mgmt_methylation, mgmt_methylation_method)

clinical_tumor_gold_mgmt = clinical_tumor_gold %>% 
  mutate(sample_bar_a = substr(tumor_barcode_a, 1, 15),
         sample_bar_b = substr(tumor_barcode_b, 1, 15)) %>% 
  left_join(surgery_mgmt, by=c("sample_bar_a"="sample_barcode")) %>% 
  left_join(surgery_mgmt, by=c("sample_bar_b"="sample_barcode"))

# There is MGMT data available for roughly half the GLASS dataset.
sum(!is.na(clinical_tumor_gold_mgmt$mgmt_methylation.x))
# Far fewer samples that have MGMT methylation at recurrence.
sum(!is.na(clinical_tumor_gold_mgmt$mgmt_methylation.y))
sum(clinical_tumor_gold_mgmt$mgmt_methylation.x==clinical_tumor_gold_mgmt$mgmt_methylation.y, na.rm = TRUE)
# 12/75 samples that were discordant for MGMT status.
sum(clinical_tumor_gold_mgmt$mgmt_methylation.x!=clinical_tumor_gold_mgmt$mgmt_methylation.y, na.rm = TRUE)

# Is there an association between survival and hypermutation while including MGMT status?
hyper_cox_mgmt = coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + alk_assoc_hypermutator_status + idh_codel_subtype + mgmt_methylation.x, data = clinical_tumor_gold_mgmt)
summary(hyper_cox_mgmt) # Still no association with hypermutation, but there was an association with unmethylated MGMT.

