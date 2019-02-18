##############################################
# Examine the survival and PFS associations with hypermutation.
# Updated: 2019.02.05
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
library(RColorBrewer)

######################################################## 
#Establish connection with Floris' database.
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

# Combine all this information for clinical tumor pairs.
clinical_silver = silver_set %>% 
  left_join(cases, by="case_barcode") %>% 
  inner_join(clinal_tumor_pairs, by="tumor_pair_barcode") %>% 
  inner_join(all_drivers, by="tumor_pair_barcode") %>% 
  mutate(treatment = ifelse(received_tmz == 1 | received_rt == 1, "YES", "NO"),
         treatment = ifelse(is.na(treatment), NA, treatment))
clinical_silver$patient_vital = ifelse(clinical_silver$case_vital_status=="alive", 0, 1)
clinical_silver$patient_post_recur_surv = clinical_silver$case_overall_survival_mo-clinical_silver$surgical_interval

# Hypermutation is associated with a grade change.
fisher.test(table(clinical_silver$hypermutator_status, clinical_silver$grade_change))

### Generate stacked barplot for TMZ-alkylating agent treated samples.
clinical_silver_alkylating= clinical_silver %>% 
  filter(received_tmz == 1 | received_alkylating_agent== 1) %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, "IDHmut_codel" = "IDHmut codel", "IDHmut_noncodel" = "IDHmut noncodel", "IDHwt_noncodel" = "IDHwt"))
# Again, pay attention that we are restricting these analyses to TMZ-alkylating agent treated samples. 
stacked_tmz = clinical_silver_alkylating %>% 
  group_by(idh_codel_subtype, hypermutator_status) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n))
ggplot(stacked_tmz, aes(x=idh_codel_subtype, y=freq, fill=hypermutator_status)) + geom_bar(stat="identity") + scale_fill_manual(values=c("royalblue4", "tomato3")) +
  labs(fill="TMZ induced \n hypermutation") + xlab("") + theme_bw() + ylab("")

# Is there a difference in the proportion of hypermutators across these subtypes.
fisher.test(table(clinical_silver_alkylating$hypermutator_status, clinical_silver_alkylating$idh_codel_subtype))

# Restrict to only those TMZ or other alkylating agent-treated tumors. TCGA-06-0152 was treated with Carmustine. 
clinical_silver_IDHwt = clinical_silver %>% 
  filter(idh_codel_subtype == "IDHwt_noncodel") %>% 
  filter(received_tmz == 1 | received_alkylating_agent== 1)

# Just the non-codels since the codels only have two hypermutants
clinical_silver_IDHmut = clinical_silver %>% 
  filter(idh_codel_subtype == "IDHmut_noncodel") %>% 
  filter(received_tmz == 1 | received_alkylating_agent== 1)


### Survival, IDHwt ####
# Overall survival for IDHwt tumors treated with alkylating agent (n=96). 
fit_IDHwt_surv <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ hypermutator_status,
                             data = clinical_silver_IDHwt)
ggsurvplot(fit_IDHwt_surv, data = clinical_silver_IDHwt, risk.table = FALSE, pval= TRUE, pval.coord = c(100, 0.50),
           palette = c("royalblue4", "tomato3"), 
           surv.median.line = "v", ylab = "Overall survival \n probability") 
hyper_cox_model = coxph(Surv(case_overall_survival_mo, patient_vital) ~ hypermutator_status, data = clinical_silver_IDHwt)
summary(hyper_cox_model)

# Post-recurrence overall survival for IDHwt tumors.
fit_IDHwt_recur_surv <- survfit(Surv(patient_post_recur_surv, patient_vital) ~ hypermutator_status,
                          data = clinical_silver_IDHwt)
ggsurvplot(fit_IDHwt_recur_surv, data = clinical_silver_IDHwt, risk.table = FALSE, pval= TRUE, pval.coord = c(75, 0.50), 
           palette = c("royalblue4", "tomato3"), surv.median.line = "v", ylab = "Post-recurrence overall \n survival probability") 


### Survival, IDHmut ####
# Overall survival. 36 IDHmut-noncodels (n=20 non-hypermutators, n=16 hypermutators).
fit_IDHmut_surv <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ hypermutator_status,
                          data = clinical_silver_IDHmut)
ggsurvplot(fit_IDHmut_surv, data = clinical_silver_IDHmut, risk.table = FALSE, pval= TRUE, pval.method = FALSE, pval.coord = c(150, 0.50),
           surv.median.line = "v", palette = c("royalblue4", "tomato3"))

# Post-recurrence overall-survival.
fit_IDHmut_recur_surv <- survfit(Surv(patient_post_recur_surv, patient_vital) ~ hypermutator_status,
                          data = clinical_silver_IDHmut)
ggsurvplot(fit_IDHmut_recur_surv, data = clinical_silver_IDHmut, risk.table = FALSE, pval= TRUE, pval.method = TRUE, pval.coord = c(40, 0.50),
           surv.median.line = "v", palette = c("royalblue4", "tomato3"), ylab = "Post-recurrence overall \n survival probability") 

# Create-legend for survival plot.
ggsurvplot(fit_IDHmut_recur_surv, data = clinical_silver_IDHmut, risk.table = FALSE, pval= TRUE, pval.method = TRUE, pval.coord = c(40, 0.50),
           surv.median.line = "v", palette = c("royalblue4", "tomato3"), ylab = "Post-recurrence overall survival \n probability", 
           legend = "right", legend.labs =c("Non-hypermutatnt", "Hypermutant")) 


#### Progression #####
# Restrict to only those sequential samples to reflect progression-free survival.
# As you can see in the section below including the small percentage of non-sequential samples does have an impact.
fit_IDHwt_int <- survfit(Surv(surgical_interval, rep(1, length(clinical_silver_IDHwt$surgical_interval))) ~ hypermutator_status,
                          data = clinical_silver_IDHwt)
ggsurvplot(fit_IDHwt_int, data = clinical_silver_IDHwt, risk.table = FALSE, pval= TRUE, pval.method = TRUE) + ylab("Progression-Free probability")
wilcox.test(clinical_silver_IDHwt$surgical_interval~clinical_silver_IDHwt$hypermutator_status)

fit_IDHmut_int <- survfit(Surv(surgical_interval, rep(1, length(clinical_silver_IDHmut$surgical_interval))) ~ hypermutator_status,
                         data = clinical_silver_IDHmut)
ggsurvplot(fit_IDHmut_int, data = clinical_silver_IDHmut, risk.table = TRUE, pval= TRUE, pval.method = TRUE)
wilcox.test(clinical_silver_IDHmut$surgical_interval~clinical_silver_IDHmut$hypermutator_status)


# True progression-free survival (sequentially sampled tumor samples)
clinical_silver_seq = clinical_silver %>% 
mutate(sample_barcode_a = substr(tumor_barcode_a, 1, 15),
       sample_barcode_b = substr(tumor_barcode_b, 1, 15)) %>% 
  left_join(surgeries, by=c("sample_barcode_a" = "sample_barcode")) %>% 
  left_join(surgeries, by=c("sample_barcode_b" = "sample_barcode")) %>% 
  mutate(surgery_pair = paste(surgery_number.x, surgery_number.y, sep="-")) %>% 
  filter(received_tmz == 1 | received_alkylating_agent== 1, surgery_pair%in%c("1-2", "2-3","3-4"))

# Divide into IDHwt and IDHmut noncodel
clinical_silver_IDHwt_seq = clinical_silver_seq %>% 
  filter(idh_codel_subtype == "IDHwt_noncodel") %>% 
  filter(received_tmz == 1 | received_alkylating_agent== 1)
clinical_silver_IDHmut_seq = clinical_silver_seq %>% 
  filter(idh_codel_subtype == "IDHmut_noncodel") %>% 
  filter(received_tmz == 1 | received_alkylating_agent== 1)

# IDHwt - PFS. n = 89 (non-hypermutants = 76, hypermutants = 13).
fit_IDHwt_seq_int <- survfit(Surv(surgical_interval, rep(1, length(clinical_silver_IDHwt_seq$surgical_interval))) ~ hypermutator_status,
                         data = clinical_silver_IDHwt_seq)
ggsurvplot(fit_IDHwt_seq_int, data = clinical_silver_IDHwt_seq, risk.table = FALSE, pval= TRUE, pval.method = FALSE, pval.coord = c(40, 0.50), 
           surv.median.line = "v", palette = c("royalblue4", "tomato3"), 
           ylab = "Progression-free survival \n probability")
# To restrict to a specific interval/window use: break.time.by = 10, xlim = c(0,36).

# Univariate Cox model for IDHwt
hyper_cox_model_IDHwt = coxph(Surv(surgical_interval, rep(1, length(clinical_silver_IDHwt_seq$surgical_interval))) ~ hypermutator_status, data = clinical_silver_IDHwt_seq)
summary(hyper_cox_model_IDHwt)


# IDHmut noncodels. n = 26 (n = 14, n =12).
fit_IDHmut_seq_int <- survfit(Surv(surgical_interval, rep(1, length(clinical_silver_IDHmut_seq$surgical_interval))) ~ hypermutator_status,
                             data = clinical_silver_IDHmut_seq)
ggsurvplot(fit_IDHmut_seq_int, data = clinical_silver_IDHmut_seq, risk.table = FALSE, pval= TRUE, pval.method = FALSE, pval.coord = c(40, 0.75),
           surv.median.line = "v", palette = c("royalblue4", "tomato3"), 
           ylab = "Progression-free survival \n probability")

# Univariate Cox model for IDHmut:
hyper_cox_model_IDHmut = coxph(Surv(surgical_interval, rep(1, length(clinical_silver_IDHmut_seq$surgical_interval))) ~ hypermutator_status, data = clinical_silver_IDHmut_seq)
summary(hyper_cox_model_IDHmut)



