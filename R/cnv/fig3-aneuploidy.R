##################################################
# Determine changes in aneuploidy results based on GATK CNV
# Updated: 2019.05.15
# Author: Kevin J.
##################################################

# Working directory for this analysis in the GLASS-analysis project. 
mybasedir = "/Volumes/verhaak-lab/GLASS-analysis/"
setwd(mybasedir)

##################################################

# Necessary packages:
library(tidyverse)
library(DBI)
library(openxlsx)
library(survminer)
library(survival)
library(pheatmap)
library(RColorBrewer)
library(cowplot)

##################################################
# Establish connection with database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB2")

# Load in the necessary data tables as they relate to aneuploidy.
mutation_freq = dbGetQuery(con, "SELECT * FROM analysis.mut_freq")
gold_set = dbReadTable(con,  Id(schema="analysis", table="gold_set"))
gatk_aneuploidy = dbReadTable(con,  Id(schema="analysis", table="gatk_aneuploidy"))
pairs = dbReadTable(con,  Id(schema="analysis", table="pairs"))
blocklist = dbReadTable(con,  Id(schema="analysis", table="blocklist"))
subtypes = dbReadTable(con,  Id(schema="clinical", table="subtypes"))
cases = dbReadTable(con,  Id(schema="clinical", table="cases"))
seqz_params = dbReadTable(con,  Id(schema="variants", table="seqz_params"))
titan_params = dbReadTable(con,  Id(schema="variants", table="titan_params"))

# Construct a table that provides subject-level information about clinical variables between two timepoints.
clinical_tumor_pairs_query = read_file("/Users/johnsk/Documents/Life-History/glass-analyses/sql/clinical-tumor-pairs-db2.sql")
clinical_tumor_pairs <- dbGetQuery(con, clinical_tumor_pairs_query)

# Define aneuploidy pairs from the silver set.
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

# Define aneuploidy pairs by IDH subtype.
aneuploidy_pairs_IDHwt = aneuploidy_pairs_clin %>% 
  filter(idh_codel_subtype == "IDHwt")
aneuploidy_pairs_IDH_codel = aneuploidy_pairs_clin %>% 
  filter(idh_codel_subtype == "IDHmut-codel")
aneuploidy_pairs_IDH_noncodel = aneuploidy_pairs_clin %>% 
  filter(idh_codel_subtype == "IDHmut-noncodel")

# Fraction of the genome with copy number alterations.
wilcox.test(aneuploidy_pairs_IDH_codel$aneuploidy_a, aneuploidy_pairs_IDH_codel$aneuploidy_b, paired = TRUE)

# High level of significance between primary and recurrence for IDHmut noncodels.
wilcox.test(aneuploidy_pairs_IDH_noncodel$aneuploidy_a, aneuploidy_pairs_IDH_noncodel$aneuploidy_b, paired = TRUE)
wilcox.test(aneuploidy_pairs_IDHwt$aneuploidy_a, aneuploidy_pairs_IDHwt$aneuploidy_b, paired = TRUE)

# When comparing aneuploidy SCOREs (i.e., arm-level events) between the primary and recurrent, there are significant difference for both IDHmut groups.
wilcox.test(aneuploidy_pairs_IDH_codel$aneuploidy_score_a, aneuploidy_pairs_IDH_codel$aneuploidy_score_b, paired = TRUE)
wilcox.test(aneuploidy_pairs_IDH_noncodel$aneuploidy_score_a, aneuploidy_pairs_IDH_noncodel$aneuploidy_score_b, paired = TRUE)
wilcox.test(aneuploidy_pairs_IDHwt$aneuploidy_score_a, aneuploidy_pairs_IDHwt$aneuploidy_score_b, paired = TRUE)

# Arm-level gains. Statistically significant change in arm-level events.
wilcox.test(aneuploidy_pairs_IDHwt$aneuploidy_amp_score_a, aneuploidy_pairs_IDHwt$aneuploidy_amp_score_b, paired = TRUE)
# Arm-level losses. No changes.
wilcox.test(aneuploidy_pairs_IDHwt$aneuploidy_del_score_a, aneuploidy_pairs_IDHwt$aneuploidy_del_score_b, paired = TRUE)

# Combine the aneuploidy values for each tumor pair.
aneuploid_plot_value = aneuploidy_pairs_clin %>% 
  gather(sample, aneuploidy, c(aneuploidy_a, aneuploidy_b)) %>% 
  mutate(sample = recode(sample,  "aneuploidy_a" = "Primary",  "aneuploidy_b" = "Recurrence"))

# Ladder plot for paired gold set samples.
pdf(file = "/Users/johnsk/Documents/aneuploidy-value-diff-goldset.pdf", height = 4, width = 8, bg = "transparent", useDingbats = FALSE)
ggplot(aneuploid_plot_value, aes(x = sample, y = aneuploidy, group = tumor_pair_barcode, color = sample)) +
  geom_line(linetype="solid", size=1, alpha = 0.3, color ="black") + ylab("Aneuploidy value") + xlab("") + theme_bw() +
  geom_point(size=2) + scale_color_manual(values = c("Primary" = "#a6611a", "Recurrence" = "#018571"), name = "Sample Type") + facet_grid(~idh_codel_subtype)
dev.off()

# Combine the aneuploidy scores for each tumor pair.
aneuploid_plot_score = aneuploidy_pairs_clin %>% 
  gather(sample, aneuploidy_score, c(aneuploidy_score_a, aneuploidy_score_b)) %>% 
  mutate(sample = recode(sample,  "aneuploidy_score_a" = "Primary",  "aneuploidy_score_b" = "Recurrence"))

# Ladder plot for aneuploidy score in gold set samples.
pdf(file = "/Users/johnsk/Documents/aneuploidy-scores-diff-goldset.pdf", height = 4, width = 8, bg = "transparent", useDingbats = FALSE)
ggplot(aneuploid_plot_score, aes(x = sample, y = aneuploidy_score, group = tumor_pair_barcode, color=sample)) +
  geom_line(linetype="solid", size=1, alpha = 0.3, color ="black") + ylab("Aneuploidy score") + xlab("") + theme_bw() +
  geom_point(size=2) + scale_color_manual(values = c("Primary" = "#a6611a", "Recurrence" = "#018571"), name = "Sample Type") + facet_grid(~idh_codel_subtype)
dev.off()

# We know that IDHmut-noncodel samples have a difference in aneuploidy between the two timepoints.
# How is the difference distributed?
aneuploidy_pairs_IDH_noncodel = aneuploidy_pairs_IDH_noncodel %>% 
  mutate(aneuploidy_diff = aneuploidy_b-aneuploidy_a)
p1 = ggplot(aneuploidy_pairs_IDH_noncodel, aes(x = 1, y = aneuploidy_diff)) + geom_boxplot() + ylab("delta aneuploidy") + xlab("") +
  coord_flip() +  ylim(min(aneuploidy_pairs_IDH_noncodel$aneuploidy_diff), 0.75) + theme(axis.text.y=element_blank())
p2 = ggplot(aneuploidy_pairs_IDH_noncodel, aes(aneuploidy_diff)) + geom_histogram() + xlab("delta aneuploidy") +
  xlim(min(aneuploidy_pairs_IDH_noncodel$aneuploidy_diff), 0.75) + geom_vline(xintercept = 0.12, color="red")
pdf(file = "/Users/johnsk/Documents/idhmut-noncodel-aneuploidy-difference.pdf", height = 4, width = 8, bg = "transparent", useDingbats = FALSE)
plot_grid(p1, p2, labels = c("", ""), align = "v", ncol = 1)
dev.off()

# What is the upper quartile of acquired aneuploidy?
summary(aneuploidy_pairs_IDH_noncodel$aneuploidy_diff)

# Test aneuploidy with post-recurrence survival.
aneuploidy_idh_noncodel = aneuploidy_pairs_IDH_noncodel %>% 
  mutate(post_recurrence_surv = case_overall_survival_mo-surgical_interval,
         acquired_aneuploidy = ifelse(aneuploidy_diff >= 0.12285, 1, 0),
         aneuploidy_levels = ntile(aneuploidy_diff, 3),
         patient_vital = ifelse(case_vital_status=="alive", 0, 1)) %>% 
  left_join(mutation_freq, by=c("tumor_barcode_b" = "aliquot_barcode")) %>% 
  mutate(hypermutation = ifelse(coverage_adj_mut_freq > 10, 1, 0)) %>% 
  dplyr::select(tumor_pair_barcode, tumor_barcode_a, tumor_barcode_b, aneuploidy_a, aneuploidy_b, aneuploidy_diff, acquired_aneuploidy, hypermutation, mutation_count, patient_vital, case_overall_survival_mo, post_recurrence_surv) 

# Are hypermutant tumor also the acquired_aneuploidy tumors?
table(aneuploidy_idh_noncodel$hypermutation)

# Aneuploidy and hypermutation do not seem to be associated.
fisher.test(table(aneuploidy_idh_noncodel$acquired_aneuploidy, aneuploidy_idh_noncodel$hypermutation))

# Overall survival (months) with the designation of "acquired_aneuploidy".
noncodel_aneuploidy_surv <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ factor(acquired_aneuploidy),
                                    data = aneuploidy_idh_noncodel)
pdf(file = "/Users/johnsk/Documents/f3c-kcj.pdf", height = 5, width = 7, bg = "transparent", useDingbats = FALSE)
ggsurvplot(noncodel_aneuploidy_surv, data = aneuploidy_idh_noncodel, risk.table = FALSE, pval= TRUE, pval.method = FALSE, pval.coord = c(150, 0.75), 
           surv.median.line = "v", ylab = "Overall survival \n probability", xlab = "Time (months)", palette = c("#27408B", "#CD4F39")) 
dev.off()

# Does this change when age is included as a variable? 
cases_age = cases %>% select(case_barcode, case_age_diagnosis_years)
aneuploidy_idh_noncodel = aneuploidy_idh_noncodel %>% 
  mutate(case_barcode = substr(tumor_pair_barcode, 1, 12)) %>% 
  inner_join(cases_age, by="case_barcode")
# Treated as a factor or continuous?
aneuploidy_cox_model = coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + factor(acquired_aneuploidy), data = aneuploidy_idh_noncodel)
summary(aneuploidy_cox_model) # Yes, P = 0.004.
aneuploidy_cox_model = coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + aneuploidy_diff, data = aneuploidy_idh_noncodel)
summary(aneuploidy_cox_model) # No, P = 0.32.

#################################
# Tangential analysis: inspect TITAN estimated purity and ploidy differences between two timepoints.
#################################
# Aside: Inspect TITAN ploidy distribution.
titan_aliquot = titan_params %>% 
  inner_join(pairs, by="pair_barcode") %>% 
  left_join(blocklist, by="aliquot_barcode") %>% 
  left_join(blocklist, by="aliquot_barcode") %>% 
  filter(cnv_exclusion.x %in%c("allow ", "review"), cnv_exclusion.y %in% c("allow ", "review")) %>% 
  dplyr::select(aliquot_barcode = tumor_barcode, purity, ploidy)

# Define silver_set for TITAN with purity and ploidy metrics.
gold_titan = gold_set %>% 
  left_join(titan_aliquot, by=c("tumor_barcode_a"="aliquot_barcode")) %>% 
  left_join(titan_aliquot, by=c("tumor_barcode_b"="aliquot_barcode")) %>% 
  left_join(subtypes, by="case_barcode") %>% 
  select(tumor_pair_barcode:tumor_barcode_b, idh_codel_subtype, purity_a = purity.x, ploidy_a = ploidy.x, purity_b = purity.y, ploidy_b = ploidy.y)

# Define TITAN silver set by IDH subtype.
titan_pairs_IDHwt = gold_titan %>% 
  filter(idh_codel_subtype == "IDHwt")
titan_pairs_IDH_codel = gold_titan %>% 
  filter(idh_codel_subtype == "IDHmut-codel")
titan_pairs_IDH_noncodel = gold_titan %>% 
  filter(idh_codel_subtype == "IDHmut-noncodel")

# Plot both TITAN purity and ploidy for silver set.
titan_plot_purity = gold_titan %>% 
  gather(sample, purity, c(purity_a, purity_b)) %>% 
  mutate(sample = recode(sample,  "purity_a" = "Primary",  "purity_b" = "Recurrence"))
# How many samples have all the data? n = 220.
sum(complete.cases(gold_titan))
pdf(file = "/Users/johnsk/Documents/TITAN_purity_goldset_n220.pdf", height = 4, width = 8, bg = "transparent", useDingbats = FALSE)
ggplot(titan_plot_purity, aes(x=idh_codel_subtype, y=purity, fill=sample)) + 
  geom_boxplot()  + ylab("TITAN Purity") + xlab("") + theme_bw() + scale_fill_manual(values = c("Primary" = "#a6611a", "Recurrence" = "#018571"))
dev.off()

# Ploidy
titan_plot_ploidy = gold_titan %>% 
  gather(sample, ploidy, c(ploidy_a, ploidy_b)) %>% 
  mutate(sample = recode(sample,  "ploidy_a" = "Primary",  "ploidy_b" = "Recurrence"))
pdf(file = "/Users/johnsk/Documents/TITAN_ploidy_goldset_n220.pdf", height = 4, width = 8, bg = "transparent", useDingbats = FALSE)
ggplot(titan_plot_ploidy, aes(x=idh_codel_subtype, y=ploidy, fill=sample)) + 
  geom_boxplot()  + ylab("TITAN Ploidy") + xlab("") + theme_bw() + scale_fill_manual(values = c("Primary" = "#CA6C18", "Recurrence" = "#005496"))
dev.off()

# Test purity for statistical difference.
wilcox.test(titan_pairs_IDH_codel$purity_a, titan_pairs_IDH_codel$purity_b, paired = TRUE)
wilcox.test(titan_pairs_IDH_noncodel$purity_a, titan_pairs_IDH_noncodel$purity_b, paired = TRUE)
wilcox.test(titan_pairs_IDHwt$purity_a, titan_pairs_IDHwt$purity_b, paired = TRUE)
# Testing ploidy for statistical difference.
wilcox.test(titan_pairs_IDHwt$ploidy_a, titan_pairs_IDHwt$ploidy_b, paired = TRUE)
wilcox.test(titan_pairs_IDH_codel$ploidy_a, titan_pairs_IDH_codel$ploidy_b, paired = TRUE)
wilcox.test(titan_pairs_IDH_noncodel$ploidy_a, titan_pairs_IDH_noncodel$ploidy_b, paired = TRUE)


# It looks like there is a right tail after 2.5 ploidy.
titan_wgd_aneu = gold_titan %>% 
  mutate(wgd_a = ifelse(as.numeric(ploidy_a) >= 3, 1, 0),
         wgd_b = ifelse(as.numeric(ploidy_b) >= 3, 1, 0)) 
table(titan_wgd_aneu$wgd_a, titan_wgd_aneu$wgd_b)

##################################################################################
## Temporary code below this line:
## Examine which samples are different between the submission and resubmission sets:
## The decisions we made, appear to correct. Samples that should be blocked are now blocked.

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

# Load in the blocklist to implement CN thresholds:
tumor_blocklist = dbReadTable(con,  Id(schema="analysis", table="blocklist"))

tmp = blocklist %>% 
  filter(aliquot_barcode %in%c(silver_set$tumor_barcode_a, silver_set$tumor_barcode_b)) %>% 
  mutate(case_barcode = substr(aliquot_barcode, 1, 12)) %>% 
  left_join(subtypes, by="case_barcode") %>% 
  filter(idh_codel_subtype == "IDHmut-noncodel", !grepl("-NB-", aliquot_barcode)) %>% 
  inner_join(tumor_blocklist, by="aliquot_barcode") %>% 
  dplyr::select(aliquot_barcode, cnv_exclusion.x, cnv_exclusion_reason.x, cnv_exclusion.y, cnv_exclusion_reason.y)

# Manually review purity estimates from TITAN and Sequenza.
tmp2 = tmp[tmp$cnv_exclusion.x!=tmp$cnv_exclusion.y,]




