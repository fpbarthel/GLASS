##################################################
# Visualize aneuploidy results based on GATK CNV
# Updated: 2019.01.16
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

##################################################
# Establish connection with database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

# Load in the blocklist to implement CN thresholds:
tumor_blocklist = dbReadTable(con,  Id(schema="analysis", table="blocklist"))
pairs = dbReadTable(con,  Id(schema="analysis",table="pairs"))
tumor_pairs = dbReadTable(con,  Id(schema="analysis", table="tumor_pairs"))
aliquots = dbReadTable(con,  Id(schema="biospecimen", table="aliquots"))
surgeries = dbReadTable(con,  Id(schema="clinical", table="surgeries"))
cases = dbReadTable(con,  Id(schema="clinical", table="cases"))
subtypes = dbReadTable(con,  Id(schema="clinical", table="subtypes"))
mut_freq = dbReadTable(con,  Id(schema="analysis",table="mutation_freq"))
silver_set = dbGetQuery(con, "SELECT * FROM analysis.silver_set")
clinal_tumor_pairs = dbGetQuery(con,"SELECT * FROM analysis.clinical_by_tumor_pair")
# Use TITAN estimates of whole genome doubling (WGD) with ploidy values.
titan_param = dbReadTable(con,  Id(schema="analysis",table="titan_params"))


# Floris has created two aneuploidy metrics 
# 1. to represent the fraction of the genome with copy number alterations 
# 2. An aneuploidy score to mimic that of Taylor et al. Cancer Cell 2018.

# Query that should return both aneuploidy estimates.
aneuploidy_final = dbGetQuery(con, "SELECT ta.aliquot_barcode,aneuploidy_score,aneuploidy,idh_codel_subtype,surgery_number
FROM analysis.taylor_aneuploidy ta
                              LEFT JOIN analysis.aneuploidy aa ON aa.aliquot_barcode = ta.aliquot_barcode
                              LEFT JOIN analysis.blocklist bl ON bl.aliquot_barcode = ta.aliquot_barcode
                              LEFT JOIN biospecimen.aliquots al ON al.aliquot_barcode = ta.aliquot_barcode
                              LEFT JOIN clinical.surgeries su ON su.sample_barcode = al.sample_barcode
                              WHERE bl.coverage_exclusion = 'allow' AND bl.cnv_exclusion <> 'block' AND bl.fingerprint_exclusion = 'allow'")

# Add metadata to merge with other public datasets.
aneuploidy_final = aneuploidy_final %>% 
  mutate_if(bit64::is.integer64, as.double) %>% 
  mutate(sample_barcode = substr(aliquot_barcode, 1, 15),
         seq_type = substr(aliquot_barcode, 21, 23))

# Examine the relationship between the two metrics for all samples.
cor.test(aneuploidy_final$aneuploidy_score, aneuploidy_final$aneuploidy, method = "spearman")
# Report p-value and correlation coefficient.
ggplot(aneuploidy_final, aes(x = aneuploidy, y = aneuploidy_score)) + geom_point() + theme_bw() +ggtitle("rho = 0.73, P-val < 2.2E-16 (n=522)") +
  xlab("GLASS - aneuploidy value") + ylab("GLASS - aneuploidy score")

# Taylor et al. provided aneuploidy scores as part of their paper. Extract available LGG and GBM samples that overlap
# with GLASS. Note: The aneuploidy data in the paper was generated using ABSOLUTE for array-based data in TCGA.
taylor_data = readWorkbook("/Users/johnsk/Documents/Life-History/titan-analyses/data/taylor-cancer-cell-2018.xlsx", sheet = 1, rowNames = F, colNames = TRUE)

# Extract only the iniformation for GBM and LGG. Recode some of their metadata to fit our own.
glioma_dat = taylor_data %>% 
  filter(Type%in%c("GBM", "LGG")) %>% 
  mutate(sample_type_num = substr(Sample, 14, 15), 
         sample_type = recode(sample_type_num, "01" = "TP", "02" = "R1"),
         sample_barcode = paste(substr(Sample, 1, 12), sample_type, sep ="-"))

# Examine how the AneuploidyScore data determined by ABSOLUTE match with our estimates of aneuploidy.
glass_taylor_aneuploidy = aneuploidy_final %>%
  inner_join(glioma_dat, by="sample_barcode")

# Create plot of two aneuploidy scores vs. one another. Make sure both metrics consider all the same chromosomal arms.
cor.test(glass_taylor_aneuploidy$`AneuploidyScore(AS)`, glass_taylor_aneuploidy$aneuploidy_score, method = "spearman")
ggplot(glass_taylor_aneuploidy, aes(x = `AneuploidyScore(AS)`, y = aneuploidy_score, color = seq_type)) + geom_point() + theme_bw() +
  ylab("GLASS - aneuploidy score") + xlab("Taylor et al. AneuploidyScore(AS)") + ggtitle("rho = 0.81, P-val < 2.2E-16 (n=81)") + labs(color = "Seq platform") + xlim(0, 39) + ylim(0,39)

# Examine the relationship between the fraction of the genome altered aneuploidy score and AneuploidyScores.
cor.test(glass_taylor_aneuploidy$`AneuploidyScore(AS)`, glass_taylor_aneuploidy$aneuploidy, method = "spearman")
ggplot(glass_taylor_aneuploidy, aes(x = `AneuploidyScore(AS)`, y = aneuploidy, color = seq_type)) + geom_point() + theme_bw() +
  ylab("GLASS - aneuploidy value") + xlab("Taylor et al. AneuploidyScore(AS)") + ggtitle("rho = 0.77, P-val < 2.2E-16 (n=81)") + labs(color = "Seq platform") + xlim(0, 39) + ylim(0,1)

# Use this information in the context of the silver set for aneuploidy pairs. Note: that this includes some poor CN data.
aneuploidy_pairs = dbGetQuery(con, "SELECT tumor_pair_barcode, tumor_barcode_a, tumor_barcode_b, a1.aneuploidy AS aneuploidy_a, a2.aneuploidy AS aneuploidy_b, t1.aneuploidy_score::integer AS aneuploidy_score_a,t2.aneuploidy_score::integer AS aneuploidy_score_b,idh_codel_subtype
                              FROM analysis.silver_set ss
                              LEFT JOIN analysis.aneuploidy a1 ON a1.aliquot_barcode = ss.tumor_barcode_a
                              LEFT JOIN analysis.aneuploidy a2 ON a2.aliquot_barcode = ss.tumor_barcode_b
                              LEFT JOIN analysis.taylor_aneuploidy t1 ON t1.aliquot_barcode = ss.tumor_barcode_a
                              LEFT JOIN analysis.taylor_aneuploidy t2 ON t2.aliquot_barcode = ss.tumor_barcode_b
                              LEFT JOIN biospecimen.aliquots al ON al.aliquot_barcode = ss.tumor_barcode_a
                              LEFT JOIN clinical.surgeries su ON su.sample_barcode = al.sample_barcode")


# Filter out any poor performing samples on CNV blocklist. 197 patients because I merged silver set with allow|review CN.
# Gold set has two tumor pairs (different combination of tumors) that are not in silver set.
aneuploidy_pairs_filtered = aneuploidy_pairs %>% 
  inner_join(tumor_blocklist, by=c("tumor_barcode_a"="aliquot_barcode")) %>% 
  inner_join(tumor_blocklist, by=c("tumor_barcode_b"="aliquot_barcode")) %>% 
  # Note db introduced some trailing whitespace.
  filter(cnv_exclusion.x %in%c("allow ", "review")) %>% 
  filter(cnv_exclusion.y %in% c("allow ", "review")) %>% 
  mutate(aneuploidy_diff = aneuploidy_b - aneuploidy_a) %>% 
  left_join(clinal_tumor_pairs, by="tumor_pair_barcode") %>%
  # Relabel subtypes for aesthetics.
  mutate(idh_codel_subtype = recode(idh_codel_subtype, "IDHmut_codel" = "IDHmut codel", "IDHmut_noncodel" = "IDHmut noncodel", "IDHwt_noncodel" = "IDHwt"))


# Define aneuploidy pairs by IDH subtype.
aneuploidy_pairs_IDHwt = aneuploidy_pairs_filtered %>% 
  filter(idh_codel_subtype == "IDHwt")
aneuploidy_pairs_IDH_codel = aneuploidy_pairs_filtered %>% 
  filter(idh_codel_subtype == "IDHmut codel")
aneuploidy_pairs_IDH_noncodel = aneuploidy_pairs_filtered %>% 
  filter(idh_codel_subtype == "IDHmut noncodel")

# Interesting, the fraction of the genome altered approach is significantly different at recurrence, but aneuploidy score is not. 
wilcox.test(aneuploidy_pairs_filtered$aneuploidy_a, aneuploidy_pairs_filtered$aneuploidy_b, paired = TRUE)
wilcox.test(aneuploidy_pairs_filtered$aneuploidy_score_a, aneuploidy_pairs_filtered$aneuploidy_score_b, paired = TRUE)

# Subtype specific analyses that now have filtered samples.
wilcox.test(aneuploidy_pairs_IDHwt$aneuploidy_a, aneuploidy_pairs_IDHwt$aneuploidy_b, paired = TRUE)
wilcox.test(aneuploidy_pairs_IDH_codel$aneuploidy_a, aneuploidy_pairs_IDH_codel$aneuploidy_b, paired = TRUE)
wilcox.test(aneuploidy_pairs_IDH_noncodel$aneuploidy_a, aneuploidy_pairs_IDH_noncodel$aneuploidy_b, paired = TRUE)

# When comparing aneuploidy scores between the primary and recurrent, there are significant difference for both IDHmut groups.
wilcox.test(aneuploidy_pairs_IDHwt$aneuploidy_score_a, aneuploidy_pairs_IDHwt$aneuploidy_score_b, paired = TRUE)
wilcox.test(aneuploidy_pairs_IDH_codel$aneuploidy_score_a, aneuploidy_pairs_IDH_codel$aneuploidy_score_b, paired = TRUE)
wilcox.test(aneuploidy_pairs_IDH_noncodel$aneuploidy_score_a, aneuploidy_pairs_IDH_noncodel$aneuploidy_score_b, paired = TRUE)

# Combine the aneuploidy values for each tumor pair.
aneuploid_plot_value = aneuploidy_pairs_filtered %>% 
  gather(sample, aneuploidy, c(aneuploidy_a, aneuploidy_b)) 
ggplot(aneuploid_plot_value, aes(x = sample, y = aneuploidy, group = tumor_pair_barcode)) +
  geom_line(linetype="solid", size=1, alpha = 0.3) + ylab("GLASS - aneuploidy value") + xlab("Gold set pairs (n=197)") +
  geom_point(color="black", size=2) + theme_bw() + facet_grid(~idh_codel_subtype)

# Combine the aneuploidy scores for each tumor pair.
aneuploid_plot_score = aneuploidy_pairs_filtered %>% 
  gather(sample, aneuploidy_score, c(aneuploidy_score_a, aneuploidy_score_b)) 
ggplot(aneuploid_plot_score, aes(x = sample, y = aneuploidy_score, group = tumor_pair_barcode)) +
  geom_line(linetype="solid", size=1, alpha = 0.3) + ylab("GLASS - aneuploidy score") + xlab("Gold set pairs (n=197)") +
  geom_point(color="black", size=2) + theme_bw() + facet_grid(~idh_codel_subtype)

# Inspect TITAN ploidy distribution.
ggplot(titan_param, aes(x=ploidy)) + geom_histogram() + theme_bw()
# It looks like there is a right tail after 2.5 ploidy.
titan_wgd_aneu = titan_param %>% 
  left_join(pairs, by="pair_barcode") %>% 
  mutate(wgd = ifelse(as.numeric(ploidy) >= 3, 1, 0)) %>% 
  inner_join(aneuploidy_final, by=c("tumor_barcode"="aliquot_barcode")) %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, "IDHmut_codel" = "IDHmut codel", "IDHmut_noncodel" = "IDHmut noncodel", "IDHwt_noncodel" = "IDHwt"))

# While we don't have actual whole genome doubling calls (like ABSOLUTE) we can test whether higher ploidy samples have more aneuploidy.
cor.test(titan_wgd_aneu$ploidy, titan_wgd_aneu$aneuploidy_score, method = "spearman")
ggplot(titan_wgd_aneu, aes(x=ploidy, y= aneuploidy_score, color=as.factor(wgd))) +
  geom_point() + theme_bw() + facet_grid(~idh_codel_subtype) +
  labs(color = "Ploidy > 3") + xlab("TITAN estimated ploidy") + ylab("GLASS - aneuploidy score") + ggtitle("N = 522")

### Associations between aneuploidy and survival. 
case_level_subtype = surgeries %>% 
  select(case_barcode, idh_codel_subtype) %>% 
  distinct() %>% 
  filter(!is.na(idh_codel_subtype)) 

# Test aneuploidy with post-recurrence survival.
aneuploidy_idh_noncodel = aneuploidy_pairs_IDH_noncodel %>% 
  left_join(cases, by="case_barcode") %>% 
  mutate(post_recurrence_surv = case_overall_survival_mo-surgical_interval,
         # summary(aneuploidy_idh_noncodel$aneuploidy_diff) = 0.12980, upper quartile vs. bottom three quartiles.
         # I tried mean as well as median for thresholds. Mean or upper quartile seem most appropriate based on distribution.
         acquired_aneuploidy = ifelse(aneuploidy_diff > 0.12980, 1, 0),
         patient_vital = ifelse(case_vital_status=="alive", 0, 1))

# NOTE: Make sure you are using the gold set (Or silver set filtered through for CN allow|review.
ggplot(aneuploidy_idh_noncodel, aes(aneuploidy_diff)) + geom_histogram() + theme_bw() + geom_vline(xintercept = 0.12980, color="red") +
  xlab("Aneuploidy value difference") + ggtitle("Gold set (n=65)")

# Overall survival (months) with the designation of "acquired_aneuploidy".
noncodel_aneuploidy_surv <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ acquired_aneuploidy,
                                    data = aneuploidy_idh_noncodel)
ggsurvplot(noncodel_aneuploidy_surv, data = aneuploidy_idh_noncodel, risk.table = TRUE, pval= TRUE) + ggtitle("Overall survival (high aneuploidy, Gold set)")

# Post recurrence survival (months) with the designation of "acquired_aneuploidy".
noncodel_aneuploidy_post_recur <- survfit(Surv(post_recurrence_surv, patient_vital) ~ acquired_aneuploidy,
                                          data = aneuploidy_idh_noncodel)
ggsurvplot(noncodel_aneuploidy_post_recur, data = aneuploidy_idh_noncodel, risk.table = TRUE, pval= TRUE) + ggtitle("Post recurrence survival (high aneuploidy, Gold set)")

# Overall survival (months) with the designation of "acquired_aneuploidy".
cox_fit_noncodel_1 = coxph(Surv(case_overall_survival_mo, patient_vital) ~ aneuploidy_diff, data= aneuploidy_idh_noncodel)
summary(cox_fit_noncodel_1)
# aneuploidy_diff treated as a continuous variable.
cox_fit_noncodel_2 = coxph(Surv(case_overall_survival_mo, patient_vital) ~ acquired_aneuploidy, data= aneuploidy_idh_noncodel)
summary(cox_fit_noncodel_2)

# Post recurrence survival (months) with the designation of "acquired_aneuploidy".
cox_fit_noncodel_3 = coxph(Surv(post_recurrence_surv, patient_vital) ~ aneuploidy_diff, data= aneuploidy_idh_noncodel)
summary(cox_fit_noncodel_3)
# aneuploidy_diff treated as a continuous variable.
cox_fit_noncodel_4 = coxph(Surv(post_recurrence_surv, patient_vital) ~ acquired_aneuploidy, data= aneuploidy_idh_noncodel)
summary(cox_fit_noncodel_4)

# Is there an association between patients that received TMZ and aneuploidy difference (across all three subtypes)?
wilcox.test(aneuploidy_pairs_IDHwt$aneuploidy_diff~aneuploidy_pairs_IDHwt$received_tmz)
wilcox.test(aneuploidy_pairs_IDH_codel$aneuploidy_diff~aneuploidy_pairs_IDH_codel$received_tmz)
wilcox.test(aneuploidy_pairs_IDH_noncodel$aneuploidy_diff~aneuploidy_pairs_IDH_noncodel$received_tmz)

# Is aneuploidy_difference between tumor_a and tumor_b related with radiation therapy?
wilcox.test(aneuploidy_pairs_IDHwt$aneuploidy_diff~aneuploidy_pairs_IDHwt$received_rt)
wilcox.test(aneuploidy_pairs_IDH_codel$aneuploidy_diff~aneuploidy_pairs_IDH_codel$received_rt)
wilcox.test(aneuploidy_pairs_IDH_noncodel$aneuploidy_diff~aneuploidy_pairs_IDH_noncodel$received_rt)

# Is aneuploidy_difference between tumor_a and tumor_b related with hypermutation?
wilcox.test(aneuploidy_pairs_IDHwt$aneuploidy_diff~aneuploidy_pairs_IDHwt$hypermutator_status)
wilcox.test(aneuploidy_pairs_IDH_codel$aneuploidy_diff~aneuploidy_pairs_IDH_codel$hypermutator_status)
wilcox.test(aneuploidy_pairs_IDH_noncodel$aneuploidy_diff~aneuploidy_pairs_IDH_noncodel$hypermutator_status)


### Arm-level ###
             
### Investigate the arm-level calls, compare with Taylor et al.
arm_level_calls = dbGetQuery(con, "SELECT * FROM analysis.cnv_by_arm_gatk")

# Double check to make sure the arm calls don't have any normals and are passable copy number estimates.
arm_level_calls_filtered = arm_level_calls %>% 
  filter(!grepl("-NB-", aliquot_barcode)) %>% 
  inner_join(tumor_blocklist, by="aliquot_barcode") %>% 
  # Note remove trailing whitespace.
  filter(cnv_exclusion %in%c("allow ", "review"))

# Creating summary variables for each sample. AneuploidyScore for amplifications and deletions.
arm_level_summary_1 =  arm_level_calls_filtered %>% 
  group_by(aliquot_barcode, arm_call) %>% 
  summarise(events = n()) %>% 
  filter(arm_call != 0) %>% 
  spread(arm_call, events) %>% 
  mutate(sample_barcode = substr(aliquot_barcode, 1, 15)) %>% 
  inner_join(surgeries, by="sample_barcode") %>% 
  select(aliquot_barcode, sample_barcode, surgery_number, idh_codel_subtype, AS_amp = `1`, AS_del = `-1`) 

# Reduce set to only those arm-level analyses considered by Taylor et al.
arm_level_summary2 = arm_level_calls_filtered %>% 
  select(-chrom) %>% 
  spread(arm, arm_call) %>% 
  select(aliquot_barcode, `1p`, `1q`, `2p`, `2q`, `3p`, `3q`, `4p`, `4q`, `5p`, `5q`, `6p`, `6q`, `7p`, `7q`, `8p`,
         `8q`, `9p`, `9q`, `10p`,`10q`, `11p`, `11q`, `12p`, `12q`, `13q`, `14q`, `15q`, `16q`, `16p`, `17p`, `17q`,
         `18p`, `18q`, `19p`, `19q`, `20p`,`20q`, `21q`, `22q`) 

# The CN values should be filtered so merge with the silver set.
silver_set_merge = silver_set %>% 
  gather(sample_type, tumor_barcode, c(tumor_barcode_a, tumor_barcode_b))

# Combine two tables.
arm_level_full = arm_level_summary_1 %>% 
  inner_join(arm_level_summary2, by="aliquot_barcode") %>% 
  inner_join(silver_set_merge, by=c("aliquot_barcode"="tumor_barcode")) %>% 
  select(aliquot_barcode:`22q`) %>% 
  mutate(AS_amp = ifelse(is.na(AS_amp), 0, AS_amp),
         AS_del = ifelse(is.na(AS_del), 0, AS_del))
  
# Write out table to be used in downstream analyses or uploaded as a table to the database.
# write.table(arm_level_full, file = "/Users/johnsk/Documents/glass-arm-level-20190116.txt", sep="\t", row.names = F, col.names = T, quote = F)

# Compare the arm level calls with those from **Taylor et al. Cancer Cell 2018**.
arm_level_merge = arm_level_full %>% 
  mutate(dataset = "glass") %>% 
  select(aliquot_barcode, sample_barcode, dataset, AS_del, AS_amp, `1p`:`22q`) %>% 
  filter(sample_barcode%in%glioma_dat$sample_barcode) 
# Prepare the Taylor data to be merged with GLASS data.
taylor_data_merged = glioma_dat %>% 
  mutate(dataset = "taylor") %>% 
  select(aliquot_barcode = Sample, sample_barcode, dataset, AS_del, AS_amp, `1p`:`22q`) %>% 
  filter(sample_barcode%in%arm_level_merge$sample_barcode) %>% 
  bind_rows(arm_level_merge)
# After writing out table and sorting by sample_barcode, taylor and glass arm-level calls were comparable even for NAs.
# write.table(taylor_data_merged, file = "/glass-taylor-arm-level-20190116.txt", sep="\t", row.names = F, col.names = T, quote = F)

# Finally, analyze the Aneuploidy score at the arm-level.
taylor_glass = glioma_dat %>% 
  mutate(dataset = "taylor") %>% 
  select(aliquot_barcode = Sample, sample_barcode, dataset, AS_del, AS_amp, `1p`:`22q`) %>% 
  filter(sample_barcode%in%arm_level_merge$sample_barcode) %>% 
  inner_join(arm_level_merge, by="sample_barcode")

# Plot the AS_amp scores vs. one another.
cor.test(taylor_glass$AS_amp.x, taylor_glass$AS_amp.y, method = "spearman")
ggplot(taylor_glass, aes(x = AS_amp.x, y = AS_amp.y)) + geom_point() + theme_bw() +
  xlab("Taylor - AS score amplification") + ylab("GLASS - AS score amplification") + ggtitle("rho = 0.67, P-val < 1.0E-04 (n=37)") + xlim(0, 39) + ylim(0,39)
# Plot the AS_del scores vs. one another.
cor.test(taylor_glass$AS_del.x, taylor_glass$AS_del.y, method = "spearman")
ggplot(taylor_glass, aes(x = AS_del.x, y = AS_del.y)) + geom_point() + theme_bw() +
  xlab("Taylor - AS score deletion") + ylab("GLASS - AS score deletion") + ggtitle("rho = 0.91, P-val < 2.1E-15 (n=37)") + xlim(0, 39) + ylim(0, 39)