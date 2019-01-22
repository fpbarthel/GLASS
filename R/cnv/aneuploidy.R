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
library(pheatmap)
library(RColorBrewer)

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
cnv_by_gene_gatk = dbReadTable(con,  Id(schema="analysis",table="cnv_by_gene_gatk"))
# cnv_by_gene_gatk = read_csv("/Users/johnsk/Downloads/analysis.cnv_by_gene_gatk_2019-01-19T13_11_45+0000.csv")
pathway_gene_list = dbReadTable(con,  Id(schema="ref",table="driver_genes"))

# Use TITAN estimates of whole genome doubling (WGD) with ploidy values.
titan_param = dbReadTable(con,  Id(schema="analysis",table="titan_params"))

# Quickly get the frequency of each gene-level event in the silver set of recurrence.
cnv_example = cnv_by_gene_gatk %>% 
  filter(!grepl("-NB-", aliquot_barcode)) %>% 
  mutate(sample_barcode =  substr(aliquot_barcode, 1, 15)) %>% 
  inner_join(surgeries, by = "sample_barcode") %>% 
  inner_join(silver_set, by=c("aliquot_barcode"="tumor_barcode_b")) %>% 
  group_by(gene_symbol, hlvl_call, idh_codel_subtype) %>% 
  summarise(alt_breakdown = n()) 


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

# Test whether the aneuploidy change in IDHmut noncodels is associated with a gain in cell cycle genes.
cell_cycle_cnv = pathway_gene_list %>% 
  filter(pathway=="Cell cycle", has_cnv==1)
cnv_gatk_cell_cycle = cnv_by_gene_gatk %>% 
  inner_join(cell_cycle_cnv, by="gene_symbol") %>% 
  select(-pathway, -has_cnv, -has_mut, -wcr) %>% 
  spread(gene_symbol, hlvl_call) %>% 
  filter(!grepl("-NB-", aliquot_barcode))

# Need to include mutation information.
genedata <- dbGetQuery(con, read_file("/Users/johnsk/Documents/Life-History/glass-analyses/scripts/heatmap-mutation.sql"))
cell_cycle_muts = genedata %>% 
  filter(gene_symbol %in% c("RB1", "TP53"), variant_call=="R") %>% 
  select(case_barcode, gene_symbol, variant_classification) %>% 
  spread(gene_symbol, variant_classification) %>% 
  select(case_barcode, RB1_mut = RB1, TP53_mut = TP53)

cell_cycle_df_noncodel = aneuploidy_idh_noncodel %>% 
  select(-tumor_barcode_a.y, -tumor_barcode_b.y) %>% 
  inner_join(cnv_gatk_cell_cycle, by=c("tumor_barcode_a.x"="aliquot_barcode")) %>% 
  inner_join(cnv_gatk_cell_cycle, by=c("tumor_barcode_b.x"="aliquot_barcode")) %>% 
  arrange(idh_codel_subtype, tumor_pair_barcode) %>% 
  left_join(cell_cycle_muts, by="case_barcode")
colnames(cell_cycle_df_noncodel) = gsub("\\.x", "_P", colnames(cell_cycle_df_noncodel))
colnames(cell_cycle_df_noncodel) = gsub("\\.y", "_R", colnames(cell_cycle_df_noncodel)) 


cell_cycle_df_noncodel = cell_cycle_df_noncodel %>% 
  mutate(CCND2 = ifelse(CCND2_R==2 & CCND2_P!=2, "+CCND2", NA),
         CDK4 = ifelse(CDK4_R==2 & CDK4_P!=2, "+CDK4", NA),
         CDK6 = ifelse(CDK6_R==2 & CDK6_P!=2, "+CDK6", NA),
         CDKN2A = ifelse(CDKN2A_R==-2 & CDKN2A_P!=-2, "-CDKN2A", NA),
         MDM2 = ifelse(MDM2_R==2 & MDM2_P!=2, "+MDM2", NA),
         MDM4 = ifelse(MDM4_R==2 & MDM4_P!=2, "+MDM4", NA),
         RB1 = ifelse(RB1_R==-2 & RB1_P!=-2, "-RB1", NA),
         TP53 = ifelse(TP53_R==-2 & TP53_P!=-2, "-TP53", NA),
         cell_cycle = ifelse(is.na(CCND2) & is.na(CDK4) & is.na(CDK6) & is.na(CDKN2A) & is.na(MDM2) & is.na(MDM4) & is.na(RB1) & is.na(TP53) & is.na(RB1_mut) & is.na(TP53_mut), "cell_cycle_stable", "cell_cycle_alt"))
fisher.test(table(cell_cycle_df_noncodel$cell_cycle, cell_cycle_df_noncodel$acquired_aneuploidy))
wilcox.test(cell_cycle_df_noncodel$aneuploidy_diff~cell_cycle_df_noncodel$cell_cycle)

cell_cycle_df_noncodel %>% 
  group_by(cell_cycle) %>% 
  summarise(median_aneu =  median(aneuploidy_diff),
            mean_aneu = mean(aneuploidy_diff))


# Replace column names.
test2 = (test[ , c(45:52)])
rownames(test2) = test$tumor_pair_barcode
test2 = t(test2)
test3 = (test[ , c(53:60)])
rownames(test3) = test$tumor_pair_barcode
test3 = t(test3)
annotation = data.frame(aneuploidy = as.factor(test$acquired_aneuploidy))
rownames(annotation) <- colnames(test2)

pheatmap(test2, annotation = annotation, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE)
pheatmap(test3, annotation = annotation, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE)


# Must be able to enumerate the number of acquireed aneuploidy and cell cycle gain.


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
#  filter(arm_call != 0) %>% 
#  filter(!is.na(arm_call)) %>% 
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
  xlab("Taylor - AS score amplification") + ylab("GLASS - AS score amplification") + ggtitle("rho = 0.70, P-val < 2.4E-10 (n=43)") + xlim(0, 39) + ylim(0,39)
# Plot the AS_del scores vs. one another.
cor.test(taylor_glass$AS_del.x, taylor_glass$AS_del.y, method = "spearman")
ggplot(taylor_glass, aes(x = AS_del.x, y = AS_del.y)) + geom_point() + theme_bw() +
  xlab("Taylor - AS score deletion") + ylab("GLASS - AS score deletion") + ggtitle("rho = 0.90, P-val < 2.2E-16 (n=43)") + xlim(0, 39) + ylim(0, 39)


# Produce the same plots as before for aneuploidy_score, but use AS_del and AS_amp data.
complete_aneuploidy = aneuploidy_pairs %>% 
  left_join(arm_level_full, by=c("tumor_barcode_a"="aliquot_barcode")) %>% 
  inner_join(arm_level_full, by=c("tumor_barcode_b"="aliquot_barcode")) %>% 
  select(tumor_pair_barcode:AS_del.x, idh_codel_subtype:AS_del.y)  %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, "IDHmut_codel" = "IDHmut codel", "IDHmut_noncodel" = "IDHmut noncodel", "IDHwt_noncodel" = "IDHwt"))

# Create a plot for AneuploidyScore for just the deletions.
AS_del_plot = complete_aneuploidy %>% 
  gather(sample, AS_del, c(AS_del.x, AS_del.y)) %>% 
  mutate(sample =  recode(sample, "AS_del.x" = "primary", "AS_del.y" = "recurrence"))
ggplot(AS_del_plot, aes(x = sample, y = AS_del, group = tumor_pair_barcode)) +
  geom_line(linetype="solid", size=1, alpha = 0.3) + ylab("GLASS - Aneuploidy score (deletions)") + xlab("Gold set pairs (n=199)") +
  geom_point(color="black", size=2) + theme_bw() + facet_grid(~idh_codel_subtype)

# Create a plot for AneuploidyScore for just the amplifications.
AS_amp_plot = complete_aneuploidy %>% 
  gather(sample, AS_amp, c(AS_amp.x, AS_amp.y)) %>% 
  mutate(sample =  recode(sample, "AS_amp.x" = "primary", "AS_amp.y" = "recurrence"))
ggplot(AS_amp_plot, aes(x = sample, y = AS_amp, group = tumor_pair_barcode)) +
  geom_line(linetype="solid", size=1, alpha = 0.3) +  ylab("GLASS - Aneuploidy score (amplifications)")+ xlab("Gold set pairs (n=199)") +
  geom_point(color="black", size=2) + theme_bw() + facet_grid(~idh_codel_subtype)

# Group by IDH subtype.
AS_pairs_IDHwt = complete_aneuploidy %>% 
  filter(idh_codel_subtype == "IDHwt")
AS_pairs_IDH_codel = complete_aneuploidy %>% 
  filter(idh_codel_subtype == "IDHmut codel")
AS_pairs_IDH_noncodel = complete_aneuploidy %>% 
  filter(idh_codel_subtype == "IDHmut noncodel")

# Subtype specific analyses that now have filtered samples for Aneuploidy Score specific to deletions.
wilcox.test(AS_pairs_IDHwt$AS_del.x, AS_pairs_IDHwt$AS_del.y, paired = TRUE)
wilcox.test(AS_pairs_IDH_codel$AS_del.x, AS_pairs_IDH_codel$AS_del.y, paired = TRUE)
wilcox.test(AS_pairs_IDH_noncodel$AS_del.x, AS_pairs_IDH_noncodel$AS_del.y, paired = TRUE)

# And those specific to amplifications.
wilcox.test(AS_pairs_IDHwt$AS_amp.x, AS_pairs_IDHwt$AS_amp.y, paired = TRUE)
wilcox.test(AS_pairs_IDH_codel$AS_amp.x, AS_pairs_IDH_codel$AS_amp.y, paired = TRUE)
wilcox.test(AS_pairs_IDH_noncodel$AS_amp.x, AS_pairs_IDH_noncodel$AS_amp.y, paired = TRUE)

# Generate an ordered heatmap or grid for most commonly altered chromosome arms.
aneuploidy_heatmap = aneuploidy_pairs %>% 
  left_join(arm_level_full, by=c("tumor_barcode_a"="aliquot_barcode")) %>% 
  inner_join(arm_level_full, by=c("tumor_barcode_b"="aliquot_barcode")) %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, "IDHmut_codel" = "IDHmut codel", "IDHmut_noncodel" = "IDHmut noncodel", "IDHwt_noncodel" = "IDHwt")) %>% 
  arrange(idh_codel_subtype, tumor_pair_barcode)
# Replace column names.
colnames(aneuploidy_heatmap) = gsub("\\.x", "_P", colnames(aneuploidy_heatmap))
colnames(aneuploidy_heatmap) = gsub("\\.y", "_R", colnames(aneuploidy_heatmap))

# Rough sketch heatmap for the primary and recurrent tumors, separated by a whitespace for the three subtypes.
pheatmap(aneuploidy_heatmap[ , c(14:52)], cluster_rows = FALSE, cluster_cols = FALSE, gaps_row=c(24,90), show_rownames = FALSE)
pheatmap(aneuploidy_heatmap[ , c(58:96)], cluster_rows = FALSE, cluster_cols = FALSE, gaps_row=c(24,90), show_rownames = FALSE)

# 
aneuploidy_heatmap_noncodel_p = aneuploidy_heatmap %>% 
  filter(idh_codel_subtype=="IDHmut noncodel") %>% 
  select(-idh_codel_subtype_P, idh_codel_subtype_R) %>% 
  gather(primary_arm, primary_arm_call, c(`1p_P`:`22q_P`), -c(`1p_R`:`22q_R`)) 

aneuploidy_heatmap_noncodel_r = aneuploidy_heatmap %>% 
  filter(idh_codel_subtype=="IDHmut noncodel") %>% 
  select(-idh_codel_subtype_P, idh_codel_subtype_R) %>% 
  gather(recurrent_arm, recurrent_arm_call, c(`1p_R`:`22q_R`), -c(`1p_P`:`22q_P`)) 


# 66 samples that are IDHmut noncodels.
tmp1 = aneuploidy_heatmap_noncodel_p %>% 
  group_by(primary_arm, primary_arm_call) %>% 
  summarise(sample_number =  n()) 
tmp2 = aneuploidy_heatmap_noncodel_r %>% 
  group_by(recurrent_arm, recurrent_arm_call) %>% 
  summarise(sample_number =  n()) 



