##################################################
# Heatmap of arm-level copy number calls - GLASS
# Updated: 2019.05.19
# Author: Kevin J.
##################################################

# Working directory for this analysis in the GLASS-analysis project. 
mybasedir = "/Volumes/verhaak-lab/GLASS-analysis/"
setwd(mybasedir)

##################################################

# Necessary packages:
library(tidyverse)
library(DBI)
library(pheatmap)
library(RColorBrewer)
library(openxlsx)

##################################################
# Establish connection with database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB2")

# Investigate the arm-level calls, compare with Taylor et al.
arm_level_calls = dbGetQuery(con, "SELECT * FROM analysis.gatk_cnv_by_arm")
gold_set = dbReadTable(con,  Id(schema="analysis", table="gold_set"))
surgeries = dbReadTable(con,  Id(schema="clinical", table="surgeries"))

# Define aneuploidy pairs from the gold set.
aneuploidy_pairs <- dbGetQuery(con, "SELECT tumor_pair_barcode, case_barcode, tumor_barcode_a, tumor_barcode_b, a1.prop_aneuploidy AS aneuploidy_a, a1.aneuploidy_score::integer AS aneuploidy_score_a, 
        a1.aneuploidy_amp_score::integer AS aneuploidy_amp_score_a, a1.aneuploidy_del_score::integer AS aneuploidy_del_score_a, 
                               a2.prop_aneuploidy AS aneuploidy_b, a2.aneuploidy_score::integer AS aneuploidy_score_b, a2.aneuploidy_amp_score::integer AS aneuploidy_amp_score_b, a2.aneuploidy_del_score::integer AS aneuploidy_del_score_b
                               FROM analysis.gold_set gs
                               LEFT JOIN analysis.gatk_aneuploidy a1 ON a1.aliquot_barcode = gs.tumor_barcode_a
                               LEFT JOIN analysis.gatk_aneuploidy a2 ON a2.aliquot_barcode = gs.tumor_barcode_b")


# Only use arm calls from the gold set. Note individual chromosome arms may not be resolved.
arm_level_calls_filtered = arm_level_calls %>% 
  filter(aliquot_barcode %in% c(gold_set$tumor_barcode_a, gold_set$tumor_barcode_b))

# Creating summary variables for each sample. AneuploidyScore for amplifications and deletions.
arm_level_summary_1 =  arm_level_calls_filtered %>% 
  group_by(aliquot_barcode, arm_call) %>% 
  summarise(events = n()) %>% 
  spread(arm_call, events) %>% 
  mutate(sample_barcode = substr(aliquot_barcode, 1, 15)) %>% 
  inner_join(surgeries, by="sample_barcode") %>% 
  select(aliquot_barcode, sample_barcode, surgery_number, idh_codel_subtype, AS_amp = `1`, AS_del = `-1`) 

# Reduce set to only those arm-level analyses considered by Taylor et al.
arm_level_summary2 = arm_level_calls_filtered %>% 
  select(-chrom, -arm_num_seg, -arm_cr_wmean, -arm_cr_wsd) %>% 
  spread(arm, arm_call) %>% 
  select(aliquot_barcode, `1p`, `1q`, `2p`, `2q`, `3p`, `3q`, `4p`, `4q`, `5p`, `5q`, `6p`, `6q`, `7p`, `7q`, `8p`,
         `8q`, `9p`, `9q`, `10p`,`10q`, `11p`, `11q`, `12p`, `12q`, `13q`, `14q`, `15q`, `16q`, `16p`, `17p`, `17q`,
         `18p`, `18q`, `19p`, `19q`, `20p`,`20q`, `21q`, `22q`) 

# The CN values should be filtered so merge with the gold_set. Make data long.
gold_set_merge = gold_set %>% 
  gather(sample_type, tumor_barcode, c(tumor_barcode_a, tumor_barcode_b))

# Combine two tables.
arm_level_full = arm_level_summary_1 %>% 
  inner_join(arm_level_summary2, by="aliquot_barcode") %>% 
  inner_join(gold_set_merge, by=c("aliquot_barcode"="tumor_barcode")) %>% 
  select(aliquot_barcode:`22q`) %>% 
  mutate(AS_amp = ifelse(is.na(AS_amp), 0, AS_amp),
         AS_del = ifelse(is.na(AS_del), 0, AS_del))

# Taylor et al. provided aneuploidy scores as part of their paper. Extract available LGG and GBM samples that overlap
# with GLASS. Note: The aneuploidy data in the paper was generated using ABSOLUTE for array-based data in TCGA.
taylor_data = readWorkbook("data/published/taylor-cancer-cell-2018.xlsx", sheet = 1, rowNames = F, colNames = TRUE)

# Extract only the iniformation for GBM and LGG. Recode some of their metadata to fit our own.
glioma_dat = taylor_data %>% 
  filter(Type%in%c("GBM", "LGG")) %>% 
  mutate(sample_type_num = substr(Sample, 14, 15), 
         sample_type = recode(sample_type_num, "01" = "TP", "02" = "R1"),
         sample_barcode = paste(substr(Sample, 1, 12), sample_type, sep ="-"))

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

# Finally, analyze the Aneuploidy score at the arm-level.
taylor_glass = glioma_dat %>% 
  mutate(dataset = "taylor") %>% 
  select(aliquot_barcode = Sample, sample_barcode, dataset, AS_del, AS_amp, `1p`:`22q`) %>% 
  filter(sample_barcode%in%arm_level_merge$sample_barcode) %>% 
  inner_join(arm_level_merge, by="sample_barcode")

# Plot the **AS_amp scores vs. one another.
cor.test(taylor_glass$AS_amp.x, taylor_glass$AS_amp.y, method = "spearman")
ggplot(taylor_glass, aes(x = AS_amp.x, y = AS_amp.y)) + geom_point() + theme_bw() +
  xlab("Taylor - AS score amplification") + ylab("GLASS - AS score amplification") + ggtitle("rho = 0.77, P = 2.6E-09 (n=41)") + xlim(0, 39) + ylim(0,39)
# Plot the **AS_del scores vs. one another.
cor.test(taylor_glass$AS_del.x, taylor_glass$AS_del.y, method = "spearman")
ggplot(taylor_glass, aes(x = AS_del.x, y = AS_del.y)) + geom_point() + theme_bw() +
  xlab("Taylor - AS score deletion") + ylab("GLASS - AS score deletion") + ggtitle("rho = 0.88, P = 2.2E-14 (n=41)") + xlim(0, 39) + ylim(0, 39)

# Produce the same plots as before for aneuploidy_score, but use AS_del and AS_amp data.
complete_aneuploidy = aneuploidy_pairs %>% 
  left_join(arm_level_full, by=c("tumor_barcode_a"="aliquot_barcode")) %>% 
  inner_join(arm_level_full, by=c("tumor_barcode_b"="aliquot_barcode")) %>% 
  select(-idh_codel_subtype.y) %>% 
  select(tumor_pair_barcode:surgery_number.x, idh_codel_subtype = idh_codel_subtype.x, AS_amp.x:`22q.y`)

# Replace column names with names representative of primary ("_a") and recurrence ("_b").
colnames(complete_aneuploidy) = gsub("\\.x", "_a", colnames(complete_aneuploidy))
colnames(complete_aneuploidy) = gsub("\\.y", "_b", colnames(complete_aneuploidy))

# Create a plot for AneuploidyScore for just the deletions.
AS_del_plot = complete_aneuploidy %>% 
  gather(sample, AS_del, c(AS_del_a, AS_del_b)) %>% 
  mutate(sample =  recode(sample, "AS_del_a" = "primary", "AS_del_b" = "recurrence"))
ggplot(AS_del_plot, aes(x = sample, y = AS_del, group = tumor_pair_barcode)) +
  geom_line(linetype="solid", size=1, alpha = 0.3) + ylab("GLASS - Aneuploidy score (deletions)") + xlab("") +
  geom_point(color="black", size=2) + theme_bw() + facet_grid(~idh_codel_subtype)

# Create a plot for AneuploidyScore for just the amplifications.
AS_amp_plot = complete_aneuploidy %>% 
  gather(sample, AS_amp, c(AS_amp_a, AS_amp_b)) %>% 
  mutate(sample =  recode(sample, "AS_amp_a" = "primary", "AS_amp_b" = "recurrence"))
ggplot(AS_amp_plot, aes(x = sample, y = AS_amp, group = tumor_pair_barcode)) +
  geom_line(linetype="solid", size=1, alpha = 0.3) +  ylab("GLASS - Aneuploidy score (amplifications)")+ xlab("") +
  geom_point(color="black", size=2) + theme_bw() + facet_grid(~idh_codel_subtype)

# Generate an ordered heatmap or grid for most commonly altered chromosome arms.
aneuploidy_heatmap = complete_aneuploidy %>% 
  arrange(idh_codel_subtype, tumor_pair_barcode)

# Rough heatmap for the primary and recurrent tumors, separated by a whitespace for the three subtypes.
# Samples will be ordered in the same way across primary and recurrence.
pdf(file = "/Users/johnsk/Documents/cnv-arm-heatmap-primary.pdf", height = 5, width = 7, bg = "transparent", useDingbats = FALSE)
pheatmap(aneuploidy_heatmap[ , c(18:56)], cluster_rows = FALSE, cluster_cols = FALSE, gaps_row=c(25,89), show_rownames = FALSE)
dev.off()
pdf(file = "/Users/johnsk/Documents/cnv-arm-heatmap-recurrence.pdf", height = 5, width = 7, bg = "transparent", useDingbats = FALSE)
pheatmap(aneuploidy_heatmap[ , c(61:99)], cluster_rows = FALSE, cluster_cols = FALSE, gaps_row=c(25,89), show_rownames = FALSE)
dev.off()


