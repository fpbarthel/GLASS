##################################################
# Visualize aneuploidy results based on GATK CNV
# Updated: 2019.01.14
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

##################################################
# Establish connection with database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

# Load in the blocklist to understand CN thresholds:
tumor_blocklist = dbReadTable(con,  Id(schema="analysis", table="blocklist"))

# Floris has created two aneuploidy metrics, 1. to represent the fraction of the genome with copy number alterations AND
# 2. An aneuploidy score to mimic that of Taylor et al. Cancer Cell 2018.

# Query that should return both aneuploidy estimates.
aneuploidy_final = dbGetQuery(con, "WITH t1 AS (
         SELECT gs.aliquot_barcode,
                              gs.chrom,
                              ca.arm,
                              ca.pos * gs.pos AS pos,
                              gs.log2_copy_ratio::numeric AS log2_copy_ratio,
                              gs.call,
                              upper(ca.pos * gs.pos)::numeric - lower(ca.pos * gs.pos)::numeric - 1::numeric AS seg_size,
                              sum(upper(ca.pos * gs.pos)::numeric - lower(ca.pos * gs.pos)::numeric - 1::numeric) OVER w AS arm_size,
                              lag(gs.call) OVER w2 = gs.call AS lag_call
                              FROM analysis.gatk_seg gs
                              JOIN ref.chr_arms ca ON ca.chrom = gs.chrom::bpchar AND ca.pos && gs.pos
                              WHERE gs.chrom::text <> ALL (ARRAY['X'::character varying, 'Y'::character varying]::text[])
                              WINDOW w AS (PARTITION BY gs.aliquot_barcode, ca.arm), w2 AS (PARTITION BY gs.aliquot_barcode, ca.arm ORDER BY (ca.pos * gs.pos))
), t2 AS (
                              SELECT t1.aliquot_barcode,
                              t1.arm,
                              sum(
                              CASE
                              WHEN t1.call = '+'::bpchar THEN t1.seg_size / t1.arm_size
                              ELSE 0::numeric
                              END) AS prop_amp,
                              sum(
                              CASE
                              WHEN t1.call = '-'::bpchar THEN t1.seg_size / t1.arm_size
                              ELSE 0::numeric
                              END) AS prop_del,
                              sum(
                              CASE
                              WHEN NOT t1.lag_call THEN 1
                              ELSE NULL::integer
                              END) AS n_call_changes
                              FROM t1
                              GROUP BY t1.aliquot_barcode, t1.arm
                              ORDER BY t1.aliquot_barcode, t1.arm
), t3 AS (
                              SELECT t2.aliquot_barcode,
                              sum(
                              CASE
                              WHEN (t2.prop_amp > 0.80 OR t2.prop_del > 0.80) THEN 1
                              ELSE 0
                              END) AS aneuploidy_score
                              FROM t2
                              GROUP BY t2.aliquot_barcode
)

SELECT ta.aliquot_barcode,aneuploidy_score,aneuploidy,idh_codel_subtype
FROM t3 ta
LEFT JOIN analysis.aneuploidy aa ON aa.aliquot_barcode = ta.aliquot_barcode
LEFT JOIN analysis.blocklist bl ON bl.aliquot_barcode = ta.aliquot_barcode
LEFT JOIN biospecimen.aliquots al ON al.aliquot_barcode = ta.aliquot_barcode
LEFT JOIN clinical.surgeries su ON su.sample_barcode = al.sample_barcode
WHERE bl.coverage_exclusion = 'allow' AND bl.cnv_exclusion <> 'block' AND bl.fingerprint_exclusion = 'allow'")

# Add metadata to merge with other datasets.
aneuploidy_final = aneuploidy_final %>% 
  mutate_if(bit64::is.integer64, as.double) %>% 
  mutate(sample_barcode = substr(aliquot_barcode, 1, 15),
         seq_type = substr(aliquot_barcode, 21, 23))

# Examine the relationship between the two metrics for all samples.
cor.test(aneuploidy_final$aneuploidy_score, aneuploidy_final$aneuploidy, method = "spearman")
ggplot(aneuploidy_final, aes(x = aneuploidy, y = aneuploidy_score)) + geom_point() + theme_bw() +ggtitle("rho = 0.84")

# Taylor et al. provided aneuploidy scores as part of their paper. Extract available LGG and GBM samples that overlap
# with GLASS. Note: The aneuploidy data in the paper was generated using ABSOLUTE for array-based data.
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

# Create plot of two aneuploidy scores vs. one another. Double-check whether they consider all the same chromosomal arms.
cor.test(glass_taylor_aneuploidy$`AneuploidyScore(AS)`, glass_taylor_aneuploidy$aneuploidy_score, method = "spearman")
ggplot(glass_taylor_aneuploidy, aes(x = `AneuploidyScore(AS)`, y = aneuploidy_score, color = seq_type)) + geom_point() + theme_bw() +
  ylab("GLASS - aneuploidy score") + xlab("Taylor et al. AneuploidyScore(AS)") + ggtitle("rho = 0.86")

# Examine the relationship between the fraction of the genome altered aneuploidy score and AneuploidyScores.
cor.test(glass_taylor_aneuploidy$`AneuploidyScore(AS)`, glass_taylor_aneuploidy$aneuploidy, method = "spearman")
ggplot(glass_taylor_aneuploidy, aes(x = `AneuploidyScore(AS)`, y = aneuploidy, color = seq_type)) + geom_point() + theme_bw() +
  ylab("GLASS - aneuploidy value") + xlab("Taylor et al. AneuploidyScore(AS)") + ggtitle("rho = 0.77")

# Use this information in the context of the silver set for aneuploidy pairs. Note: that this includes some poor CN data.
aneuploidy_pairs = dbGetQuery(con, "SELECT tumor_pair_barcode, tumor_barcode_a, tumor_barcode_b, a1.aneuploidy AS aneuploidy_a, a2.aneuploidy AS aneuploidy_b, t1.aneuploidy_score::integer AS aneuploidy_score_a,t2.aneuploidy_score::integer AS aneuploidy_score_b,idh_codel_subtype
                              FROM analysis.silver_set ss
                              LEFT JOIN analysis.aneuploidy a1 ON a1.aliquot_barcode = ss.tumor_barcode_a
                              LEFT JOIN analysis.aneuploidy a2 ON a2.aliquot_barcode = ss.tumor_barcode_b
                              LEFT JOIN analysis.taylor_aneuploidy t1 ON t1.aliquot_barcode = ss.tumor_barcode_a
                              LEFT JOIN analysis.taylor_aneuploidy t2 ON t2.aliquot_barcode = ss.tumor_barcode_b
                              LEFT JOIN biospecimen.aliquots al ON al.aliquot_barcode = ss.tumor_barcode_a
                              LEFT JOIN clinical.surgeries su ON su.sample_barcode = al.sample_barcode")

# Filter out any poor performing samples on CNV blocklist.
aneuploidy_pairs_filtered = aneuploidy_pairs %>% 
  inner_join(tumor_blocklist, by=c("tumor_barcode_a"="aliquot_barcode")) %>% 
  inner_join(tumor_blocklist, by=c("tumor_barcode_b"="aliquot_barcode")) %>% 
  # Note remove trailing whitespace.
  filter(cnv_exclusion.x %in%c("allow ", "review")) %>% 
  filter(cnv_exclusion.y %in% c("allow ", "review")) 

# Define aneuploidy pairs by IDH subtype.
aneuploidy_pairs_IDHwt = aneuploidy_pairs_filtered %>% 
  filter(idh_codel_subtype == "IDHwt_noncodel")
aneuploidy_pairs_IDH_codel = aneuploidy_pairs_filtered %>% 
  filter(idh_codel_subtype == "IDHmut_codel")
aneuploidy_pairs_IDH_noncodel = aneuploidy_pairs_filtered %>% 
  filter(idh_codel_subtype == "IDHmut_noncodel")

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
  geom_line(linetype="solid", size=1, alpha = 0.3) + ylab("Aneuploidy value") + xlab("Gold set pairs (n=197)") +
  geom_point(color="black", size=2) + theme_bw() + facet_grid(~idh_codel_subtype)

# Combine the aneuploidy scores for each tumor pair.
aneuploid_plot_score = aneuploidy_pairs_filtered %>% 
  gather(sample, aneuploidy_score, c(aneuploidy_score_a, aneuploidy_score_b)) 
ggplot(aneuploid_plot_score, aes(x = sample, y = aneuploidy_score, group = tumor_pair_barcode)) +
  geom_line(linetype="solid", size=1, alpha = 0.3) + ylab("Aneuploidy score") + xlab("Gold set pairs (n=197)") +
  geom_point(color="black", size=2) + theme_bw() + facet_grid(~idh_codel_subtype)