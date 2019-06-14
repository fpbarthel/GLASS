####################################################
# CDKN2A, Aneuploidy timing analysis
####################################################

####################################################
# Do necessary settings
####################################################
# Load necessary packages
library(tidyverse)
library(DBI)
library(ggpubr)

# Create connection to database and load necessary tables
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

# Load in variables from the database.
gold_set = dbReadTable(con,  Id(schema="analysis", table="gold_set"))

#modify goldset
gold <- gold_set %>% gather(key = class, value = aliquot_barcode, -tumor_pair_barcode, -case_barcode) %>% mutate(class = ifelse(class == "tumor_barcode_a", "P", "R"))
# load in CCF for CDKN2A
CCF_query <- ("WITH
seg_wccf AS
                       (
                       SELECT
                       aliquot_barcode,
                       COUNT(*) AS num_seg,
                       sum(CASE WHEN cellular_prevalence IS NOT NULL THEN (upper(pos) - lower(pos) -1) * cellular_prevalence END) / sum(CASE WHEN cellular_prevalence IS NOT NULL THEN upper(pos) - lower(pos) -1 END)::decimal AS wccf,
                       max(clonal_cluster) AS num_clusters,
                       sum(CASE WHEN copy_number <> 2 THEN upper(pos) - lower(pos) -1 END) / sum(upper(pos) - lower(pos) -1)::decimal AS titan_aneuploidy
                       FROM variants.titan_seg
                       GROUP BY 1
                       ),
                       cdkn2a_call AS
                       (
                       SELECT aliquot_barcode, hlvl_call AS cdkn2a_call, cellular_prevalence AS cdkn2a_ccf
                       FROM analysis.gatk_cnv_by_gene
                       WHERE gene_symbol = 'CDKN2A'
                       ),
                       selected_regions AS
                       (
                       SELECT '10q25-26' AS region, * FROM ref.cytobands WHERE chrom = 10 AND substring(cytoband from 1 for 3) IN ('q25','q26')
                       ),
                       gene_seg_intersect AS
                       (
                       SELECT aliquot_barcode, region, gs.chrom, (upper(t0.pos * gs.pos) - lower(t0.pos * gs.pos) -1) AS w, 2^log2_copy_ratio::decimal As cr
                       FROM variants.gatk_seg gs
                       INNER JOIN selected_regions t0 ON t0.chrom = gs.chrom AND t0.pos && gs.pos
                       ),
                       gene_sample_call AS
                       (
                       SELECT aliquot_barcode, region, 
                       sum(w * cr) / sum(w) AS wcr
                       FROM gene_seg_intersect
                       GROUP BY aliquot_barcode, region
                       ),
                       seg_stats_optimized AS
                       (
                       SELECT
                       gs.aliquot_barcode,
                       LEAST(0.9, neu_fwmean - 2 * neu_fwsd) AS del_thres,
                       GREATEST(1.1, neu_fwmean + 2 * neu_fwsd) AS amp_thres,
                       (CASE
                       WHEN max_loss_arm_wmean < 0.9 AND max_loss_arm_n >= 3 THEN GREATEST(0,max_loss_arm_wmean - 2 * max_loss_arm_wsd)
                       WHEN del_fwmean < 0.9 AND del_n >= 3 THEN GREATEST(0,del_fwmean - 2 * del_fwsd)
                       ELSE NULL
                       END) AS hldel_thres,
                       (CASE
                       WHEN max_gain_arm_wmean > 1.1 AND max_gain_arm_n >= 3 THEN max_gain_arm_wmean + 2 * max_gain_arm_wsd
                       WHEN amp_fwmean > 1.1 AND amp_n >= 3 THEN amp_fwmean + 2 * amp_fwsd
                       ELSE NULL
                       END) AS hlamp_thres
                       FROM analysis.gatk_seg_stats gs
                       LEFT JOIN analysis.gatk_aneuploidy gsa ON gsa.aliquot_barcode = gs.aliquot_barcode
                       ),
                       gene_cp AS
                       (
                       SELECT ts.aliquot_barcode, region, ts.chrom, (upper(t0.pos * ts.pos) - lower(t0.pos * ts.pos) -1) AS w, cellular_prevalence As cp
                       FROM variants.titan_seg ts
                       INNER JOIN selected_regions t0 ON t0.chrom = ts.chrom AND t0.pos && ts.pos
                       ),
                       gene_cp_agg AS
                       (
                       SELECT aliquot_barcode, region, 
                       COALESCE(sum(w * cp) / NULLIF(sum(w),0),NULL) AS wcp
                       FROM gene_cp
                       GROUP BY 1, 2
                       ),
                       calls AS
                       (
                       SELECT
                       gc.aliquot_barcode,
                       gc.region,
                       (CASE
                       WHEN gc.wcr >= del_thres AND gc.wcr <= amp_thres THEN 0
                       WHEN gc.wcr < hldel_thres THEN -2
                       WHEN gc.wcr < del_thres THEN -1
                       WHEN gc.wcr > hlamp_thres THEN 2
                       WHEN gc.wcr > amp_thres THEN 1
                       ELSE NULL
                       END) c10q25_26_call,
                       gc.wcr,
                       wcp AS c10q25_26_ccf
                       FROM gene_sample_call gc
                       LEFT JOIN seg_stats_optimized ss ON ss.aliquot_barcode = gc.aliquot_barcode
                       LEFT JOIN gene_cp_agg cp ON cp.aliquot_barcode = gc.aliquot_barcode AND cp.region = gc.region
                       ORDER BY 3
                       )
                       SELECT seg.aliquot_barcode, idh_codel_subtype, num_seg, num_clusters, wccf, cnv_exclusion, cdkn2a_call, cdkn2a_ccf, titan_aneuploidy, prop_aneuploidy, aneuploidy_score, c10q25_26_call, c10q25_26_ccf
                       FROM seg_wccf seg
                       INNER JOIN cdkn2a_call cc ON cc.aliquot_barcode = seg.aliquot_barcode
                       INNER JOIN analysis.gatk_aneuploidy ga ON ga.aliquot_barcode = seg.aliquot_barcode
                       INNER JOIN clinical.subtypes su ON su.case_barcode = substring(seg.aliquot_barcode from 1 for 12)
                       INNER JOIN analysis.blocklist bl ON bl.aliquot_barcode = seg.aliquot_barcode
                       INNER JOIN calls ca ON ca.aliquot_barcode = seg.aliquot_barcode
                       ORDER BY 1")
CCF_df <- dbGetQuery(con, CCF_query)


CCF <- CCF_df %>% 
  mutate(cnv_exclusion = trimws(cnv_exclusion, "both"), aneuploidy_score = as.numeric(aneuploidy_score), num_clusters = as.character(num_clusters)) %>% 
  inner_join(gold, by = "aliquot_barcode")

CCF_cdkn2a_aneuploidy <- CCF %>% 
  filter(cnv_exclusion != "block", cdkn2a_call == -2, !is.na(cdkn2a_ccf), idh_codel_subtype != "IDHmut-codel") %>% 
  mutate(diff1 = cdkn2a_ccf - wccf,
         cdkn2a_timing_v1 = case_when(cdkn2a_ccf - wccf <= -0.05  ~ "aneuploidy_precedes", 
                                      cdkn2a_ccf - wccf >= 0.05 ~ "cdkn2a_precedes",
                                      is.na(cdkn2a_ccf) ~ "NA",
                                      TRUE ~ "simultaneous"),
          cdkn2a_clonality = ifelse(cdkn2a_ccf >= 0.5,"clonal","subclonal"), 
         genome_clonality = ifelse(wccf >= 0.5, "clonal", "subclonal"), 
         cdkn2a_timing_v2 = case_when(cdkn2a_clonality == "clonal" & genome_clonality == "clonal" ~ "both_clonal", 
                                      cdkn2a_clonality == "clonal" & genome_clonality == "subclonal" ~ "early", 
                                      cdkn2a_clonality == "subclonal" & genome_clonality == "clonal" ~ "late", 
                                      cdkn2a_clonality == "subclonal" & genome_clonality == "subclonal" ~ "both_subclonal")) 

#tables
table(CCF_cdkn2a_aneuploidy$idh_codel_subtype, CCF_cdkn2a_aneuploidy$cdkn2a_timing_v1)

#plot differences
ggplot(CCF_cdkn2a_aneuploidy,aes(x = diff1)) +
  geom_density() +
  facet_grid(~idh_codel_subtype)


# build subtypes
CCF_cdkn2a_aneuploidy_noncodel <- CCF_cdkn2a_aneuploidy %>% filter(idh_codel_subtype == "IDHmut-noncodel")
CCF_cdkn2a_aneuploidy_wt <- CCF_cdkn2a_aneuploidy %>% filter(idh_codel_subtype == "IDHwt")
CCF_cdkn2a_aneuploidy_wt_initial <- CCF_cdkn2a_aneuploidy_wt %>% filter(class == "P")
CCF_cdkn2a_aneuploidy_wt_recurrence <- CCF_cdkn2a_aneuploidy_wt %>% filter(class == "R")
CCF_cdkn2a_aneuploidy_wt_shared <- gold_set %>% rename(aliquot_barcode = tumor_barcode_a) %>% full_join(CCF_cdkn2a_aneuploidy_wt, by = c("case_barcode", "aliquot_barcode")) %>% 
  rename(tumor_barcode_a = aliquot_barcode, aliquot_barcode = tumor_barcode_b) %>% full_join(CCF_cdkn2a_aneuploidy_wt, by = c("case_barcode", "aliquot_barcode", "idh_codel_subtype")) %>%
  filter(idh_codel_subtype == "IDHwt", cdkn2a_call.x == -2 & cdkn2a_call.y == -2) %>% rename(tumor_barcode_b = aliquot_barcode)

# statistical tests separately for subtypes
wilcox.test(CCF_cdkn2a_aneuploidy_noncodel$cdkn2a_ccf, CCF_cdkn2a_aneuploidy_noncodel$wccf, paired = T)
wilcox.test(CCF_cdkn2a_aneuploidy_wt$cdkn2a_ccf, CCF_cdkn2a_aneuploidy_wt$wccf, paired = T)
wilcox.test(CCF_cdkn2a_aneuploidy_wt_initial$cdkn2a_ccf, CCF_cdkn2a_aneuploidy_wt_initial$wccf, paired = T) 
wilcox.test(CCF_cdkn2a_aneuploidy_wt_recurrence$cdkn2a_ccf, CCF_cdkn2a_aneuploidy_wt_recurrence$wccf, paired = T)
wilcox.test(CCF_cdkn2a_aneuploidy_wt_shared$cdkn2a_ccf.x, CCF_cdkn2a_aneuploidy_wt_shared$wccf.x, paired = T) 
wilcox.test(CCF_cdkn2a_aneuploidy_wt_shared$cdkn2a_ccf.y, CCF_cdkn2a_aneuploidy_wt_shared$wccf.y, paired = T)

# plot changes: cdkn2a vs aneuploidy
nc1 <- CCF_cdkn2a_aneuploidy_noncodel %>% 
  select(aliquot_barcode, idh_codel_subtype, cdkn2a_timing_v1, wccf, cdkn2a_ccf) %>% 
  gather(key=type, value = ccf, wccf:cdkn2a_ccf) %>% 
  ggplot(aes(x = type, y = ccf, color = cdkn2a_timing_v1, group = aliquot_barcode)) +
  geom_line(linetype="solid", size=1, alpha = 0.3) + ylab("Cancer Cell Fraction") + xlab("") + ggtitle("IDHmut-noncodel")+ theme_bw() +
  geom_point(size=2) + annotate("text", x=1.5, y=1.03, label= "P = 0.017", fontface= "italic", size = 5) +
  scale_x_discrete(labels=c("cdkn2a_ccf" = "CDKN2A\nhomdel", "wccf" = "CNV\ngenome-wide")) +
  scale_color_manual(values=c("cdkn2a_precedes" = "#FF7100", "aneuploidy_precedes" = "#00B945", "simultaneous" = "black", "distal10q_precedes" = "#008EFF"), 
                      name="Order of events",
                      labels=c("cdkn2a_precedes" = "CDKN2A-homdel precedes", "aneuploidy_precedes" = "Aneuploidy precedes", "simultaneous" = "Simultaneous", "distal10q_precedes" = "10q25-26-del precedes"))

wt1 <- CCF_cdkn2a_aneuploidy_wt %>% 
  select(aliquot_barcode, idh_codel_subtype, cdkn2a_timing_v1, wccf, cdkn2a_ccf) %>% 
  gather(key=type, value = ccf, wccf:cdkn2a_ccf) %>% 
  ggplot(aes(x = type, y = ccf, color = cdkn2a_timing_v1, group = aliquot_barcode)) +
  geom_line(linetype="solid", size=1, alpha = 0.3) + ylab("Cancer Cell Fraction") + xlab("") + ggtitle("IDHwt")+ theme_bw() +
  geom_point(size=2) + annotate("text", x=1.5, y=1.03, label= "P = 1.957e-07", fontface= "italic", size = 5) +
  scale_x_discrete(labels=c("cdkn2a_ccf" = "CDKN2A\nhomdel", "wccf" = "CNV\ngenome-wide")) +
  scale_color_manual(values=c("cdkn2a_precedes" = "#FF7100", "aneuploidy_precedes" = "#00B945", "simultaneous" = "black", "distal10q_precedes" = "#008EFF"), 
                     name="Order of events",
                     labels=c("cdkn2a_precedes" = "CDKN2A-homdel precedes", "aneuploidy_precedes" = "Aneuploidy precedes", "simultaneous" = "Simultaneous", "distal10q_precedes" = "10q25-26-del precedes"))

wt1_initial <- CCF_cdkn2a_aneuploidy_wt_initial %>% 
  select(aliquot_barcode, idh_codel_subtype, cdkn2a_timing_v1, wccf, cdkn2a_ccf) %>% 
  gather(key=type, value = ccf, wccf:cdkn2a_ccf) %>% 
  ggplot(aes(x = type, y = ccf, color = cdkn2a_timing_v1, group = aliquot_barcode)) +
  geom_line(linetype="solid", size=1, alpha = 0.3) + ylab("Cancer Cell Fraction") + xlab("") + ggtitle("IDHwt - Initial")+ theme_bw() +
  geom_point(size=2) + annotate("text", x=1.5, y=1.03, label= "P = 3.921e-03", fontface= "italic", size = 5) +
  scale_x_discrete(labels=c("cdkn2a_ccf" = "CDKN2A\nhomdel", "wccf" = "CNV\ngenome-wide")) +
  scale_color_manual(values=c("cdkn2a_precedes" = "#FF7100", "aneuploidy_precedes" = "#00B945", "simultaneous" = "black", "distal10q_precedes" = "#008EFF"), 
                     name="Order of events",
                     labels=c("cdkn2a_precedes" = "CDKN2A-homdel precedes", "aneuploidy_precedes" = "Aneuploidy precedes", "simultaneous" = "Simultaneous", "distal10q_precedes" = "10q25-26-del precedes"))

wt1_recurrence <- CCF_cdkn2a_aneuploidy_wt_recurrence %>% 
  select(aliquot_barcode, idh_codel_subtype, cdkn2a_timing_v1, wccf, cdkn2a_ccf) %>% 
  gather(key=type, value = ccf, wccf:cdkn2a_ccf) %>% 
  ggplot(aes(x = type, y = ccf, color = cdkn2a_timing_v1, group = aliquot_barcode)) +
  geom_line(linetype="solid", size=1, alpha = 0.3) + ylab("Cancer Cell Fraction") + xlab("") + ggtitle("IDHwt - Recurrence")+ theme_bw() +
  geom_point(size=2) + annotate("text", x=1.5, y=1.03, label= "P = 1.608e-05", fontface= "italic", size = 5) +
  scale_x_discrete(labels=c("cdkn2a_ccf" = "CDKN2A\nhomdel", "wccf" = "CNV\ngenome-wide")) +
  scale_color_manual(values=c("cdkn2a_precedes" = "#FF7100", "aneuploidy_precedes" = "#00B945", "simultaneous" = "black", "distal10q_precedes" = "#008EFF"), 
                     name="Order of events",
                     labels=c("cdkn2a_precedes" = "CDKN2A-homdel precedes", "aneuploidy_precedes" = "Aneuploidy precedes", "simultaneous" = "Simultaneous", "distal10q_precedes" = "10q25-26-del precedes"))

wt1_shared_initial <- CCF_cdkn2a_aneuploidy_wt_shared %>% 
  select(tumor_barcode_a, idh_codel_subtype, cdkn2a_timing_v1.x, wccf.x, cdkn2a_ccf.x) %>% 
  gather(key=type, value = ccf, wccf.x:cdkn2a_ccf.x) %>% 
  ggplot(aes(x = type, y = ccf, color = cdkn2a_timing_v1.x, group = tumor_barcode_a)) +
  geom_line(linetype="solid", size=1, alpha = 0.3) + ylab("Cancer Cell Fraction") + xlab("") + ggtitle("IDHwt - Initial (only shared)")+ theme_bw() +
  geom_point(size=2) + annotate("text", x=1.5, y=1.03, label= "P = 0.01365", fontface= "italic", size = 5) +
  scale_x_discrete(labels=c("cdkn2a_ccf.x" = "CDKN2A\nhomdel", "wccf.x" = "CNV\ngenome-wide")) +
  scale_color_manual(values=c("cdkn2a_precedes" = "#FF7100", "aneuploidy_precedes" = "#00B945", "simultaneous" = "black", "distal10q_precedes" = "#008EFF"), 
                     name="Order of events",
                     labels=c("cdkn2a_precedes" = "CDKN2A-homdel precedes", "aneuploidy_precedes" = "Aneuploidy precedes", "simultaneous" = "Simultaneous", "distal10q_precedes" = "10q25-26-del precedes"))

wt1_shared_recurrence <- CCF_cdkn2a_aneuploidy_wt_shared %>% 
  select(tumor_barcode_b, idh_codel_subtype, cdkn2a_timing_v1.y, wccf.y, cdkn2a_ccf.y) %>% 
  gather(key=type, value = ccf, wccf.y:cdkn2a_ccf.y) %>% 
  ggplot(aes(x = type, y = ccf, color = cdkn2a_timing_v1.y, group = tumor_barcode_b)) +
  geom_line(linetype="solid", size=1, alpha = 0.3) + ylab("Cancer Cell Fraction") + xlab("") + ggtitle("IDHwt - recurrence (only shared)")+ theme_bw() +
  geom_point(size=2) + annotate("text", x=1.5, y=1.03, label= "P = 3.297e-03", fontface= "italic", size = 5) +
  scale_x_discrete(labels=c("cdkn2a_ccf.y" = "CDKN2A\nhomdel", "wccf.y" = "CNV\ngenome-wide")) +
  scale_color_manual(values=c("cdkn2a_precedes" = "#FF7100", "aneuploidy_precedes" = "#00B945", "simultaneous" = "black", "distal10q_precedes" = "#008EFF"), 
                     name="Order of events",
                     labels=c("cdkn2a_precedes" = "CDKN2A-homdel precedes", "aneuploidy_precedes" = "Aneuploidy precedes", "simultaneous" = "Simultaneous", "distal10q_precedes" = "10q25-26-del precedes"))


######## arrange multiple plots into one figure / export figures
pdf(file = "/Users/c-kocake/Box Sync/Projects/GLASS/Timing/timing_cdkn2a_noncodel.pdf", height = 5, width = 5, bg = "transparent", useDingbats = FALSE)
nc1
dev.off()

pdf(file = "/Users/c-kocake/Box Sync/Projects/GLASS/Timing/timing_cdkn2a_wt.pdf", height = 5, width = 9, bg = "transparent", useDingbats = FALSE)
ggarrange(wt1_initial, wt1_recurrence, common.legend = T, legend = "right")
dev.off()

pdf(file = "/Users/c-kocake/Box Sync/Projects/GLASS/Timing/timing_cdkn2a_shared_wt.pdf", height = 5, width = 9, bg = "transparent", useDingbats = FALSE)
ggarrange(wt1_shared_initial, wt1_shared_recurrence, common.legend = T, legend = "right")
dev.off()

###### END ######