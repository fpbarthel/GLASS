library(tidyverse)
library(odbc)
library(DBI)
library(ggthemes)
library(ggplot2)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

q = "SELECT tmc.*, sa.sample_type AS sample_type_a, sb.sample_type AS sample_type_b, ua.surgery_number AS surgery_number_a, ub.surgery_number AS surgery_number_b, ua.surgical_interval_mo AS surgical_interval_a, ub.surgical_interval_mo AS surgical_interval_b, ua.idh_codel_subtype,
ma.cumulative_coverage AS coverage_a, mb.cumulative_coverage AS coverage_b, ma.coverage_adj_mut_freq AS mf_a, mb.coverage_adj_mut_freq AS mf_b
FROM analysis.tumor_mut_comparison tmc
INNER JOIN biospecimen.aliquots aa ON tmc.tumor_a_barcode = aa.aliquot_barcode
INNER JOIN biospecimen.aliquots ab ON tmc.tumor_b_barcode = ab.aliquot_barcode
INNER JOIN biospecimen.samples sa ON aa.sample_barcode = sa.sample_barcode
INNER JOIN biospecimen.samples sb ON ab.sample_barcode = sb.sample_barcode
LEFT JOIN clinical.surgeries ua ON ua.sample_barcode = sa.sample_barcode
LEFT JOIN clinical.surgeries ub ON ub.sample_barcode = sb.sample_barcode
LEFT JOIN analysis.mutation_freq ma ON ma.aliquot_barcode = aa.aliquot_barcode
LEFT JOIN analysis.mutation_freq mb ON mb.aliquot_barcode = ab.aliquot_barcode"

qres <- dbGetQuery(con, q)

df = qres %>% 
  filter(mf_a > 0, mf_b > 0) %>% ## (1) Filter out hypermutator samples
  filter(sample_type_a == "TP", sample_type_b == "R1") %>%
  mutate(prop_a = setdiff_a / union_ab,
         prop_b = setdiff_b / union_ab,
         prop_ab = intersection_ab / union_ab,
         adj_a = 1e6*(setdiff_a / coverage_a),
         adj_b = 1e6*(setdiff_b / coverage_b),
         adj_ab = 1e6*(intersection_ab / ((coverage_a+coverage_b)/2)),
         timediff = abs(surgical_interval_b - surgical_interval_a),
         pertime_a = setdiff_a / adj_a / timediff,
         pertime_b = setdiff_b / adj_b / timediff) %>%
  select(tumor_pair_barcode, idh_codel_subtype, prop_a, prop_b, prop_ab, adj_a, adj_b, adj_ab, timediff, pertime_a, pertime_b)

ggdat = df %>%
  select(tumor_pair_barcode, idh_codel_subtype, prop_a, prop_b) %>%
  gather(starts_with("prop"), key = "Variant", value = "Proportion") %>%
  arrange(Proportion)

ggplot(ggdat, aes(x=tumor_pair_barcode, y=Proportion)) +
  geom_bar(aes(fill=Variant), stat="identity") +
  facet_wrap(~idh_codel_subtype, scales = "free_x")

ggdat = df %>%
  select(tumor_pair_barcode, idh_codel_subtype, adj_a, adj_b, adj_ab, timediff) %>%
  gather(starts_with("adj"), key = "Variant", value = "VariantsMB") %>%
  mutate(Variant = factor(Variant, labels=c("Unique in Primary", "Shared", "Unique in Recurrence")))

ggplot(ggdat, aes(x=idh_codel_subtype, y=VariantsMB)) +
  geom_boxplot(aes(color=Variant)) + 
  coord_cartesian(ylim=c(0,10))

ggplot(ggdat, aes(x=timediff, y=VariantsMB)) +
  geom_smooth(aes(color=Variant), method="lm") + 
  facet_wrap(~idh_codel_subtype, scales = "free_x")
  
ggdat = df %>%
  select(tumor_pair_barcode, idh_codel_subtype, pertime_a, pertime_b) %>%
  gather(starts_with("pertime"), key = "ab", value = "MutMo") %>%
  arrange(MutMo)

ggplot(ggdat, aes(x=idh_codel_subtype, y=MutMo)) +
  geom_boxplot(aes(color=ab)) + 
  scale_y_log10()

ggplot(ggdat, aes(x=idh_codel_subtype, y=MutMo)) +
  geom_boxplot(aes(color=ab)) + 
  coord_cartesian(ylim = c(0,50))

##########################

ggplot(ggdat, aes(TP_vaf, R1_vaf)) +
  geom_point(aes(color=case_sex)) + 
  geom_abline(slope=1, alpha=0.2, linetype=2) +
  facet_wrap(~gene_symbol_label) +
  labs(x="Primary", y="First Recurrence", color = "Variant Classification") +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  theme_bw(base_size = 18) +
  theme(axis.text=element_text(size=10))

ggplot(df, aes(TP_vaf,R1_vaf)) +
  geom_point(aes(color=paste(gene_symbol, hgvs_p), shape=dp)) + 
  geom_abline(slope=1, alpha=0.2, linetype=2) +
  labs(x="Primary", y="First Recurrence", color = "Protein Change", shape = "Coverage") +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  theme_bw()

 
########

#   filter(any_called) %>%
#   group_by(sample_barcode,gene_symbol) %>%
#   mutate(n=n()) %>%
#   ungroup()
#   
#   
#     mutate(severity_index = order(severity_score))
#   filter(sample_type %in% c("TP", "R1")) %>%
#   mutate(vaf = alt_count/read_depth, var = sprintf("%s:%s-%s_%s", chrom, start, end, alt)) %>%
#   select(case_barcode, gene_symbol, hgvs_p, var, sample_type, vaf, dp = read_depth) %>%
#   gather(variable, value, vaf, dp) %>%
#   unite(temp, sample_type, variable) %>%
#   spread(temp, value) %>%
#   mutate(dp = factor(case_when(TP_dp > 14 & R1_dp > 14 ~ "DP TP/R1 > 14",
#                         TP_dp > 14 ~ "DP TP>14",
#                         R1_dp > 14 ~ "DP R1>14",
#                         TRUE ~ "DP TP/R1 < 15")))
# 
# #%>%
#  filter(complete.cases(dp))


qres <- dbGetQuery(con, "SELECT DISTINCT ON (1,2,3,4,5) v.chrom, v.start, v.end, v.alt, s.sample_barcode, variant_classification, variant_type, gt.aliquot_barcode, v.hgvs_p, v.hgvs_c, ca.case_barcode, sample_type, gene_symbol, ref_count, alt_count, read_depth, aliquot_portion
  FROM analysis.snvs v, analysis.snv_genotypes gt, biospecimen.aliquots al, biospecimen.samples s, clinical.cases ca
                   WHERE v.chrom = gt.chrom AND v.start = gt.start AND v.end = gt.end AND v.alt = gt.alt AND
                   gt.aliquot_barcode = al.aliquot_barcode AND
                   al.sample_barcode = s.sample_barcode AND
                   ca.case_barcode = s.case_barcode AND
                   al.aliquot_analysis_type = 'WGS' AND
                   v.variant_classification NOT IN ('IGR','Intron') AND
                   v.gene_symbol = 'TERT'
                   ORDER BY 1,2,3,4,5") 

df = qres %>% #filter(aliquot_portion == 1, sample_type != 'NB') %>%
  filter(sample_type %in% c("TP", "R1")) %>%
  mutate(vaf = alt_count/read_depth, var = sprintf("%s:%s-%s_%s", chrom, start, end, alt)) %>%
  select(case_barcode, gene_symbol, chrom, start, end, hgvs_p, hgvs_c, variant_type, variant_classification, var, sample_type, vaf, dp = read_depth) %>%
  gather(variable, value, vaf, dp) %>%
  unite(temp, sample_type, variable) %>%
  spread(temp, value) %>%
  mutate(dp = factor(case_when(TP_dp > 14 & R1_dp > 14 ~ "DP TP/R1 > 14",
                               TP_dp > 14 ~ "DP TP>14",
                               R1_dp > 14 ~ "DP R1>14",
                               TRUE ~ "DP TP/R1 < 15")))

ggplot(subset(df, start>1295200), aes(TP_vaf,R1_vaf)) +
  geom_point(aes(color=factor(start), shape=dp)) + 
  geom_abline(slope=1, alpha=0.2, linetype=2) +
  labs(x="Primary", y="First Recurrence", color = "Protein Change", shape = "Coverage") +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  theme_bw()


#### 


qres <- dbGetQuery(con, "SELECT DISTINCT ON (1,2,3,4,5) v.chrom, v.start, v.\"end\", v.alt, s.sample_barcode, variant_classification, variant_type, gt.aliquot_barcode, v.hgvs_p, ca.case_barcode, sample_type, gene_symbol, ref_count, alt_count, read_depth, aliquot_portion, gt.called
  FROM analysis.snvs v, analysis.snv_genotypes gt, biospecimen.aliquots al, biospecimen.samples s, clinical.cases ca
                   WHERE v.chrom = gt.chrom AND v.start = gt.start AND v.end = gt.end AND v.alt = gt.alt AND
                   variant_classification = 'Missense_Mutation' AND
                   gt.aliquot_barcode = al.aliquot_barcode AND
                   al.sample_barcode = s.sample_barcode AND
                   ca.case_barcode = s.case_barcode AND
                   al.aliquot_analysis_type = 'WXS' AND
                   v.gene_symbol IN ('IDH1','IDH2') AND
                   hgvs_p IN ('p.R132H','p.R132C','p.R132G','p.R132S','p.R172K')
                   ORDER BY 1,2,3,4,5")

df = qres %>% 
  filter(called=="1") %>%
  mutate(vaf = alt_count/read_depth, var = sprintf("%s:%s-%s_%s", chrom, start, end, alt)) %>%
  select(case_barcode, gene_symbol, hgvs_p, var, sample_type, vaf, dp = read_depth) 


ggplot(df, aes(x=sample_type,y=vaf)) + geom_boxplot(aes(color=hgvs_p)) 


###

qres <- dbGetQuery(con, "SELECT DISTINCT ON (1,2,3,4,5) v.chrom, v.start, v.end, v.alt, s.sample_barcode, variant_classification, variant_type, gt.aliquot_barcode, v.hgvs_p, v.hgvs_c, ca.case_barcode, sample_type, gene_symbol, ref_count, alt_count, read_depth, aliquot_portion, gt.called
  FROM analysis.snvs v, analysis.snv_genotypes gt, biospecimen.aliquots al, biospecimen.samples s, clinical.cases ca
                   WHERE v.chrom = gt.chrom AND v.start = gt.start AND v.end = gt.end AND v.alt = gt.alt AND
                   gt.aliquot_barcode = al.aliquot_barcode AND
                   al.sample_barcode = s.sample_barcode AND
                   ca.case_barcode = s.case_barcode AND
                   al.aliquot_analysis_type = 'WXS' AND
                   v.variant_classification NOT IN ('IGR','Intron') AND
                   v.gene_symbol = 'TP53'
                   ORDER BY 1,2,3,4,5") 

df = qres %>% 
  filter(called=="1") %>%
  mutate(vaf = alt_count/read_depth, var = sprintf("%s:%s-%s_%s", chrom, start, end, alt)) %>%
  select(case_barcode, gene_symbol, hgvs_p, var, variant_classification, variant_type, sample_type, vaf, dp = read_depth) 


ggplot(df, aes(x=sample_type,y=vaf)) + geom_boxplot(aes(color=variant_type)) 
