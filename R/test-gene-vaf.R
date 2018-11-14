library(tidyverse)
library(odbc)
library(DBI)
library(ggthemes)
library(ggplot2)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

# qres <- dbGetQuery(con, "SELECT DISTINCT ON (1,2,3,4,5) v.chrom, v.start, v.\"end\", v.alt, s.sample_barcode, variant_classification, variant_type, gt.aliquot_barcode, v.hgvs_p, ca.case_barcode, sample_type, gene_symbol, ref_count, alt_count, read_depth, aliquot_portion
#   FROM analysis.snvs v, analysis.snv_genotypes gt, biospecimen.aliquots al, biospecimen.samples s, clinical.cases ca
#   WHERE v.chrom = gt.chrom AND v.start = gt.start AND v.end = gt.end AND v.alt = gt.alt AND
#   variant_classification = 'Missense_Mutation' AND
#   gt.aliquot_barcode = al.aliquot_barcode AND
#   al.sample_barcode = s.sample_barcode AND
#   ca.case_barcode = s.case_barcode AND
#   al.aliquot_analysis_type = 'WXS' AND
#   v.gene_symbol IN ('IDH1','IDH2') AND
#   hgvs_p IN ('p.R132H','p.R132C','p.R132G','p.R132S','p.R172K')
#   ORDER BY 1,2,3,4,5") 

q = "SELECT ts.case_barcode, case_sex, s.idh_status, s.codel_status, v.gene_symbol, ts.sample_type, ts.sample_barcode, 
v.chrom, v.start, v.end, v.alt, v.variant_classification, v.variant_type, 
gt.aliquot_barcode, v.hgvs_p, gt.ref_count, gt.alt_count, gt.read_depth, sift, polyphen, gt.called, mf.coverage_adj_mut_freq
FROM analysis.snvs v
FULL JOIN analysis.snv_genotypes gt ON v.chrom = gt.chrom AND v.start = gt.start AND v.end = gt.end AND v.alt = gt.alt
INNER JOIN biospecimen.aliquots ta ON ta.aliquot_barcode = gt.aliquot_barcode
INNER JOIN biospecimen.samples ts ON ts.sample_barcode = ta.sample_barcode
LEFT JOIN analysis.mutation_freq mf ON mf.aliquot_barcode = gt.aliquot_barcode
LEFT JOIN clinical.surgeries s ON s.sample_barcode = ts.sample_barcode
INNER JOIN clinical.cases ca ON ts.case_barcode = ca.case_barcode
WHERE ts.sample_type IN ('TP','R1') AND ((v.variant_classification IN ('In_Frame_Ins','In_Frame_Del','Frame_Shift_Del','Frame_Shift_Ins','Missense_Mutation','Nonstop_Mutation','Nonsense_Mutation')
AND v.gene_symbol IN ('TP53','ATRX','IDH1','EGFR','PTEN','NF1','CIC','FUBP1','NOTCH1','PIK3CA','PIK3R1')) 
OR (v.variant_classification = '5''Flank' AND v.gene_symbol = 'TERT'))"

qres <- dbGetQuery(con, q)

df = qres %>% 
  filter(coverage_adj_mut_freq < 8, complete.cases(idh_status, codel_status)) %>% ## (1) Filter out hypermutator samples
  group_by(case_barcode) %>%
  mutate(idh_codel_grp = ifelse(any(idh_status == "IDH.mt") && any(codel_status == "codel"), "IDHmut-codel", 
                                ifelse(any(idh_status == "IDH.mt"), "IDHmut",
                                       ifelse(any(idh_status == "IDH.wt"), "IDHwt", NA)))) %>%
  ungroup() %>%
  mutate(var = sprintf("%s:%s-%s_%s", chrom, start, end, alt),
         called = called == "1",
         severity_score = case_when(variant_classification == "Nonsense" ~ 0,
                                    variant_classification %in% c("Frame_Shift_Del","Frame_Shift_Del") ~ 1,
                                    variant_classification %in% c("In_Frame_Del","In_Frame_Ins") ~ 2,
                                    variant_classification == "Missense_Mutation" ~ 3,
                                    variant_classification == "5'Flank" ~ 4)) %>%
  group_by(sample_barcode, var) %>%
  mutate(optimal_variant = order(read_depth, decreasing = T)) %>%
  ungroup() %>%
  filter(optimal_variant == 1) %>% ## (2) For each sample/variant combination, select the variant with the highest read depth, eg. when a variant was profiled across multple sectors or with both WGS and WES
  group_by(case_barcode, var) %>%
  mutate(sufficient_dp = all(ref_count + alt_count > 14)) %>%
  ungroup() %>%
  filter(sufficient_dp) %>% ## (3) For each patient/variant combination filter out variants that do not have >14x coverage across all subsamples
  group_by(case_barcode,gene_symbol) %>%
  mutate(any_called = any(called, na.rm=T),
         num_samples = n_distinct(sample_barcode)) %>%
  ungroup() %>%
  group_by(gene_symbol) %>%
  mutate(num_patient = n_distinct(case_barcode),
         gene_symbol_label = sprintf("%s (n=%s)", gene_symbol, num_patient)) %>%
  ungroup() %>%
  filter(any_called, num_samples > 1) %>% ## (4) a. filter out gene/patient combinations where that gene was not called as mutant in any subsample and b. filter out singletons (unpaired samples)
  arrange(severity_score) %>%
  group_by(case_barcode, gene_symbol) %>%
  mutate(keep = chrom == chrom[which(called)[1]] & 
           start == start[which(called)[1]] & 
           end == end[which(called)[1]] & 
           alt == alt[which(called)[1]]) %>%
  ungroup() %>%
  filter(keep) ## (5) For each patient/gene combination, keep only those specific variants that were called in at least one subsample

ggdat = df %>%
  mutate(vaf = alt_count/(alt_count+ref_count)) %>%
  select(case_barcode, idh_codel_grp, case_sex, gene_symbol_label, chrom, start, end, hgvs_p, variant_type, variant_classification, var, sample_type, vaf, dp = read_depth) %>%
  gather(variable, value, vaf, dp) %>%
  unite(temp, sample_type, variable) %>%
  spread(temp, value) 

ggplot(ggdat, aes(TP_vaf, R1_vaf)) +
  geom_point(aes(color=variant_classification, shape=idh_codel_grp)) + 
  geom_abline(slope=1, alpha=0.2, linetype=2) +
  facet_wrap(~gene_symbol_label) +
  labs(x="Primary", y="First Recurrence", color = "Variant Classification", shape = "Glioma subtype") +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  theme_bw(base_size = 18) +
  theme(axis.text=element_text(size=10))


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
