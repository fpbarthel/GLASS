library(tidyverse)
library(odbc)
library(DBI)
library(ggthemes)
library(ggplot2)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

qres <- dbGetQuery(con, "SELECT DISTINCT ON (1,2,3,4,5) v.chrom, v.start, v.\"end\", v.alt, s.sample_barcode, variant_classification, variant_type, gt.aliquot_barcode, v.hgvs_p, ca.case_barcode, sample_type, gene_symbol, ref_count, alt_count, read_depth, aliquot_portion
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

df = qres %>% #filter(aliquot_portion == 1, sample_type != 'NB') %>%
  filter(sample_type %in% c("TP", "R1")) %>%
  mutate(vaf = alt_count/read_depth, var = sprintf("%s:%s-%s_%s", chrom, start, end, alt)) %>%
  select(case_barcode, gene_symbol, hgvs_p, var, sample_type, vaf, dp = read_depth) %>%
  gather(variable, value, vaf, dp) %>%
  unite(temp, sample_type, variable) %>%
  spread(temp, value) %>%
  mutate(dp = factor(case_when(TP_dp > 14 & R1_dp > 14 ~ "DP TP/R1 > 14",
                        TP_dp > 14 ~ "DP TP>14",
                        R1_dp > 14 ~ "DP R1>14",
                        TRUE ~ "DP TP/R1 < 15")))

#%>%
#  filter(complete.cases(dp))

ggplot(df, aes(TP_vaf,R1_vaf)) +
  geom_point(aes(color=paste(gene_symbol, hgvs_p), shape=dp)) + 
  geom_abline(slope=1, alpha=0.2, linetype=2) +
  labs(x="Primary", y="First Recurrence", color = "Protein Change", shape = "Coverage") +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  theme_bw()

 
########

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
