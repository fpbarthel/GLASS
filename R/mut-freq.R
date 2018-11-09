library(tidyverse)
library(odbc)
library(DBI)
library(ggthemes)
library(ggplot2)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

qres <- dbGetQuery(con, "SELECT mf.aliquot_barcode, case_barcode, sample_type, coverage_adj_mut_freq, aliquot_analysis_type
                   FROM analysis.mutation_freq mf
                   LEFT JOIN biospecimen.aliquots al ON al.aliquot_barcode = mf.aliquot_barcode
                   LEFT JOIN biospecimen.samples s ON s.sample_barcode = al.sample_barcode")

ggplot(qres, aes(x=coverage_adj_mut_freq)) + geom_density(aes(color=aliquot_analysis_type)) +
  theme_bw()

ggplot(qres, aes(x=coverage_adj_mut_freq)) + geom_density(aes(color=aliquot_analysis_type)) +
  theme_bw() + coord_cartesian(xlim = c(0,50))

ggplot(qres, aes(y=coverage_adj_mut_freq)) + geom_boxplot(aes(x=sample_type, color=sample_type)) +
  theme_bw() + scale_y_log10()

df = qres %>% filter(aliquot_analysis_type == "WXS") %>%
  mutate(sample_type = case_when(sample_type == "TP" ~ 0,
                                 sample_type == "R1" ~ 1,
                                 sample_type == "R2" ~ 2,
                                 sample_type == "R3" ~ 3))

ggplot(df, aes(y=coverage_adj_mut_freq)) + geom_line(aes(x=sample_type, color=case_barcode)) +
  guides(color=FALSE) + theme_bw()


ggplot(df, aes(y=coverage_adj_mut_freq)) + geom_line(aes(x=sample_type, color=case_barcode)) +
  guides(color=FALSE) + theme_bw() + coord_cartesian(ylim = c(0,50))
