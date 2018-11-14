library(tidyverse)
library(odbc)
library(DBI)
library(ggthemes)
library(ggplot2)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

q <- "SELECT mf.aliquot_barcode, s.case_barcode, s.sample_barcode, sample_type, idh_status, codel_status, treatment_tmz, treatment_tmz_cycles_6, coverage_adj_mut_freq, cumulative_coverage, aliquot_analysis_type
FROM analysis.mutation_freq mf
INNER JOIN biospecimen.aliquots al ON al.aliquot_barcode = mf.aliquot_barcode
INNER JOIN biospecimen.samples s ON s.sample_barcode = al.sample_barcode
INNER JOIN clinical.surgeries su ON s.sample_barcode = su.sample_barcode
INNER JOIN clinical.cases ca ON s.case_barcode = ca.case_barcode"

qres <- dbGetQuery(con, q)

ggplot(qres, aes(x=coverage_adj_mut_freq)) + geom_density(aes(color=aliquot_analysis_type)) +
  theme_bw()

ggplot(qres, aes(x=coverage_adj_mut_freq)) + geom_density(aes(color=aliquot_analysis_type)) +
  theme_bw() + coord_cartesian(xlim = c(0,50))

ggplot(qres, aes(y=coverage_adj_mut_freq)) + geom_boxplot(aes(x=sample_type, color=sample_type)) +
  theme_bw() + scale_y_log10()

df <- qres %>%
  filter(sample_type %in% c("TP","R1"),
         cumulative_coverage > 60000000) %>%
  group_by(case_barcode) %>%
  mutate(num_samples = n_distinct(sample_barcode),
         tmz = case_when("TP" %in% sample_type && treatment_tmz[sample_type=="TP"]=="1" ~ "TMZ post-primary",
                         "TP" %in% sample_type && treatment_tmz[sample_type=="TP"]=="0" ~ "No TMZ",
                         TRUE ~ NA_character_),
         idh_codel_grp = ifelse(any(idh_status == "IDH.mt") && any(codel_status == "codel"), "IDHmut-codel", 
                                ifelse(any(idh_status == "IDH.mt"), "IDHmut",
                                       ifelse(any(idh_status == "IDH.wt"), "IDHwt", NA)))) %>%
  
  ungroup() %>%
  filter(num_samples > 1) %>% ## (1) Drop singletons
  group_by(sample_barcode) %>%
  mutate(optimal_coverage = order(cumulative_coverage, decreasing = T)) %>%
  ungroup() %>%
  filter(optimal_coverage == 1) %>% ## (2) For each sample select mutation rate from aliquot with highest coverage
  select(case_barcode,sample_type,tmz,idh_codel_grp,coverage_adj_mut_freq) %>%
  spread(sample_type,coverage_adj_mut_freq)
  
gg <- ggplot(df, aes(x = TP, y = R1)) + 
  geom_vline(xintercept = 10, linetype=2, alpha=0.6) +
  geom_hline(yintercept = 10, linetype=2, alpha=0.6) +
  geom_abline(slope=1, linetype=2, alpha=0.6) +
  geom_point(aes(color = tmz)) +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim = c(0.1,100), ylim = c(0.1,100)) +
  labs(title = "Coverage adjusted mutation frequency", #expression(paste("Coverage adjusted mutation frequency \n", italic("mutations/covered MB"))),
       x="Primary", y="First Recurrence", color = "TMZ treatment") +
  theme_bw(base_size = 18) +
  theme(axis.text=element_text(size=10))
  
gg
gg + facet_wrap(~idh_codel_grp)

df = qres %>% filter(aliquot_analysis_type == "WXS") %>%
  mutate(sample_type = case_when(sample_type == "TP" ~ 0,
                                 sample_type == "R1" ~ 1,
                                 sample_type == "R2" ~ 2,
                                 sample_type == "R3" ~ 3))

ggplot(df, aes(y=coverage_adj_mut_freq)) + geom_line(aes(x=sample_type, color=case_barcode)) +
  guides(color=FALSE) + theme_bw()


ggplot(df, aes(y=coverage_adj_mut_freq)) + geom_line(aes(x=sample_type, color=case_barcode)) +
  guides(color=FALSE) + theme_bw() + coord_cartesian(ylim = c(0,50))
