library(tidyverse)
library(odbc)
library(DBI)
library(ggthemes)
library(ggplot2)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

q = "SELECT mf.aliquot_barcode, s.case_barcode, treatment_tmz, who_classification, ts.sample_type, ta.aliquot_portion AS portion, ta.aliquot_analysis_type AS analysis_type,
surgery_number, surgical_interval_mo, histology, grade, idh_status, codel_status, coverage_adj_mut_freq, cumulative_coverage
FROM analysis.mutation_freq mf
INNER JOIN biospecimen.aliquots ta ON ta.aliquot_barcode = mf.aliquot_barcode
INNER JOIN biospecimen.samples ts ON ts.sample_barcode = ta.sample_barcode
LEFT JOIN clinical.surgeries s ON s.sample_barcode = ts.sample_barcode"

qres <- dbGetQuery(con, q)

df <- qres %>% 
  filter(case_barcode != "GLSS-MD-0018",
         sample_type %in% c("TP","R1"),
         cumulative_coverage > 60000000,
         complete.cases(histology, grade, coverage_adj_mut_freq, surgical_interval_mo, idh_status, codel_status, grade, histology, who_classification)) %>%
  group_by(case_barcode, surgery_number) %>%
  summarize(coverage_adj_mut_freq = mean(coverage_adj_mut_freq),
            surgical_interval_mo = surgical_interval_mo[1],
            idh_status = idh_status[1],
            codel_status = codel_status[1],
            histology = histology[1],
            grade = grade[1],
            treatment_tmz = treatment_tmz[1], 
            who_classification = who_classification[1]) %>%
  ungroup() %>%
  group_by(case_barcode) %>%
  mutate(idh_codel_grp = ifelse(any(idh_status == "IDH.mt") && any(codel_status == "codel"), "IDHmut-codel", 
                                ifelse(any(idh_status == "IDH.mt"), "IDHmut",
                                       ifelse(any(idh_status == "IDH.wt"), "IDHwt", NA))),
         hypermut_grp = ifelse(any(coverage_adj_mut_freq > 8), "hypermut", "non-hypermut"),
         n = n()) %>%
  ungroup() %>% filter(n>1)

df_hyper = df %>% filter(hypermut_grp == "hypermut", complete.cases(idh_codel_grp))
df_nonhyper = df %>% filter(hypermut_grp == "non-hypermut", complete.cases(idh_codel_grp))

ggplot(df_hyper, aes(x=surgical_interval_mo, y=coverage_adj_mut_freq)) + 
  geom_point(aes(shape=who_classification)) + geom_line(aes(group=case_barcode, color=treatment_tmz)) + 
  facet_wrap(~ idh_codel_grp) + scale_y_log10()

ggplot(df_nonhyper, aes(x=surgical_interval_mo, y=coverage_adj_mut_freq)) + 
  geom_point(aes(shape=who_classification)) + geom_line(aes(group=case_barcode, color=treatment_tmz)) +
  facet_wrap(~ idh_codel_grp)

ggplot(df, aes(x=surgical_interval_mo, y=coverage_adj_mut_freq)) + 
  geom_smooth(method="lm", aes(color = hypermut_grp)) +
  scale_y_log10() +
  facet_wrap(~ idh_codel_grp, scales="free_x")

ggplot(df_hyper, aes(x=surgical_interval_mo, y=coverage_adj_mut_freq)) + 
  geom_smooth(method="lm") +
  facet_wrap(~ idh_codel_grp, scales="free_x")

library(nlme)
lmList(coverage_adj_mut_freq~surgical_interval_mo|idh_codel_grp, data = df_hyper[,c("surgical_interval_mo","coverage_adj_mut_freq","idh_codel_grp")])
lmList(coverage_adj_mut_freq~surgical_interval_mo|idh_codel_grp, data = df_nonhyper[,c("surgical_interval_mo","coverage_adj_mut_freq","idh_codel_grp")])

ggplot(df_nonhyper, aes(x=surgical_interval_mo, y=coverage_adj_mut_freq)) + 
  geom_smooth(method="lm") +
  facet_wrap(~ idh_codel_grp, scales="free_x")
  
  geom_point(aes(shape=who_classification)) + geom_line(aes(group=case_barcode, color=treatment_tmz)) +
  facet_wrap(~ idh_codel_grp)

df_idhmut_codel = df %>% filter(idh_codel_grp == "IDHmut-codel")
df_idhmut = df %>% filter(idh_codel_grp == "IDHmut")
df_idhwt = df %>% filter(idh_codel_grp == "IDHwt")

ggplot(df_idhwt, aes(x = surgical_interval_mo, y=case_barcode)) +
  geom_line() + geom_point(aes(size=coverage_adj_mut_freq, color = who_classification))

ggplot(df_idhmut, aes(x = surgical_interval_mo, y=case_barcode)) +
  geom_line() + geom_point(aes(size=coverage_adj_mut_freq, color = who_classification))


ggplot(df_idhmut_codel, aes(x = surgical_interval_mo, y=case_barcode)) +
  geom_line() + geom_point(aes(size=coverage_adj_mut_freq, color = who_classification))


