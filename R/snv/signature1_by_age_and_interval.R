library(DBI)
library(tidyverse)
library(ggplot2)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")
q <- "SELECT ms.tumor_pair_barcode, tp.case_barcode, mutation_status, ms.relative_contribution, cs.case_age_diagnosis_years, tp.surgical_interval_mo, su.idh_codel_subtype
FROM analysis.mutsig_private_vs_shared ms
LEFT JOIN analysis.tumor_pairs tp ON tp.tumor_pair_barcode = ms.tumor_pair_barcode
LEFT JOIN clinical.cases cs ON cs.case_barcode = tp.case_barcode
LEFT JOIN clinical.surgeries su ON su.sample_barcode = substring(tp.tumor_pair_barcode from 1 for 15)
WHERE signature = 'Signature.1' AND mut_count >= 100"

qres <- dbGetQuery(con, q)

ggplot(qres, aes(x=case_age_diagnosis_years, y = relative_contribution)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~mutation_status)
ggplot(qres, aes(x=case_age_diagnosis_years, y = relative_contribution)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~mutation_status + idh_codel_subtype)

cor.test(~ relative_contribution + case_age_diagnosis_years, data = qres, subset = qres$mutation_status=="shared")
cor.test(~ relative_contribution + case_age_diagnosis_years, data = qres, subset = qres$mutation_status=="primary")
cor.test(~ relative_contribution + case_age_diagnosis_years, data = qres, subset = qres$mutation_status=="recurrent")

ggplot(qres, aes(x=surgical_interval_mo, y = relative_contribution)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  facet_wrap(~mutation_status + idh_codel_subtype, scales = "free_x")
       