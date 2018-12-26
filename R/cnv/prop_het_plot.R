library(DBI)
library(tidyverse)
library(ggplot2)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")
q <- "SELECT ts.*, cl.idh_codel_subtype
FROM analysis.titan_seg_paired_comparison ts
LEFT JOIN biospecimen.aliquots al ON al.aliquot_barcode = ts.tumor_barcode_a
LEFT JOIN clinical.surgeries cl ON cl.sample_barcode = al.sample_barcode"

qres <- dbGetQuery(con, q)

ggplot(qres, aes(x=delta_prop_het, color = idh_codel_subtype)) + geom_density() + coord_cartesian(xlim=c(-1,1)) + labs(x="Heterozygous proportion of the genome\n(recurrence-primary)")
ggplot(qres, aes(x=prop_delta_eq, color = idh_codel_subtype)) + geom_density() + labs(x="Proportion of the genome with identical copy number states\n(recurrence-primary)")

+ geom_smooth(method = "lm") + facet_wrap(~mutation_status)
ggplot(qres, aes(x=case_age_diagnosis_years, y = relative_contribution)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~mutation_status + idh_codel_subtype)

cor.test(~ relative_contribution + case_age_diagnosis_years, data = qres, subset = qres$mutation_status=="shared")
cor.test(~ relative_contribution + case_age_diagnosis_years, data = qres, subset = qres$mutation_status=="primary")
cor.test(~ relative_contribution + case_age_diagnosis_years, data = qres, subset = qres$mutation_status=="recurrent")

ggplot(qres, aes(x=surgical_interval_mo, y = relative_contribution)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  facet_wrap(~mutation_status + idh_codel_subtype, scales = "free_x")
       