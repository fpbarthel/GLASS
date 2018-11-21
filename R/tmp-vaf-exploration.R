#######################################################
# Useful commands compiled by Floris for interacting with the database.
# Date: 2018.11.20
# Author: Kevin J.
#######################################################

# Essential packages:
library(tidyverse)
library(DBI)

#######################################################
# Make connection with the database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

# Query the database to retrieve specific variants.
qres <- dbGetQuery(con, "SELECT mf.aliquot_barcode, case_barcode, sample_type, cumulative_coverage, coverage_adj_mut_freq, aliquot_analysis_type
                   FROM analysis.mutation_freq mf
                   LEFT JOIN biospecimen.aliquots al ON al.aliquot_barcode = mf.aliquot_barcode
                   LEFT JOIN biospecimen.samples s ON s.sample_barcode = al.sample_barcode")


qres2 = dbGetQuery(con, "SELECT *FROM analysis.called_genotypes gt
           INNER JOIN biospecimen.aliquots al ON al.aliquot_barcode = gt.aliquot_barcode
           INNER JOIN biospecimen.samples s ON s.sample_barcode = al.sample_barcode
           INNER JOIN clinical.cases ca ON ca.case_barcode = s.case_barcode
           WHERE ca.case_barcode = 'GLSS-MD-0002'")

# Merge mutation frequencies:
qres3 = dbGetQuery(con, "SELECT tmc.*, sa.sample_type AS sample_type_a, sb.sample_type AS sample_type_b, ua.surgery_number AS surgery_number_a, ub.surgery_number AS surgery_number_b, ua.surgical_interval_mo AS surgical_interval_a, ub.surgical_interval_mo AS surgical_interval_b, ua.idh_codel_subtype,
ma.cumulative_coverage AS coverage_a, mb.cumulative_coverage AS coverage_b, ma.coverage_adj_mut_freq AS mf_a, mb.coverage_adj_mut_freq AS mf_b
FROM analysis.tumor_mut_comparison tmc
INNER JOIN biospecimen.aliquots aa ON tmc.tumor_a_barcode = aa.aliquot_barcode
INNER JOIN biospecimen.aliquots ab ON tmc.tumor_b_barcode = ab.aliquot_barcode
INNER JOIN biospecimen.samples sa ON aa.sample_barcode = sa.sample_barcode
INNER JOIN biospecimen.samples sb ON ab.sample_barcode = sb.sample_barcode
LEFT JOIN clinical.surgeries ua ON ua.sample_barcode = sa.sample_barcode
LEFT JOIN clinical.surgeries ub ON ub.sample_barcode = sb.sample_barcode
LEFT JOIN analysis.mutation_freq ma ON ma.aliquot_barcode = aa.aliquot_barcode
LEFT JOIN analysis.mutation_freq mb ON mb.aliquot_barcode = ab.aliquot_barcode")


plot(qres3$surgical_interval_b, qres3$intersection_ab)

tmp = qres3 %>% mutate_if(bit64::is.integer64, as.double)

ggplot(tmp, aes(x=surgical_interval_b, y=intersection_ab)) + geom_point() + facet_grid(~idh_codel_subtype)
