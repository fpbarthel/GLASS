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
#qres3 <- dbGetQuery(con, 'SELECT sa.sample_type AS sample_type_a, sb.sample_type AS sample_type_b, ua.surgery_number AS surgery_number_a, ub.surgery_number AS surgery_number_b, ua.surgical_interval_mo AS surgical_interval_a, ub.surgical_interval_mo AS surgical_interval_b, ua.idh_codel_subtype,ma.cumulative_coverage AS coverage_a, mb.cumulative_coverage AS coverage_b, ma.coverage_adj_mut_freq AS mf_a, mb.coverage_adj_mut_freq AS mf_b FROM analysis.tumor_mut_comparison tmc INNER JOIN biospecimen.aliquots aa ON tmc.tumor_a_barcode = aa.aliquot_barcode INNER JOIN biospecimen.aliquots ab ON tmc.tumor_b_barcode = ab.aliquot_barcode INNER JOIN biospecimen.samples sa ON aa.sample_barcode = sa.sample_barcode INNER JOIN biospecimen.samples sb ON ab.sample_barcode = sb.sample_barcode LEFT JOIN clinical.surgeries ua ON ua.sample_barcode = sa.sample_barcode LEFT JOIN clinical.surgeries ub ON ub.sample_barcode = sb.sample_barcode LEFT JOIN analysis.mutation_freq ma ON ma.aliquot_barcode = aa.aliquot_barcode LEFT JOIN analysis.mutation_freq mb ON mb.aliquot_barcode = ab.aliquot_barcode')

#plot(qres3$surgical_interval_b, qres3$intersection_ab)

#tmp = qres3 %>% mutate_if(bit64::is.integer64, as.double)

#ggplot(tmp, aes(x=surgical_interval_b, y=intersection_ab)) + geom_point() + facet_grid(~idh_codel_subtype)

qres3 <- dbGetQuery(con, 'SELECT * FROM analysis.tumor_mut_comparison_anno')
df <- qres3 %>% 
  mutate_if(bit64::is.integer64, as.double) %>%
  gather(starts_with("adj"), key = "VariantType", value = "VariantsMB") %>%
  mutate(VariantType = factor(VariantType,
                              levels = c("adj_a","adj_b","adj_intersection","adj_union"),
                              labels = c("Unique to Primary", "Unique To Recurrence", "Shared", "Union")))

ggplot(df, aes(y = VariantsMB, fill = VariantType)) +
  geom_boxplot() + 
  scale_y_log10() +
  facet_wrap(~idh_codel_subtype) + 
  labs(y = "Variants at 14x coverage / megabase > 14x coverage") + 
  theme_bw(base_size = 18)

g <- ggplot(df, aes(y = VariantsMB, x = surgical_interval_b - surgical_interval_a, color = VariantType)) + 
  geom_smooth(method = "lm") +
  geom_point(alpha=0.6) + 
  scale_y_log10() +
  labs(y = "Variants at 14x coverage / megabase > 14x coverage", x = "Surgical Interval (months)") + 
  theme_bw(base_size = 18)

g
g + facet_wrap(~idh_codel_subtype, scales = "free_x") 

df2 <- df %>%
  filter(mf_a < 10, mf_b < 10)

ggplot(df2, aes(y = VariantsMB, fill = VariantType)) +
  geom_boxplot() + 
  facet_wrap(~idh_codel_subtype) + 
  labs(y = "Variants at 14x coverage / megabase > 14x coverage") + 
  theme_bw(base_size = 18)

g <- ggplot(df2, aes(y = VariantsMB, x = surgical_interval_b - surgical_interval_a, color = VariantType)) + 
  geom_smooth(method = "lm") +
  geom_point(alpha=0.6) + 
  labs(y = "Variants at 14x coverage / megabase > 14x coverage", x = "Surgical Interval (months)") + 
  theme_bw(base_size = 18)

g + coord_cartesian(ylim = c(0,7.5))
g + facet_wrap(~idh_codel_subtype, scales = "free_x") + coord_cartesian(ylim = c(0,7.5))
