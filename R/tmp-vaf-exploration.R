#######################################################
# Examine the mutation data for the 
# Date: 2018.11.08
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

## Create a stacked barplot for `GLSS-MD-0002` for shared and private mutations.
colnames(qres2)

qres2 %>% 
  group_by(sample_type) %>% 
  summarise(mutation_count = tally())
