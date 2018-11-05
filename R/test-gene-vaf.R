library(tidyverse)
library(odbc)
library(DBI)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

df <- dbGetQuery(con, "SELECT DISTINCT ON (1,2,3,4,5) v.chrom, v.start, v.\"end\", v.alt, s.sample_barcode, variant_classification, variant_type, gt.aliquot_barcode, ca.case_barcode, sample_type, gene_symbol, ref_count, alt_count, read_depth, aliquot_portion
  FROM analysis.snvs v, analysis.snv_genotypes gt, biospecimen.aliquots al, biospecimen.samples s, clinical.cases ca
  WHERE v.chrom = gt.chrom AND v.start = gt.start AND v.end = gt.end AND v.alt = gt.alt AND
  variant_classification = 'Missense_Mutation' AND
  gt.aliquot_barcode = al.aliquot_barcode AND
  al.sample_barcode = s.sample_barcode AND
  ca.case_barcode = s.case_barcode AND
  al.aliquot_analysis_type = 'WGS' AND
  v.gene_symbol = 'IDH1'
  ORDER BY 1,2,3,4,5") %>%
  filter(aliquot_portion == 1, sample_type != 'NB') %>%
  mutate(vaf = alt_count/read_depth, var = sprintf("%s:%s-%s_%s", chrom, start, end, alt)) %>%
  select(case_barcode, var, sample_type, vaf) %>%
  spread(sample_type, vaf)

ggplot(df, aes(TP,R1)) + geom_point()

 
