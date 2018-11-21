
library(devtools); install_github("im3sanger/dndscv")
library(dndscv)
library(tidyverse)
library(odbc)
library(DBI)
library(ggthemes)
library(ggplot2)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

qres <- dbGetQuery(con, "SELECT case_barcode, chrom, start, ref, alt
                   FROM analysis.called_genotypes gt
                   LEFT JOIN biospecimen.aliquots al ON al.aliquot_barcode = gt.aliquot_barcode
                   LEFT JOIN biospecimen.samples s ON al.sample_barcode = s.sample_barcode
                   WHERE variant_classification NOT IN ('Intron','IGR') AND s.sample_type = 'TP'
                   ORDER BY chrom,start desc")


###

df = data.frame(sampleID = c("a","b","c"),
                chr = rep(10,3),
                pos = c(115312925,114910854,112838833),
                ref = c("T","C","C"),
                mut = c("C","T","G"))

df <- qres %>% filter(chrom %in% as.character(1:22)) %>%
  mutate(chrom = as.numeric(chrom)) %>%
  arrange(chrom)

dnds = dndscv(df, refdb = "hg19")

dnds

## END ##