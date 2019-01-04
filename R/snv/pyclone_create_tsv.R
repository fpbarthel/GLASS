library(tidyverse)
#.libPaths("/home/barthf/R/x86_64-pc-linux-gnu-library/3.3")
library(DBI)

## database connection
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

## input/output parameters
barcode   <- snakemake@wildcards[["aliquot_barcode"]]
tsv       <- snakemake@output[["tsv"]]

## Logging
message("Processing ", barcode)

## process parameters
case_barcode  <- substring(barcode,1,12)
analysis_type <- substring(barcode,21,23)

## Fetch data from DB
rs <- dbSendQuery(con,read_file("sql/pyclone_create_tsv.sql"))
dbBind(rs, list(case_barcode, analysis_type))
qres <- dbFetch(rs)
  
df <- qres %>%
  filter(aliquot_barcode == barcode) %>% 
  mutate(has_30x_in_all_subsamples = num_aliquots_variant_30x == num_aliquots_total,
         has_50x_in_all_subsamples = num_aliquots_variant_50x == num_aliquots_total,
         non_null_cn_in_all_subsamples = num_aliquots_non_null_cn == num_aliquots_total) %>%
  filter(has_30x_in_all_subsamples, non_null_cn_in_all_subsamples, complete.cases(ref_counts,var_counts,normal_cn,minor_cn,major_cn), major_cn > 0) %>%
  select(mutation_id,ref_counts,var_counts,normal_cn,minor_cn,major_cn)
  
write.table(df, file = tsv, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)