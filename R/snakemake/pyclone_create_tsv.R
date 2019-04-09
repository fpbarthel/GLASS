library(tidyverse)
library(DBI)

.libPaths("/home/barthf/R/x86_64-pc-linux-gnu-library/3.3")
## database connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv2")

## input/output parameters
barcode   <- snakemake@wildcards[["aliquot_barcode"]]
tsv       <- snakemake@output[["tsv"]]

## Logging
message("Processing ", barcode)

## process parameters
case_barcode  <- substring(barcode,1,12)

## Fetch data from DB
rs <- dbSendQuery(con,read_file("sql/pyclone/pyclone_create_tsv.sql"))
dbBind(rs, list(case_barcode))
qres <- dbFetch(rs)
  
df <- qres %>%
  filter(aliquot_barcode == barcode) %>% 
  select(mutation_id,ref_counts,var_counts,normal_cn,minor_cn,major_cn)
  
write.table(df, file = tsv, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)