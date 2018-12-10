## push blocklist to db

library(tidyverse)
library(DBI)
library(odbc)
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

tmp = read.delim("data/ref/block_review_allow_quality_lists_20181206.txt", as.is = TRUE)

df <- tmp %>%
  select(aliquot_barcode,
         fingerprint_exclusion,
         coverage_exclusion = coverage_mut_exclusion,
         cnv_exclusion = manual_cn_exlusion,
         clinical_exclusion = surgical_interval_exclusion,
         fingerprint_exclusion_reason,
         coverage_exclusion_reason = mut_exclusion_reason,
         cnv_exclusion_reason = cn_exclusion_reason,
         clinical_exclusion_reason = surgical_interval_exclusion_reason) %>%
  mutate(clinical_exclusion_reason = ifelse(clinical_exclusion_reason == "", NA, clinical_exclusion_reason))

dbWriteTable(con, Id(schema="analysis",table="blocklist"), df, append = FALSE)
