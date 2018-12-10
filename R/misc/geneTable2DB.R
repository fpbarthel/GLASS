## export gene table to db

library(tidyverse)
library(DBI)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

setwd("/Volumes/Helix-Common/GLASS-analysis/")

generef = read.delim(file = "data/ref/ncbiRefSeqCurated_hg19.tsv", as.is = TRUE)
genecov = read_tsv(file = "data/ref/gene.covariates.txt")

df = generef %>% 
  group_by(name2) %>%
  filter(row_number()==1) %>%
  ungroup() %>%
  transmute(gene_symbol = name2,
            transcript_id = name,
            chrom = gsub("chr","",chrom),
            pos = sprintf("[%s,%s)",txStart,txEnd),
            strand = strand,
            exons = exonCount,
            tx_size = txEnd - txStart,
            cds_size = sapply(mapply("-", lapply(str_split(exonEnds, ","), as.numeric), lapply(str_split(exonStarts, ","), as.numeric), SIMPLIFY = FALSE), sum, na.rm = TRUE) ) %>%
  filter(chrom %in% c(1:22,"X","Y")) %>%
  left_join(genecov, by = c("gene_symbol" = "gene")) %>%
  mutate(expr = ifelse(is.nan(expr), 0, expr),
         reptime = ifelse(is.nan(reptime), 0, reptime),
         hic = ifelse(is.nan(hic), 0, hic)) %>%
  arrange(gene_symbol)

dbWriteTable(con, Id(schema="ref",table="genes"), df, append=T)
