## export gene table to db

library(tidyverse)
library(DBI)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

setwd("/Volumes/Helix-Common/GLASS-analysis/")

cbref = read.delim(file = "data/ref/human_grch37_hg19_ucsc_cytoBand.txt", as.is = TRUE, header = FALSE)
colnames(cbref) = c("chrom", "start", "end", "cytoband", "giestain")
#genecov = read_tsv(file = "data/ref/gene.covariates.txt")

df <- cbref %>% 
  transmute(cytoband = cytoband,
            chrom = gsub("chr","",chrom),
            pos = sprintf("[%s,%s)",start,end),
            gie_stain = giestain)
  
dbWriteTable(con, Id(schema="ref",table="cytobands"), df, append=T)
