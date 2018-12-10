### push titan seg into db

library(tidyverse)
library(DBI)
library(odbc)
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

segfiles <- list.files('results/cnv/titanfinal/seg', full.names = TRUE)

lapply(segfiles, function(f){
  message(f)
  dat <- read.delim(f, as.is=T, header=T, row.names = NULL)
  df <- dat %>%
    transmute(pair_barcode = Sample,
              chrom = Chromosome,
              pos = sprintf("[%s,%s]",Start_Position.bp.,End_Position.bp.),
              num_snp = Length.snp.,
              median_ratio = Median_Ratio,
              median_logr = Median_logR,
              titan_state = TITAN_state,
              titan_call = TITAN_call,
              copy_number = Copy_Number,
              major_cn = MajorCN,
              minor_cn = MinorCN,
              clonal_cluster = Clonal_Cluster,
              cellular_prevalence = Cellular_Prevalence,
              logr_copy_number = logR_Copy_Number,
              corrected_copy_number = Corrected_Copy_Number,
              corrected_call = Corrected_Call)
  
    dbWriteTable(con, Id(schema="analysis",table="titan_seg"), df, append=T)
    Sys.sleep(1)
})