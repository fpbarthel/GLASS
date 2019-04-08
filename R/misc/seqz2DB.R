### push titan seg into db

library(tidyverse)
library(DBI)
library(odbc)
con <- DBI::dbConnect(odbc::odbc(), "GLASSv2")

seg <- read.delim("results/sequenza/glass_seqz_segments.tsv", as.is = TRUE)
pp <- read.delim("results/sequenza/glass_seqz_purity_ploidy.tsv", as.is = TRUE)

pp <- pp %>% select(pair_barcode, cellularity, ploidy, slpp = SLPP)
seg <- seg %>% transmute(pair_barcode,
                         chrom = ifelse(chromosome=='X',23,as.integer(chromosome)),
                         pos = sprintf("[%s,%s]",start.pos,end.pos),
                         baf = Bf,
                         baf_n = N.BAF,
                         baf_sd = sd.BAF,
                         ratio = depth.ratio,
                         ratio_n = N.ratio,
                         ratio_sd = sd.ratio,
                         copy_number = CNt,
                         major_cn = A,
                         minor_cn = B,
                         log_posterior_proba = LPP)

dbWriteTable(con, Id(schema="variants",table="seqz_seg"), seg, append=T)
dbWriteTable(con, Id(schema="variants",table="seqz_params"), pp, append=T)

# segfiles <- list.files('results/cnv/titanfinal/seg', full.names = TRUE)
# 
# lapply(segfiles, function(f){
#   message(f)
#   dat <- read.delim(f, as.is=T, header=T, row.names = NULL)
#   df <- dat %>%
#     transmute(pair_barcode = Sample,
#               chrom = Chromosome,
#               pos = sprintf("[%s,%s]",Start_Position.bp.,End_Position.bp.),
#               num_snp = Length.snp.,
#               median_ratio = Median_Ratio,
#               median_logr = Median_logR,
#               titan_state = TITAN_state,
#               titan_call = TITAN_call,
#               copy_number = Copy_Number,
#               major_cn = MajorCN,
#               minor_cn = MinorCN,
#               clonal_cluster = Clonal_Cluster,
#               cellular_prevalence = Cellular_Prevalence,
#               logr_copy_number = logR_Copy_Number,
#               corrected_copy_number = Corrected_Copy_Number,
#               corrected_call = Corrected_Call)
#   
#     dbWriteTable(con, Id(schema="analysis",table="titan_seg"), df, append=T)
#     Sys.sleep(1)
# })