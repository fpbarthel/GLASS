##################################################
# Prepare input files for running GISTIC seperately for primaries and recurrences
# Ignores multi-sector samples (one sample per patient and timepoint)
# Updated: 2019.04.19
# Author: Floris B
##################################################

library(DBI)
library(tidyverse)
library(ggplot2)

con <- DBI::dbConnect(odbc::odbc(), "GLASSv2")

q <- read_file("sql/cnv/gistic_prepare.sql")

qres <- dbGetQuery(con, q)
seg = qres %>% filter(complete.cases(start,end,num_snps))

seg_p = seg %>% filter(sample_type == "P") %>% select(-sample_type)
seg_r = seg %>% filter(sample_type == "R") %>% select(-sample_type)

write.table(seg_p, file = "results/gistic2/primary.seg", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(seg_r, file = "results/gistic2/recurrence.seg", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

markers = data.frame(id = 1:(2*nrow(seg)), chr = c(seg$chrom, seg$chrom), pos = c(seg$start, seg$end), stringsAsFactors = FALSE) %>% distinct()
write.table(markers, file = "results/gistic2/markers.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)