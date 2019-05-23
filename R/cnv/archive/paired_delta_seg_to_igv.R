library(DBI)
library(tidyverse)
library(ggplot2)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

q<- "SELECT
	tumor_pair_barcode,
chrom::varchar(2),
lower(pos) as \"start\",
upper(pos)-1 as \"end\",
0 as num_snps,
delta_cn
FROM analysis.titan_seg_paired_delta pa"

qres <- dbGetQuery(con, q)
seg = qres %>% filter(complete.cases(start,end,num_snps))

#write.table(seg, file = "diff.seg", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

q <- "SELECT tumor_pair_barcode, diamond_set::integer
FROM analysis.titan_seg_paired_comparison ts
LEFT JOIN biospecimen.aliquots al ON al.aliquot_barcode = ts.tumor_barcode_a
LEFT JOIN clinical.surgeries cl ON cl.sample_barcode = al.sample_barcode"

qres <- dbGetQuery(con, q)

seg <- seg %>%
  left_join(qres) %>% 
  filter(diamond_set==1) %>%
  select(-diamond_set)

write.table(seg, file = "results/cnv/gistic/input.seg", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

markers = data.frame(id = 1:(2*nrow(seg)), chr = c(seg$chrom, seg$chrom), pos = c(seg$start, seg$end), stringsAsFactors = FALSE) %>% distinct()
write.table(markers, file = "results/cnv/gistic/markers.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#write.table(qres, file = "diff.annoseg.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
