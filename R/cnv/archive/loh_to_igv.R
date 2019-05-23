library(DBI)
library(tidyverse)
library(ggplot2)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

q <- "SELECT
	pair_barcode,
chrom::varchar(2),
lower(pos) as \"start\",
upper(pos)-1 as \"end\",
num_snp,
2*(median_ratio-0.75) AS median_ratio
FROM analysis.titan_seg"

qres <- dbGetQuery(con, q)
qres = qres %>% filter(complete.cases(start,end,num_snp,median_ratio))

write.table(qres, file = "loh.seg", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
