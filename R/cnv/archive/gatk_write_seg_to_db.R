library(DBI)
library(tidyverse)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

segfiles <- list.files("results/cnv/callsegments", full.names = TRUE)
segs <- parallel::mclapply(segfiles, function(f) {
  dat <- read.delim(f, comment.char = "@", as.is= TRUE)
  dat <- dat %>%
    mutate(aliquot_barcode = substr(basename(f),1,30), pos = sprintf("[%s,%s]", START, END)) %>%
    select(aliquot_barcode, chrom = CONTIG, pos, num_points = NUM_POINTS_COPY_RATIO, log2_copy_ratio = MEAN_LOG2_COPY_RATIO, call = CALL)
  return(dat)
}, mc.cores = 8)
segs <- data.table::rbindlist(segs) %>% as.data.frame()

dbWriteTable(con, Id(schema="analysis",table="gatk_seg"), segs)
