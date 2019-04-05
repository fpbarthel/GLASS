#######################################################
# Segments to database
#######################################################

options(scipen=999)

## Parse snakemake
if(exists("snakemake")) {
  files = snakemake@input[["seg"]]
  outfn = snakemake@output[["tsv"]]
} else {
  files = list.files("results/cnv/callsegments/", recursive = T, pattern = "called.seg", full.names = T)[1:10]
}

# Necessary packages:
library(parallel)
library(tidyverse)
library(data.table)
library(DBI)

segs = lapply(files, function(f){
  dat <- read.delim(f, comment.char = "@", as.is= TRUE)
  dat <- dat %>%
    mutate(aliquot_barcode = substr(basename(f),1,30), pos = sprintf("[%s,%s]", START, END)) %>%
    select(aliquot_barcode, chrom = CONTIG, pos, num_points = NUM_POINTS_COPY_RATIO, log2_copy_ratio = MEAN_LOG2_COPY_RATIO, call = CALL)
  return(dat)
})
segs <- data.table::rbindlist(segs) %>% as.data.frame()

write.table(segs, file = outfn, quote = F, sep = "\t", row.names = FALSE, col.names = FALSE)