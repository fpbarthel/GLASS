#######################################################
# Enumerate cumulative coverage per aliquot for WGS/WXS
# Date: 2018.11.06 
# Author: Kevin J., FP Barthel
#######################################################

options(scipen=999)

## Parse snakemake
if(exists("snakemake")) {
  files = snakemake@input[["metrics"]]
  outfn = snakemake@output[["tsv"]]
} else {
  files = list.files("results/align/wgsmetrics", recursive = T, pattern = "WgsMetrics.txt", full.names = T) # list("results/align/wgsmetrics/GLSS-DK-0012-NB-01D-WXS-ABCB18.WgsMetrics.txt", "results/align/wgsmetrics/GLSS-DK-0003-TP-01D-WXS-E43D26.WgsMetrics.txt")
}

# Necessary packages:
library(parallel)
library(tidyverse)
library(data.table)
library(DBI)

# The first 10 rows of each file represent a header of additional information.
cov_dat = lapply(files, function(f){
  dat = tryCatch(read.delim(f,as.is=T, header=T, row.names = NULL, skip = 10), error=function(e) e)
  if(inherits(dat,'error')) {
    message(f, '\n', dat, '\n')
    return()
  }
  # Truncate the file name to just the sample_id.
  dat = dat %>%
    mutate(sample_id = gsub(".WgsMetrics.txt", "", basename(f)),
           high_quality_coverage_count = as.numeric(high_quality_coverage_count)) # %>%  
  #    filter(coverage!="0") # Filter out those bases with `0` coverage.
  
  return(dat)
})#, mc.cores=20)

## Combine all the samples from the GLASS cohort.
glass_cov = data.table::rbindlist(cov_dat)

# Cumulatively add the number of bases at each level:
glass_samples_cumulative_cov = glass_cov %>% 
  group_by(sample_id) %>% 
  mutate(cumulative_coverage = rev(cumsum(rev(high_quality_coverage_count)))) %>% 
  select(aliquot_barcode = sample_id, coverage, high_quality_coverage_count, cumulative_coverage)

# Write output as one table or a table for each file:
write.table(glass_samples_cumulative_cov, file = outfn, quote = F, sep = "\t", row.names = FALSE, col.names = FALSE)