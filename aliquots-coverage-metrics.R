#######################################################
# Enumerate cumulative coverage per aliquot for WGS/WXS
# Date: 2018.11.06 
# Author: Kevin J.
#######################################################

# Directory for GLASS analysis.
mybasedir = '/Volumes/fastscratch/verhaak-lab/GLASS-WG'
datadir  = 'results/align/wgsmetrics/'
pattern   = '.WgsMetrics.txt$'

#######################################################

# Necessary packages:
library(parallel)
library(tidyverse)
library(data.table)

#######################################################

## Read in an example "*.WgsMetrics.txt" file to test the calling.
files = list.files(datadir, full.names = T, pattern = pattern, recursive=T)

# If it is desirable to include the sample names.
samples = data.frame(sample_id=gsub(".WgsMetrics.txt", "", basename(files)), library_type = substring(basename(files), 21, 23))

# The first 10 rows of each file represent a header of additional information.
cov_dat = mclapply(files, function(f){
  dat = tryCatch(read.delim(f,as.is=T, header=T, row.names = NULL, skip = 10), error=function(e) e)
  if(inherits(dat,'error')) {
    message(f, '\n', dat, '\n')
    return()
  }
  # Truncate the file name to just the sample_id.
  dat = dat %>%
    mutate(sample_id = gsub(".WgsMetrics.txt", "", basename(f))) # %>%  
#    filter(coverage!="0") # Filter out those bases with `0` coverage.
  
  return(dat)
  
}, mc.cores=20)

## Combine all the samples from the GLASS cohort.
glass_cov = data.table::rbindlist(cov_dat)

# Cumulatively add the number of bases at each level:
glass_samples_cumulative_cov = glass_cov %>% 
  group_by(sample_id) %>% 
  mutate(cumulative_coverage = rev(cumsum(rev(high_quality_coverage_count)))) %>% 
  select(aliquot_id = sample_id, coverage, high_quality_coverage_count, cumulative_coverage)

# Write output as one table or a table for each file:
write.table(glass_samples_cumulative_cov, file = "/Users/johnsk/Documents/glass-cumulative-coverage.txt", sep="\t", row.names = F, col.names = T, quote = F)


cov_dat = mclapply(files, function(f){
  dat = tryCatch(read.delim(f,as.is=T, header=T, row.names = NULL, skip = 10), error=function(e) e)
  if(inherits(dat,'error')) {
    message(f, '\n', dat, '\n')
    return()
  }
  # Truncate the file name to just the sample_id.
  dat = dat %>%
    mutate(sample_id = gsub(".WgsMetrics.txt", "", basename(f))) # %>%  
  #    filter(coverage!="0") # Filter out those bases with `0` coverage.
  
  return(dat)
  
}, mc.cores=20)
