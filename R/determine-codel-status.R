#######################################################
# 1p/19q calling in the GLASS-WXS & GLASS-WG datasets
# Date: 2018.10.30 
# Author: Kevin J.
#######################################################

# Downloaded the UCSC cytoband file.
cytoband_file = "/Users/johnsk/Documents/Life-History/GLASS-WG/data/ref/human_grch37_hg19_ucsc_cytoBand.txt"

# Directory for GLASS analysis.
mybasedir = 'Volumes/verhaak-lab/GLASS-analysis/'
datadir   = 'results/cnv/callsegments'
pattern   = '.called.seg$'

#######################################################

# Necessary packages:
library(parallel)
library(tidyverse)
library(GenomicRanges)
library(data.table)

setwd(mybasedir)

#######################################################
## Read in an example "*.called.seg" file to test the calling.
files = list.files(datadir, full.names = T, pattern = pattern, recursive=T)
# If it is desirable to include the sample names.
cnsamples = data.frame(sample_id=gsub(".called.seg", "", basename(files)), library_type = substring(basename(files), 21, 23))

# The first 88 rows of each file represent a header.
cn_dat = mclapply(files, function(f){
  dat = tryCatch(read.delim(f,as.is=T, header=T, skip = 88), error=function(e) e)
  if(inherits(dat,'error')) {
    message(f, '\n', dat, '\n')
    return()
  }
# Truncate the file name to just the sample_id.
  dat = dat %>%
    mutate(sample_id = gsub(".called.seg", "", basename(f))) %>% 
    filter(CONTIG%in%c("1", "19"))

  return(dat)

  }, mc.cores=20)

## Combine all the samples from the GLASS-WG cohort.
glass_cn = data.table::rbindlist(cn_dat)

## Create a GRanges friendly chromosome identifier.
glass_cn = glass_cn %>% 
  mutate(chromosome = paste0("chr", CONTIG)) %>% 
  select(-CONTIG)

## Inspect the averages, min, and max for how the algorithm called copy number changes. Interestingly, it appears that there must be some sample-specific normalization
# as the there are negative values for `0` as well as very high `+` values. These may be classified somehow as artifacts.
glass_cn %>% 
  group_by(CALL) %>% 
  summarise(min_call = min(MEAN_LOG2_COPY_RATIO),
            max_call = max(MEAN_LOG2_COPY_RATIO),
            median_call = median(MEAN_LOG2_COPY_RATIO))

## Create a GRanges object from this data to merge with the Cytoband data.
glass_cn_gr = makeGRangesFromDataFrame(glass_cn,
                                         keep.extra.columns = T)

## hg19 cytobands looking for chromosomes 1 and 19.
cytobands = read.delim(cytoband_file, header=FALSE)

## Create a genomic range object for chr1 cytoband p.
hg19_1p = cytobands %>% 
  filter(V1=='chr1') %>% 
  filter(grepl("p", V4))
colnames(hg19_1p) <- c("chromosome", "start", "end", "cytoband", "stain")
hg19_1p_GR = makeGRangesFromDataFrame(hg19_1p)

# Create a genomic range object for chr19 cytoband q.
hg19_19q = cytobands %>% 
  filter(V1=='chr19') %>% 
  filter(grepl("q", V4))
colnames(hg19_19q) <- c("chromosome", "start", "end", "cytoband", "stain")
hg19_19q_GR = makeGRangesFromDataFrame(hg19_19q)

## Need to only use those coordinates corresponding to the chromosome arm.
## We don't want to consider segments that run on into a different arm.
# Restrict ranges to the p arm of chromosome 1:
chr1p_calls = subsetByOverlaps(glass_cn_gr, hg19_1p_GR)
chr1p_calls_restrict <- restrict(chr1p_calls, end=max(end(hg19_1p_GR)), keep.all.ranges=TRUE)
chr1p_calls_restrict <- chr1p_calls_restrict[width(chr1p_calls_restrict) != 0 | chr1p_calls_restrict == chr1p_calls]

# Restrict ranges to the q arm of chromosome 19:
chr19q_calls = subsetByOverlaps(glass_cn_gr, hg19_19q_GR)
chr19q_calls_restrict <- restrict(chr19q_calls, start=min(start(hg19_19q_GR)), keep.all.ranges=TRUE)
chr19q_calls_restrict <- chr19q_calls_restrict[width(chr19q_calls_restrict) != 0 | chr19q_calls_restrict == chr19q_calls]

## Combine the ranges:
chr1p_19q_calls_restrict <- c(chr1p_calls_restrict, chr19q_calls_restrict)
         
# Use cut-off thresholds of log2(1.1) and log2(0.9) for calling deletion or amplification.
chr1p_19q_calls_subject <- data.frame(seqnames=seqnames(chr1p_19q_calls_restrict),
                 starts=start(chr1p_19q_calls_restrict)-1,
                 ends=end(chr1p_19q_calls_restrict),
                 widths = width(chr1p_19q_calls_restrict),
                 sample_id = elementMetadata(chr1p_19q_calls_restrict)$sample_id,
                 num_points_copy_ratio = elementMetadata(chr1p_19q_calls_restrict)$NUM_POINTS_COPY_RATIO,
                 mean_log2_copy_ratio = elementMetadata(chr1p_19q_calls_restrict)$MEAN_LOG2_COPY_RATIO,
                 cn_call = elementMetadata(chr1p_19q_calls_restrict)$CALL)

## Determine the status of 1p/19q codeletion.
glass_1p19q_status = chr1p_19q_calls_subject %>%  
  group_by(seqnames, sample_id) %>% 
  summarise(weighted_log2 = sum(widths*mean_log2_copy_ratio)/sum(widths)) %>% 
  mutate(chr_call = ifelse(weighted_log2<=log2(0.9), "deletion", "intact")) %>% 
  select(-weighted_log2) %>% 
  spread(seqnames, chr_call) %>% 
  mutate(codel_status = ifelse(chr1=="deletion" & chr19=="deletion", "codel", "noncodel"))

# Briefly, enumerate the codel_status:  
table(glass_1p19q_status$codel_status)

# Write file out for deletion only calls. We did not consider the amplification, only whether it was deleted.
write.table(glass_1p19q_status, file = "/Users/johnsk/Documents/glass_1p19q_codeletion_status.txt", sep="\t", row.names = F, col.names = T, quote = F)


# Output just TCGA WGS codel_status:
TCGA_1p_19q = glass_1p19q_status %>% 
#  filter(grepl("TCGA", sample_id)) %>% 
  filter(grepl("WGS", sample_id)) %>% 
  filter(!grepl("-NB-", sample_id))

write.table(TCGA_1p_19q, file = "/Users/johnsk/Documents/tcga_wgs_1p19q_codeletion_status.txt", sep="\t", row.names = F, col.names = T, quote = F)
