#######################################################
# 1p/19q calling in the GLASS-WXS & GLASS-WG datasets.
# Date: 2018.11.15 
# Author: Kevin J.
#######################################################

# Downloaded the UCSC cytoband file for hg19.
cytoband_file = "/Users/johnsk/Documents/Life-History/GLASS-WG/data/ref/human_grch37_hg19_ucsc_cytoBand.txt"

# Directory for GLASS analysis and copy number files.
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
## Read in the "*.called.seg" files.
files = list.files(datadir, full.names = T, pattern = pattern, recursive=T)

## The first 88 rows of each ".called.seg" file represent a header.
cn_dat = mclapply(files, function(f){
  dat = tryCatch(read.delim(f,as.is=T, header=T, skip = 88), error=function(e) e)
  if(inherits(dat,'error')) {
    message(f, '\n', dat, '\n')
    return()
  }
# Truncate the file name to just the sample_id. We are only interested in 1p and 19q.
  dat = dat %>%
    mutate(sample_id = gsub(".called.seg", "", basename(f))) %>% 
    filter(CONTIG%in%c("1", "19"))

  return(dat)

  }, mc.cores=20)

## Combine all the samples from the GLASS-WG cohort.
glass_cn = data.table::rbindlist(cn_dat)

## Make sure GRanges chromosome identifiers line up with one another.
glass_cn = glass_cn %>% 
  # swap column labels.
  mutate(chr = CONTIG) %>% 
  select(-CONTIG)

## Inspect the averages, min, and max for how the algorithm called copy number changes. It appears that there must be some sample-specific normalization
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
hg19_1p_GR = makeGRangesFromDataFrame(hg19_1p, keep.extra.columns = T)

# Create a genomic range object for chr19 cytoband q.
hg19_19q = cytobands %>% 
  filter(V1=='chr19') %>% 
  filter(grepl("q", V4))
colnames(hg19_19q) <- c("chromosome", "start", "end", "cytoband", "stain")
hg19_19q_GR = makeGRangesFromDataFrame(hg19_19q, keep.extra.columns = T)

#######################
# OPTION 1 codel calling. 
# quick and dirty.
# Average over all cytobands results in false positives.
#######################
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
  # log2(0.9) represents single copy loss. -.3 is a more conservative threshold.
  mutate(chr_call = ifelse(weighted_log2<=log2(0.9), "deletion", "intact")) %>% 
  select(-weighted_log2) %>% 
  spread(seqnames, chr_call) %>% 
  mutate(codel_status = ifelse(chr1=="deletion" & chr19=="deletion", "codel", "noncodel"))

# Briefly, enumerate the codel_status:  
table(glass_1p19q_status$codel_status)


#######################
## OPTION 2 (MORE ACCURATE)
## Accounts for issues with partial 1p 19q codeletion
#######################
# For each cytoband it may be advantageous to look to see whether there are local or total deletions.
hg19_1p19q = c(hg19_1p_GR, hg19_19q_GR)

## Try to do this the data.table way.
copynum_calls = data.table(
  chr = as.numeric(as.character(seqnames(glass_cn_gr))),
  start = start(glass_cn_gr), 
  end =  end(glass_cn_gr),
  sample_id = glass_cn_gr$sample_id,
  copy_ratio = glass_cn_gr$MEAN_LOG2_COPY_RATIO,
  CALL = glass_cn_gr$CALL)
setkey(copynum_calls, chr, start, end)

# Make sure chromosomes are correct ("1" and "19").
glioma_1p19q_coord = data.table(
  chr = as.numeric(seqnames(hg19_1p19q)),
  start = start(hg19_1p19q), 
  end =  end(hg19_1p19q),
  cytoband = hg19_1p19q$cytoband)
setkey(glioma_1p19q_coord, chr, start, end)

# Find overlaps for hg19 1p19q in the GLASS samples.
cn_1p19q_overlap <- foverlaps(copynum_calls, glioma_1p19q_coord, type = "any", nomatch = 0)
cn_1p19q_overlap_df <- as.data.frame(cn_1p19q_overlap)
cn_1p19q_overlap_df$full_cytoband = paste0(cn_1p19q_overlap_df$chr, cn_1p19q_overlap_df$cytoband)
cn_1p19q_overlap_df$cytoband = NULL

## Determine the copy number ratio per gene per sample. 
glass_cytoband_ratio = cn_1p19q_overlap_df %>%
  mutate(widths = abs(end-start)) %>% 
  group_by(full_cytoband, sample_id) %>% 
  summarise(weighted_log2 = sum(widths*copy_ratio)/sum(widths)) %>% 
  # We are using more stringent log2 ratio  values for copy number to insure true deletion.
  mutate(call = ifelse(weighted_log2<=(-0.3), "deletion", ifelse(weighted_log2>=0.3, "amplification", "neutral"))) %>% 
  ungroup() %>% 
  dplyr::select(-call) %>% 
  spread(full_cytoband, weighted_log2)

## There are 46 segments across 1p 19q. Focus on enumerating the deleted segments.
glass_cytoband_call = cn_1p19q_overlap_df %>%
  mutate(widths = abs(end-start)) %>% 
  group_by(full_cytoband, sample_id) %>% 
  summarise(weighted_log2 = sum(widths*copy_ratio)/sum(widths)) %>% 
  # We are using more stringent log2 ratio  values for copy number to insure true deletion.
  mutate(call = ifelse(weighted_log2<=(-0.3), "deletion", ifelse(weighted_log2>=0.3, "amplification", "neutral"))) %>% 
  ungroup() %>% 
  dplyr::select(-weighted_log2) %>% 
  group_by(sample_id) %>% 
  filter(!grepl("-NB-", sample_id)) %>% 
  filter(call == "deletion") %>% 
  summarise(del_counts = n())

# Check inconsistencies for 1p19q:
glass_codel = glass_cytoband_call %>% 
  filter(del_counts > 38) %>% 
  mutate(subject_id = substr(sample_id, 1, 12)) %>% 
  mutate(status_1p19q = "codel") %>% 
  select(-sample_id, -del_counts) %>% 
  # Poor performing UCSF samples.
  filter(!subject_id%in%c("GLSS-SF-0081", "GLSS-SF-0006")) %>% 
  distinct()

# Write out file with codel tumors.
write.table(glass_codel, file ="/Users/johnsk/Documents/Life-History/glass-subject-level-codel-status.txt", quote=FALSE, sep='\t', row.names = F)
