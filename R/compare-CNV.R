#######################################################
# Compare CNV
# Date: 2018.08.16
# Author: Floris B
#######################################################

library(tidyverse)

setwd("/projects/barthf/GLASS-WG/")

tcgaseg = list("../PCTL/data/seg/stddata__2016_01_28/gdac.broadinstitute.org_GBM.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0/GBM.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt",
           "../PCTL/data/seg/stddata__2016_01_28/gdac.broadinstitute.org_LGG.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0/LGG.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt") %>%
  map_df(read.delim, as.is = T)  %>%
  mutate(legacy_sample_id = substr(Sample,1,16))

lowpseg = read.delim("sandbox/gdac.broadinstitute.org_GBMLGG.Merge_cna__illuminahiseq_dnaseqc__hms_harvard_edu__Level_3__segmentation__seg.Level_3.2016012800.0.0/GBMLGG.cna__illuminahiseq_dnaseqc__hms_harvard_edu__Level_3__segmentation__seg.seg.txt", as.is=T) %>%
  mutate(legacy_sample_id = substr(Sample,1,16))

samplemap = read.delim("data/manifest/tcga/samples.tsv", as.is = T)

tcgaseg = tcgaseg %>% 
  left_join(samplemap) %>%
  filter(complete.cases(sample_id)) %>%
  select(Sample = sample_id, Chromosome, Start, End, Num_Probes, Segment_Mean)

lowpseg = lowpseg %>% 
  left_join(samplemap) %>%
  filter(complete.cases(sample_id)) %>%
  select(Sample = sample_id, Chromosome, Start, End, Num_Probes, Segment_Mean)

gatkfiles = list.files("results/cnv/callsegments", pattern = "*seg", full.names = T)
gatkseg = gatkfiles %>%
  map(read.delim, as.is = T, comment.char = "@")

names(gatkseg) = gatkfiles
gatkseg = gatkseg %>% 
  map2_df(gatkfiles, ~mutate(.x, Sample = gsub("results/cnv/callsegments/(.*)-\\w{6}.called.seg", "\\1", .y))) %>%
  select(Sample, Chromosome = CONTIG, Start = START, End = END, Num_Probes = NUM_POINTS_COPY_RATIO, Segment_Mean = MEAN_LOG2_COPY_RATIO) %>%
  filter(Sample %in% tcgaseg$Sample)

tcgaseg = tcgaseg %>%
  filter(Sample %in% gatkseg$Sample)

lowpseg = lowpseg %>%
  filter(Sample %in% gatkseg$Sample)


write.table(tcgaseg, file="sandbox/tcga_tmp.seg", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(gatkseg, file="sandbox/gatk_tmp.seg", row.names = F, col.names = T, quote = F, sep = "\t")

# library(copynumber)
# 
# plotSampleByID <- function(seg, id) {
#   
# }
# 
# write.table(tcgaseg, file="sandbox/tcga_tmp.seg", row.names = F, col.names = T, quote = F, sep = "\t")
# plotCBSsegments("sandbox/tcga_tmp.seg")
# 
# 
# plotSample(segments = tcgaseg[tcgaseg$Sample=="TCGA-06-0125-NB", c(2:6)])
