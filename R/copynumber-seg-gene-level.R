#######################################################
# Use the segmented copy number calls to derive gene-level gains/losses
# Date: 2018.11.01 
# Author: Kevin J.
#######################################################

# Directory for GLASS analysis.
mybasedir = 'Volumes/verhaak-lab/GLASS-analysis/'
datadir  = 'results/cnv/callsegments'
pattern   = '.called.seg$'

#######################################################

# Necessary packages:
library(parallel)
library(tidyverse)
library(GenomicRanges)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(DBI)

#######################################################
# Establish connection with Floris' database.
# con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

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
    mutate(sample_id = gsub(".called.seg", "", basename(f)))
  
  return(dat)
  
}, mc.cores=20)

## Combine all the samples from the GLASS-WG cohort.
glass_cn = data.table::rbindlist(cn_dat)

## Create a GRanges friendly chromosome identifier.
glass_cn = glass_cn %>% 
  mutate(chr = CONTIG) %>% 
  dplyr::select(-CONTIG)

## Create a GRanges object from this data to merge with the Cytoband data.
glass_cn_gr = makeGRangesFromDataFrame(glass_cn,
                                         keep.extra.columns = T)

## Pull gene metadata for a supervised list of genes:
genes = c('MDM4', 'AKT3', 'PDGFRA', 'PTEN', 'EGFR', 'MET', 'NF1', 'VEGF', 'IDH1', 'IDH2', 'PIK3CA', 'PIK3R1', 'CDKN2A', 'CDKN2C', 'CDK4', 
          'CDK6', 'RB1', 'MGMT', 'TERT', 'MYCNP', 'GLI2', 'FGFR3/TACC3', 'MYB', 'KIAA154/BRAF', 'MYBL1', 'MYC', 'PTCH1', 'CND1', 'CCND2', 
          'MDM2', 'PARK2', 'FGFR2', 'IRS2', 'PTPRD', 'MLH1' , 'MSH2', 'MSH6', 'PMS2', 'ERBB2', 'ARF', 'MDM2', 'TP53', 'CDKN2B', 'BRAF', 
          'HIF1', 'YKL40', 'ELDT1', 'ATM', 'ATRCTLA4', 'PD1', 'H3F3A', 'DAXX', 'PARP', 'PTEN', 'STAT3')

## Define gene metadata.
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
tx = transcripts(txdb)

### Gene to TX table.
gid2gene = toTable(org.Hs.egSYMBOL)
transcript_ids = dplyr::select(as.data.frame(tx), TXID=tx_id, TXNAME=tx_name, chr=seqnames, start, end, strand) %>% 
  mutate(chr = as.character(gsub("chr", "", chr))) %>% 
  filter(chr %in% c(seq(1,22), "X"))

## Join transcripts and genes.
tx2gene = AnnotationDbi::select(txdb, tx$tx_name, "GENEID", "TXNAME") %>% 
  left_join(transcript_ids) %>% 
  left_join(dplyr::select(gid2gene, GENEID=gene_id, SYMBOL=symbol))

## Reduce to one row per gene.
gene_meta = tx2gene %>% mutate(gene = SYMBOL) %>% 
  filter(complete.cases(gene,chr), gene %in% genes) %>% 
  mutate(chr = as.numeric(chr)) %>%
  group_by(gene, chr) %>% 
  summarize(start = min(start), end = max(end), strand = paste(unique(strand), collapse = "/")) %>% 
  ungroup() %>% 
  filter(complete.cases(gene, chr)) 
  
# Generate a GRange object for glioma-related driver genes.
glioma_genes_gr = makeGRangesFromDataFrame(gene_meta,
                                       keep.extra.columns = T)

# Subset the copy number calls to only those genes of interest.
# If regions in dataset1 overlap with other regions in dataset1, they will end up being "collapsed" into a single region in the output of intersect().
# Intersect doesn't work here because it returns only the overlapping regions.

## Try to do this the data.table way
copynum_calls = data.table(
  chr = as.numeric(seqnames(glass_cn_gr)),
  start = start(glass_cn_gr), 
  end =  end(glass_cn_gr),
  sample_id = glass_cn_gr$sample_id,
  copy_ratio = glass_cn_gr$MEAN_LOG2_COPY_RATIO,
  CALL = glass_cn_gr$CALL)
setkey(copynum_calls, chr, start, end)

glioma_gene_coord = data.table(
  chr = as.numeric(seqnames(glioma_genes_gr)),
  start = start(glioma_genes_gr), 
  end =  end(glioma_genes_gr),
  gene = glioma_genes_gr$gene)
setkey(glioma_gene_coord, chr, start, end)

# Find overlaps
cn_gene_overlap <- foverlaps(copynum_calls, glioma_gene_coord, type = "any", nomatch = 0)
cn_gene_overlap_df <- as.data.frame(cn_gene_overlap)

## Determine the copy number status per gene per sample.
# Perofrm the analysis Floris mentioned with gatk_calls.
glass_gene_status = cn_gene_overlap_df %>%
  mutate(calls = recode(CALL, "-"="-1","+"="1")) %>% 
  mutate(widths = abs(end-start)) %>% 
  group_by(gene, sample_id) %>% 
  summarise(weighted_calls = sum(widths*as.numeric(calls))/sum(widths)) %>% 
  mutate(call = ifelse(weighted_calls<=-.3, "deletion", ifelse(weighted_calls>=.3, "amplification", "neutral"))) %>% 
  ungroup() %>% 
  dplyr::select(-weighted_calls) %>% 
  spread(gene, call)

# Create a table gene x sample.
rownames(glass_gene_status) <- NULL
eventcnv = glass_gene_status %>%
  gather(gene, gatk_calls, -sample_id) %>%
  spread(sample_id, gatk_calls)

## Determine the copy number ratio per gene per sample.
glass_gene_ratio = cn_gene_overlap_df %>%
  mutate(widths = abs(end-start)) %>% 
  group_by(gene, sample_id) %>% 
  summarise(weighted_log2 = sum(widths*copy_ratio)/sum(widths)) %>% 
  mutate(call = ifelse(weighted_log2<=(-0.3), "deletion", ifelse(weighted_log2>=0.3, "amplification", "neutral"))) %>% 
  ungroup() %>% 
  dplyr::select(-call) %>% 
  spread(gene, weighted_log2)

# Gene x sample table.
rownames(glass_gene_ratio) <- NULL
eventcnv2 = glass_gene_ratio %>%
  gather(gene, log2_ratio, -sample_id) %>%
  spread(sample_id, log2_ratio)


# Write out file for 
# write.table(glass_gene_ratio, file = "/Users/johnsk/Documents/glass_copy_number_status.txt", sep="\t", row.names = F, col.names = T, quote = F)

# Load the meta data derived from the aliquot_barcode.
meta = cnsamples %>% 
  mutate(subject_id = substr(sample_id, 1, 12),
         sample_type = substr(sample_id, 14, 15),
         plot_id = substr(sample_id, 1, 15),
         cohort_id = substr(sample_id, 6, 7))

# Floris provided a driver gene table.
drivergenes = openxlsx::read.xlsx('/Users/johnsk/Documents/glioma_driver_genes OLD.xlsx')

# In a diploid genome, a single-copy gain in a perfectly pure, homogeneous sample has a copy ratio of 3/2. 
# In log2 scale, this is log2(3/2) = 0.585, and a single-copy loss is log2(1/2) = -1.0.
genecnv_gatk = eventcnv %>% 
  as.data.frame() %>% 
  gather('sample_id', 'gatk_calls', -gene) %>%
  left_join(drivergenes, by="gene") %>%
  filter(complete.cases(pathway)) %>%
  mutate(CNV = factor(ifelse(effect == 'amplification' & gatk_calls == 'amplification', "+1", ifelse(effect == 'deletion' & gatk_calls == 'deletion', "-1", "0")), levels = c("-1", "0", "+1"))) %>%  
  left_join(meta, by="sample_id") %>% 
  filter(!grepl("-NB-", sample_id))

genecnv_ratio = eventcnv2 %>% 
  as.data.frame() %>% 
  gather('sample_id', 'log2_ratio', -gene) %>%
  left_join(drivergenes, by="gene") %>%
  filter(complete.cases(pathway)) %>%
  mutate(CNV = factor(ifelse(effect == 'amplification' & log2_ratio > 0.3, "+1", ifelse(effect == 'deletion' & log2_ratio < -.3, "-1", "0")), levels = c("-1", "0", "+1"))) %>%  
  left_join(meta, by="sample_id") %>% 
  filter(!grepl("-NB-", sample_id)) %>% 
  filter(!grepl("-WXS-", sample_id))

##### Plot Gene CNV - gatk status

# genecnv = genecnv %>% filter(gene != "RB1")
genecnv_gatk$pathway = factor(genecnv_gatk$pathway, levels = c("Apoptosis", "Cell cycle", "PI3K-RTK-MAPK"), labels = c("Apoptosis", "Cell cycle", "PI3K-\nRTK-MAPK"))
genecnv_ratio$pathway = factor(genecnv_ratio$pathway, levels = c("Apoptosis", "Cell cycle", "PI3K-RTK-MAPK"), labels = c("Apoptosis", "Cell cycle", "PI3K-\nRTK-MAPK"))

p_gatk = ggplot() + 
  geom_tile(data = genecnv_gatk, aes(x = plot_id, y = gene, fill = CNV)) +
  scale_fill_manual(values = c("-1" = "blue", "0" = NA, "+1" = "red")) +
  labs(y="Gene", fill = "Copy number") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme_bw()

p_gatk


##### Plot Gene CNV - log2 ratio
p_ratio = ggplot() + 
  geom_tile(data = genecnv_ratio, aes(x = plot_id, y = gene, fill = CNV)) +
  scale_fill_manual(values = c("-1" = "blue", "0" = NA, "+1" = "red")) +
  labs(y="Gene", fill = "Copy number") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

p_ratio



##### For maf tools: tab separated file with: gene name, Sample name and copy number status (either 'Amp' or 'Del').

maf_tools_input = genecnv_ratio %>% 
  dplyr::select(gene, sample_id, CNV) %>% 
  mutate(CN = recode(CNV, "0"="0","+1"="Amp", "-1"="Del")) %>% 
  dplyr::select(-CNV) %>% 
  filter(grepl("-WGS-", sample_id)) %>%
  filter(!grepl("-NB-", sample_id)) %>% 
  filter(grepl())
           
write.table(maf_tools_input, file = "/Users/johnsk/Documents/maf_tools_input_wgs.txt", sep="\t", row.names = F, col.names = T, quote = F)


# Breakdown by IDH status:
glass_wg_idh_mut = openxlsx::read.xlsx('/Users/johnsk/Documents/Life-History/tmp/glass_wg_subtypes.xlsx')

cnv_idh_status = maf_tools_input %>% 
  left_join(glass_wg_idh_mut, by=c("sample_id"="aliquot_id")) %>% 
#  spread(gene, CN) %>% 
  mutate(sample_type = substr(sample_id, 14, 15),
         plot_id = substr(sample_id, 1, 15)) 

p_ratio_subtype = ggplot() + 
  geom_tile(data = cnv_idh_status, aes(x = plot_id, y = gene, fill = CN)) +
  scale_fill_manual(values = c("Del" = "blue", "0" = NA, "Amp" = "red")) +
  labs(y="Gene", fill = "Copy number") +
  facet_grid(~Sequencing.IDH) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

p_ratio_subtype

