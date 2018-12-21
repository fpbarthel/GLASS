#######################################################
# Use the segmented copy number calls to derive gene-level gains/losses
# Date: 2018.11.01 
# Author: Kevin J.
#######################################################

# Directory for GLASS analysis.
mybasedir = 'Volumes/verhaak-lab/GLASS-analysis/'
datadir  = 'results/cnv/'
pattern   = '.called.seg$'

#######################################################

# Necessary packages:
library(parallel)
library(tidyverse)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(DBI)

#######################################################
# Establish connection with GLASS database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

# Downloaded the UCSC cytoband file for hg19.
cytoband_file = "/Users/johnsk/Documents/Life-History/GLASS-WG/data/ref/human_grch37_hg19_ucsc_cytoBand.txt"
cytobands = read.delim(cytoband_file, header=FALSE)

# Summarize the number of cytobands per chromosomal arm.
cytobands %>% 
  mutate(arm = substring(V4, 1, 1),
    chr_cyto = paste(V1, arm, sep='.')) %>% 
  group_by(chr_cyto) %>% 
  summarise(cyto_per_chr = n())


# Retrieve cytoband-specific copy number calls
cytoband = dbGetQuery(con,"SELECT * FROM analysis.cnv_by_cytoband")

# MERGE cytoband calls with tumor_purity and tumor_ploidy.
# DETERMINE status based on whether its copy number was greater, smaller, or equal to the sample's ploidy.




