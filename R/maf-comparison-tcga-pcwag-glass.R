#######################################################
# Comparisons of maf files from TCGA (PCAWG vs. GLASS-WG Snakemake)
# Date: 2018.07.26
# Author: Kevin J.
#######################################################

# Local directory for github repo.
# We are comparing filtered Mutect2 and Varscan2 calls (not yet maf files).
mybasedir = "~/mnt/scratchhelix/johnsk/GLASS_WG_floris/results/mutect2/vep"
# mybasedir = "~/mnt/scratchhelix/verhaak-lab/GLASS-WG/results/varscan2/final"
 
setwd(mybasedir)

#######################################################

# Load necessary packages.
library(tidyverse)
library(maftools)
library(data.table)

#######################################################
# I concatenated all the TCGA gzipped files into one. It's ~400 Mb.
test_vep_file = paste(mybasedir, "/TCGA_glioma_glass.maf.gz", sep="") 
glioma_tcga_maf = read.maf(maf = test_vep_file)

# We can also analyze smaller maf files if that's helful. Mutect2 output.
test_vep_primary_file = paste(mybasedir, "/TCGA-TQ-A8XE-TP-MF1R7Z.filtered2.anno.maf.gz", sep="")
test_vep_recurrence_file = paste(mybasedir, "/TCGA-TQ-A8XE-R1-M2U7J4.filtered2.anno.maf.gz", sep="")

# Laod the single sample maf file into memory.
TCGA_TQ_A8XE_TP_MF1R7Z = read.maf(maf = test_vep_primary_file)
TCGA_TQ_A8XE_R1_M2U7J4= read.maf(maf = test_vep_recurrence_file)

# Examine a few samples, look at the mutations present, scope out a few genes.
getSampleSummary(TCGA_TQ_A8XE_TP_MF1R7Z)

# Example of gene summary. Also, plotmafSummary().
getGeneSummary(TCGA_TQ_A8XE_TP_MF1R7Z)

# Example of the cBio lollipopPlot.
lollipopPlot(maf = TCGA_TQ_A8XE_TP_MF1R7Z, gene = 'IDH1', AACol = 'Protein_Change', showMutationRate = TRUE)

# Example of some important genes to glioma.
genes = c("IDH1", "TP53", "TET1", "EGFR", "PDGFRA")

#We will draw oncoplots for top ten mutated genes.
oncoplot(maf = laml, top = 10, fontSize = 12)

# It may be possible to load all of the TCGA LGG and GBM data in using this package.
