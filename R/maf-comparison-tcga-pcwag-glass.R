#######################################################
# Comparisons of maf files from TCGA (PCAWG vs. GLASS-WG Snakemake)
# Date: 2018.08.15
# Author: Kevin J.
#######################################################
# Local directory for github repo.
# Use Kadir's linker file to identify the samples that overlap between PCAWG and GLASS-WG.
mybasedir = "/Users/johnsk/Documents/Life-History/"
setwd(mybasedir)

# Kadir's linker files.
pcawg_id_file <- "PCAWG May 2016 Data Release.xlsx"
tcga_id_file <- "pcawg_specimen_histology_August2016_v8.xlsx"

# Completed life-history barcodes.
life_history_barcodes = "/Users/johnsk/Documents/Life-History/GLASS-WG/data/sequencing-information/master_life_history_uniform_naming_complete.txt"


#######################################################

# Load necessary packages.
library(tidyverse)
library(maftools)
library(data.table)
library(openxlsx)

#######################################################
# Katie Shao (Novogene provided).
pcawg_linker = readWorkbook(pcawg_id_file, sheet = 1, startRow = 1, colNames = TRUE)
tcga_linker = readWorkbook(tcga_id_file, sheet = 1, startRow = 1, colNames = TRUE)

# Gather the TCGA IDs from our GLASS-WG cohort.
Mutect2dir =  "~/mnt/scratchhelix/johnsk/GLASS_WG_floris/results/mutect2/m2filter/"
# Inspect Mutect2 filters applied to both SNVs and small indels.
setwd(Mutect2dir)

# Create list of names with the ".filtered2.vep_filters.txt". Sample "GLSS-MD-LP03" was removed due to poor quality.
filenames <- list.files(pattern = "*_filters.txt")
list_names <-substr(filenames, 1, 22)
# Link "list_names" object with original TCGA names using barcodes.
vcf_names <- substring(list_names, 1, 15)
life_history_barcode_sheet = read.delim(life_history_barcodes, as.is=T)
# Identify those samples that are used in the GLASS-WG project.
glass_tcga_samples <- life_history_barcode_sheet %>% 
  filter(Barcode%in%vcf_names) %>% 
  filter(grepl("TCGA", Original_ID)) %>% 
  select(Original_ID) %>% 
  .[["Original_ID"]]

#  11 samples from PCAWG (Mutect) are also represent in the GLASS-WG cohort.
pcawg_glass_wg_samples <- life_history_barcode_sheet %>% 
  filter(Original_ID%in%tcga_linker$submitted_sample_id) %>% 
  select(Original_ID) %>% 
  .[["Original_ID"]]

# Retrieve the 11 IDs for variants.
  tcga_linker %>% 
  filter(submitted_sample_id%in%pcawg_glass_wg_samples) %>% 
  inner_join(pcawg_linker, by=c("donor_unique_id" = "donor_unique_id")) %>% 
  inner_join(life_history_barcode_sheet, by=c("submitted_sample_id"="Original_ID")) %>% 
  select(tumor_wgs_aliquot_id, submitted_sample_id, Barcode) %>% 
  write.table(file="/Users/johnsk/Documents/Life-History/PCAWG-GLASS-sample-overlap.txt", sep="\t",
            row.names=FALSE)


