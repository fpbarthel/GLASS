#######################################################
# Map sequenza output to GLASS aliquot barcodes
# Date: 2018.11.06 
# Author: Kevin J.
#######################################################

# Directory for GLASS analysis.
mybasedir = '/Volumes/fastscratch/verhaak-lab/GLASS-WG'
hoon_sequenza_data =  "/Users/johnsk/Documents/GLASS/data/manifest/glass_analysis-paired_bam_map_wes.variant.txt"

#######################################################

# Necessary packages:
library(parallel)
library(tidyverse)
library(data.table)
library(DBI)

#######################################################

# Establish connection with Floris' database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")
extant_aliquots <- dbReadTable(con,  Id(schema="biospecimen",table="aliquots"))


# Read in Hoon's bam sample map that was constructed on a majority of the WXS samples.
glass_wxs_bam_map = read.delim(hoon_sequenza_data, as.is=T)

# For the WXS samples, link together the aliquot_barcode with the sequenza data for quick fix.
database_link = extant_aliquots %>% 
  filter(!grepl("-WGS-", aliquot_barcode)) %>% 
  inner_join(glass_wxs_bam_map, by=c("aliquot_id_legacy"="tm_samplename"))

# Upload file to github so that Floris is able to work with it.
write.table(database_link, file = "/Users/johnsk/Documents/glass-sequenza-link.txt", sep="\t", row.names = F, col.names = T, quote = F)

