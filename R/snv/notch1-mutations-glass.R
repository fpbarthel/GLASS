##############################################
# Investigation of NOTCH1 mutations
# Updated: 2019.04.09
# Author: Kevin J.
##################################################

# Working directory for this analysis in the GLASS-analysis project. 
mybasedir = "/Volumes/verhaak-lab/GLASS-analysis/"
setwd(mybasedir)

##################################################
# Necessary packages:
library(tidyverse)
library(DBI)
library(Gviz)
library(rtracklayer)
library(trackViewer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

##################################################
# Establish connection with database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB2")

# Load in the necessary information regarding purity
notch1_mutations = dbGetQuery(con, "SELECT * FROM variants.passgeno AS pg
                           LEFT JOIN variants.passanno pa ON pa.variant_id = pg.variant_id
                           WHERE gene_symbol = 'NOTCH1' AND
                           ssm2_pass_call = 'true'")

# Examine the most frequent alterations.
table(notch1_mutations$variant_classification)
# There are some duplicate columns that resulted from the LEFT JOIN, quickly adjust these columns.
colnames(notch1_mutations)[12:16] <- c("variant_id2","chrom2", "pos2", "ref2", "alt2")

# To remove variants that may have been called in both WGS|WXS. Note: I did not use the silver set here.
notch1_unique = notch1_mutations %>% 
  mutate(sample_barcode = substr(aliquot_barcode, 1, 15)) %>% 
  separate(pos, c("tmp", "start", "end", "tmp2")) %>% 
  dplyr::select(sample_barcode, variant_id, variant_classification, variant_type, chrom, start, end, ref2, alt2, gene_symbol, protein_change) %>% 
  distinct() %>% 
  # The below submission to cBioPortal failed for some reason. All required fields seemed to be provided.
  # dplyr::select(Sample_ID	= sample_barcode, Chromosome = chrom, Start_Position = start, End_Position = end, Reference_Allele = ref2, Variant_Allele = alt2)
dplyr::select(Hugo_Symbol = gene_symbol, Protein_Change = protein_change)

# Write out table for upload to cBioPortal. 
write.table(notch1_unique, "/Users/johnsk/Documents/glass_notch1_cbioportal.txt", sep="\t", row.names = F, col.names = T, quote = F)


############ Troubleshooting below this line ##################
# I tried to generate a lollipop plot using different software
# I initially restricted to SNPs because I thought it would be easier.
notch1_unique_snp =  notch1_mutations %>% 
  mutate(sample_barcode = substr(aliquot_barcode, 1, 15)) %>% 
  separate(pos, c("tmp", "start", "end", "tmp2")) %>% 
  dplyr::select(sample_barcode, variant_id, variant_classification, variant_type, chrom, start, end, ref2, alt2, gene_symbol, protein_change) %>% 
  distinct() 
SNP <- c(as.numeric(notch1_unique_snp$start))

# Use summary to get the min and max along the gene.
summary(as.numeric(notch1_unique_snp$start))
# Create a GRanges object with the SNP data.
sample.gr <- GRanges("chr9", IRanges(SNP, width=1, names=paste0("snp", SNP)))
# The commented out code below would help you stack multiple variants on the same locus.
# sample.gr$score <- sample.int(5, length(sample.gr), replace = TRUE)
# The features was the part I struggled with. I couldn't get the exons to show up.
features <- GRanges("chr9", IRanges(c(139389591, 139445147),
                                    width = c(0, 55556),
                                    names = c("start", "end")))
# Not a visually pleasing graph.
lolliplot(sample.gr, features)

# Provided example: https://bioconductor.org/packages/release/bioc/vignettes/trackViewer/inst/doc/trackViewer.html#lolliplot
SNP <- c(10, 12, 1400, 1402)
sample.gr <- GRanges("chr1", IRanges(SNP, width=1, names=paste0("snp", SNP)))
features <- GRanges("chr1", IRanges(c(1, 501, 1001), 
                                    width=c(120, 400, 405),
                                    names=paste0("block", 1:3)))
lolliplot(sample.gr, features)