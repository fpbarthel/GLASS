#This code identifies signatures and mutations specific to GLSS-SU-0002-R3, a hypermutant sample that was not treated with TMZ.
#Author: Fred Varn
#--------------------------------------------------

library(DBI)
library(odbc)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(tidyverse)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

#Download
#--------------------------------------------------
#Get mutational signature info
q = "SELECT mf.tumor_pair_barcode, tp.tumor_barcode_a, tp.tumor_barcode_b, mf.fraction, mf.signature, mf.abs_score, mf.rel_score, mf.mut_n, tp.idh_codel_subtype,cc.received_alk,cc.received_rt,cc.hypermutator_status
FROM analysis.mut_sig_fraction mf 
INNER JOIN analysis.silver_set ss ON ss.tumor_pair_barcode = mf.tumor_pair_barcode 
LEFT JOIN analysis.tumor_mut_comparison_anno tp ON mf.tumor_pair_barcode = tp.tumor_pair_barcode
LEFT JOIN analysis.tumor_clinical_comparison cc ON mf.tumor_pair_barcode = cc.tumor_pair_barcode"

mutsig <- dbGetQuery(con,q)

#Remove WGS
mutsig_table <- mutsig_table[grep("-WXS",mutsig_table[,1]),]

#Tabularize data for heatmap
#--------------------------------------------------

mutsig_table_mod <- mutsig_table
mutsig_table_mod <- mutsig_table_mod[order(nchar(mutsig_table_mod[,"signature"]),mutsig_table_mod[,"signature"]),]
mutsig_table_mod[,"signature"] <- paste(mutsig_table_mod[,"signature"],"_",mutsig_table_mod[,"fraction"],sep="")

#Lazy for loop
mutsig_matrix <- matrix(0,nrow=length(unique(mutsig_table_mod[,"tumor_pair_barcode"])),ncol=length(unique(mutsig_table_mod[,"signature"])))
rownames(mutsig_matrix) <- unique(mutsig_table_mod[,"tumor_pair_barcode"])
colnames(mutsig_matrix) <- unique(mutsig_table_mod[,"signature"])
for(i in 1:nrow(mutsig_table_mod))
{
	mutsig_matrix[mutsig_table_mod[i,"tumor_pair_barcode"],mutsig_table_mod[i,"signature"]] <- mutsig_table_mod[i,"rel_score"]
}

mysigs <- mutsig_matrix[grep("GLSS-SU-0002",rownames(mutsig_matrix)),]
mysigs_r <- mysigs[grep("_R",names(mysigs))]

#Mutations using Samir's DNA repair mutations
mutclass <- read.delim("/projects/varnf/PubDat/Gene_annotations/dna_repair_genes_mdacc_chae_2016.tsv")
mutdata <- dbGetQuery(con, read_file("/projects/verhaak-lab/GLASS-analysis/sql/fred_mutation2.sql"))

mmr <- as.character(mutclass[grep("MMR",mutclass[,"dna_repair_pathways"]),"gene_symbol"])

single_sample_mut <- mutdata[grep("GLSS-SU-0002",mutdata[,"case_barcode"]),]
genes <- single_sample_mut[,"gene_symbol"]

mutclass[which(mutclass[,"gene_symbol"] %in% genes),]

#Check for the presence of mutations in each gene in the dataset
sub_mutclass <- mutdata[which(mutdata[,"gene_symbol"]%in%genes),]

rle(as.character(sub_mutclass[,"gene_symbol"]))

hypermutators <- pair_info[which(pair_info[,"hypermutator_status"]==1 & pair_info[,"received_tmz"]==1),"case_barcode"]


sub_mutclass <- sub_mutclass[-which(sub_mutclass[,"case_barcode"]%in%hypermutators),]
rle(as.character(sub_mutclass[,"gene_symbol"]))

sub_mutclass[,c(1,2,3,ncol(sub_mutclass))]

