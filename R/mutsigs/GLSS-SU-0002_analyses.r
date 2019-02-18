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
q = "SELECT * \
FROM analysis.mutsig_fraction mf \
INNER JOIN analysis.silver_set ss ON ss.tumor_pair_barcode = mf.tumor_pair_barcode \
WHERE ss.priority = 1"

mutsig <- dbGetQuery(con,q)

#Get pair information 

q = "SELECT * FROM analysis.tumor_mut_comparison_anno"
pair_info1 <- dbGetQuery(con,q)
pair_info1 <- pair_info1[,-which(colnames(pair_info1)=="surgical_interval_mo")]

q = "SELECT * FROM analysis.clinical_by_tumor_pair"
pair_info2 <- dbGetQuery(con,q)

overlapping_names <- intersect(colnames(pair_info1),colnames(pair_info2))
overlapping_names <- overlapping_names[-1]
pair_info2 <- pair_info2[,-which(colnames(pair_info2) %in% overlapping_names)]

pair_info <- merge(pair_info1, pair_info2, by="tumor_pair_barcode")


#Clean table
#--------------------------------------------------

mutsig_table <- merge(mutsig,pair_info,by="tumor_pair_barcode")
mutsig_table <- mutsig_table[order(mutsig_table[,"tumor_pair_barcode"],mutsig_table[,"signature"],mutsig_table[,"mutation_status"]),]

#Remove WGS
mutsig_table <- mutsig_table[grep("-WXS",mutsig_table[,1]),]

#Tabularize data for heatmap
#--------------------------------------------------

mutsig_table_mod <- mutsig_table
mutsig_table_mod <- mutsig_table_mod[order(nchar(mutsig_table_mod[,"signature"]),mutsig_table_mod[,"signature"]),]
mutsig_table_mod[,"signature"] <- paste(mutsig_table_mod[,"signature"],"_",mutsig_table_mod[,"mutation_status"],sep="")

#Lazy for loop
mutsig_matrix <- matrix(0,nrow=length(unique(mutsig_table_mod[,"tumor_pair_barcode"])),ncol=length(unique(mutsig_table_mod[,"signature"])))
rownames(mutsig_matrix) <- unique(mutsig_table_mod[,"tumor_pair_barcode"])
colnames(mutsig_matrix) <- unique(mutsig_table_mod[,"signature"])
for(i in 1:nrow(mutsig_table_mod))
{
	mutsig_matrix[mutsig_table_mod[i,"tumor_pair_barcode"],mutsig_table_mod[i,"signature"]] <- mutsig_table_mod[i,"relative_contribution"]
}

mysigs <- mutsig_matrix[grep("GLSS-SU-0002",rownames(mutsig_matrix)),]
mysigs_r <- mysigs[grep("_recurrent",names(mysigs))]

colnames(mutsig_matrix) <- gsub("primary","P",colnames(mutsig_matrix))
colnames(mutsig_matrix) <- gsub("recurrent","R",colnames(mutsig_matrix))
colnames(mutsig_matrix) <- gsub("shared","S",colnames(mutsig_matrix))

#Mutations using Samir's DNA repair mutations
mutclass <- read.delim("/projects/varnf/PubDat/Gene_annotations/dna_repair_genes_mdacc_chae_2016.tsv")
mutdata <- dbGetQuery(con, read_file("/projects/verhaak-lab/GLASS-analysis/sql/fred_mutation.sql"))

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

