#This code makes a heatmap where the columns are private/shared signatures and the rows are samples.
#This visualization is used to identify clusters of patients with unique signature patterns.
#Version 2: Uses revised tumor pair data and includes analysis of DNA repair pathways
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
q = "SELECT mf.tumor_pair_barcode, mf.fraction, mf.signature, mf.signature_score, mf.mut_count \
FROM analysis.rel_mutsig_fraction mf \
INNER JOIN analysis.silver_set ss ON ss.tumor_pair_barcode = mf.tumor_pair_barcode \
WHERE ss.priority = 1"

mutsig <- dbGetQuery(con,q)

#Get pair information 

q = "SELECT ctp.*, tp.idh_codel_subtype, cs.grade \
FROM analysis.clinical_by_tumor_pair ctp \
INNER JOIN analysis.tumor_mut_comparison_anno tp ON tp.tumor_pair_barcode = ctp.tumor_pair_barcode \
INNER JOIN clinical.surgeries cs ON cs.sample_barcode = substring(ctp.tumor_barcode_a from 1 for 15)"

pair_info <- dbGetQuery(con,q)


#Clean table
#--------------------------------------------------

mutsig_table <- merge(mutsig,pair_info,by="tumor_pair_barcode")
mutsig_table <- mutsig_table[order(mutsig_table[,"tumor_pair_barcode"],mutsig_table[,"signature"],mutsig_table[,"fraction"]),]

#Remove WGS
mutsig_table <- mutsig_table[grep("-WXS",mutsig_table[,1]),]


#Tabularize data for heatmap
#--------------------------------------------------

mutsig_table_mod <- mutsig_table
mutsig_table_mod <- mutsig_table_mod[order(nchar(mutsig_table_mod[,"signature"]),mutsig_table_mod[,"signature"]),]
mutsig_table_mod[,"signature"] <- paste("Signature",mutsig_table_mod[,"signature"],mutsig_table_mod[,"fraction"],sep="_")

#Lazy for loop
mutsig_matrix <- matrix(0,nrow=length(unique(mutsig_table_mod[,"tumor_pair_barcode"])),ncol=length(unique(mutsig_table_mod[,"signature"])))
rownames(mutsig_matrix) <- unique(mutsig_table_mod[,"tumor_pair_barcode"])
colnames(mutsig_matrix) <- unique(mutsig_table_mod[,"signature"])
for(i in 1:nrow(mutsig_table_mod))
{
	mutsig_matrix[mutsig_table_mod[i,"tumor_pair_barcode"],mutsig_table_mod[i,"signature"]] <- mutsig_table_mod[i,"signature_score"]
}


pdf("/projects/varnf/GLASS/Figures/signatures/private_shared_sigs_heatmap_full_v3.pdf",width=7,height=5)
Heatmap(mutsig_matrix,
		col = brewer.pal(9,"YlGnBu"),
		cluster_rows=TRUE,
		show_row_names=FALSE,
		cluster_columns=FALSE,
		column_names_gp = gpar(fontsize = 7),
		)
dev.off()

#Subset the matrix
keep_sigs <- c(1,3,11,12,13,15,18,21,23,24,29)
keep_sigs <- paste("Signature_",keep_sigs,"_",sep="")
keep_sigs <- paste(keep_sigs,collapse="|")
sub_mutsig_matrix <- mutsig_matrix[,grep(keep_sigs,colnames(mutsig_matrix),fixed=FALSE)]


#Annotate matrix and add heatmap annotations

big_annotation_table <- mutsig_table[match(rownames(sub_mutsig_matrix),mutsig_table[,"tumor_pair_barcode"]),]
big_annotation_table[,"cohort"] <- sapply(strsplit(big_annotation_table[,"tumor_barcode_a"],"-"),function(x)paste(x[1:2],collapse="-"))

#Cluster differences
annotation_table <- big_annotation_table[,c("grade","idh_codel_subtype","received_tmz","received_rt","hypermutator_status","cohort")]

#Set colors for sidebars
myset = brewer.pal(12,"Set3")
sidebar_labels <- apply(annotation_table,2,function(x)unique(x[!is.na(x)])[order(unique(x[!is.na(x)]))])
#sample_type <- c("black","white"); names(sample_type) = sidebar_labels[[1]]
grade = myset[1:3]; names(grade) = sidebar_labels[[1]]
idh_codel_subtype = myset[4:6]; names(idh_codel_subtype) = sidebar_labels[[2]]
received_tmz = c("white","black"); names(received_tmz) = sidebar_labels[[3]]
received_rt = c("white","black"); names(received_rt) = sidebar_labels[[4]]
hypermutator_status = c("white","black"); names(hypermutator_status) = sidebar_labels[[5]]
cohort = myset[1:12]; names(cohort) = sidebar_labels[[6]]
annotation_colors = list(grade=grade, idh_codel_subtype=idh_codel_subtype,
					received_tmz=received_tmz, received_rt=received_rt, 
					hypermutator_status=hypermutator_status, cohort=cohort)
ha = HeatmapAnnotation(df = annotation_table,which="row",
	 col=annotation_colors)

pdf("/projects/varnf/GLASS/Figures/signatures/private_shared_sigs_heatmap_subset_v3.pdf",width=7,height=7)
hm = Heatmap(sub_mutsig_matrix,
		col = brewer.pal(9,"YlGnBu"),
		cluster_rows=TRUE,
		show_row_names=FALSE,
		cluster_columns=FALSE,
		column_names_gp = gpar(fontsize = 8),
		) + ha
hm
dev.off()

#Identifying characteristics specific to different clusters
#--------------------------------------------------

#Pull out interesting samples:
draw_hm <- draw(hm)
row_order_hm <- row_order(draw_hm)

clustered_samps <- rownames(sub_mutsig_matrix)[row_order_hm[[1]]]

clustered_big_annotation <- big_annotation_table[match(clustered_samps,big_annotation_table[,"tumor_pair_barcode"]),]
clustered_mutsig_matrix <- sub_mutsig_matrix[clustered_samps,]

#Clinical information
#----------------------
#Signature 29-high samples
sig29 <- rownames(clustered_mutsig_matrix)[which(clustered_mutsig_matrix[,"Signature_29_P" ]>0.6)]
clustered_big_annotation[match(sig29,clustered_big_annotation[,"tumor_pair_barcode"]),]
#No obvious similarities

#Signature 15-high samples
sig15 <- rownames(clustered_mutsig_matrix)[which(clustered_mutsig_matrix[,"Signature_15_P" ]>0.4 | clustered_mutsig_matrix[,"Signature_15_R" ]>0.4)]
clustered_big_annotation[match(sig15,clustered_big_annotation[,"tumor_pair_barcode"]),]

sig15r <- rownames(clustered_mutsig_matrix)[which(clustered_mutsig_matrix[,"Signature_15_R" ]>0.4 & clustered_mutsig_matrix[,"Signature_15_P" ]<0.4)]
clustered_big_annotation[match(sig15r,clustered_big_annotation[,"tumor_pair_barcode"]),]
#No obvious similarities

sig3 <- rownames(clustered_mutsig_matrix)[which(clustered_mutsig_matrix[,"Signature_3_R" ]>0)]
clustered_big_annotation[match(sig3,clustered_big_annotation[,"tumor_pair_barcode"]),]

#Mutations using Samir's DNA repair genes
#----------------------
mutclass <- read.delim("/projects/varnf/PubDat/Gene_annotations/dna_repair_genes_mdacc_chae_2016.tsv")
mutdata <- dbGetQuery(con, read_file("/projects/verhaak-lab/GLASS-analysis/sql/fred_mutation.sql"))

mmr <- as.character(mutclass[grep("MMR",mutclass[,"dna_repair_pathways"]),"gene_symbol"])
ner <- as.character(mutclass[grep("NER",mutclass[,"dna_repair_pathways"]),"gene_symbol"])
hr <- as.character(mutclass[grep("HR",mutclass[,"dna_repair_pathways"]),"gene_symbol"])

mmrdata <- mutdata[which(mutdata[,"gene_symbol"] %in% mmr),]
nerdata <- mutdata[which(mutdata[,"gene_symbol"] %in% ner),]
hrdata <- mutdata[which(mutdata[,"gene_symbol"] %in% hr),]

sig29_case_barcode <- sapply(strsplit(sig29,"-"),function(x)paste(x[1:3],collapse="-"))
sig29_muts <- nerdata[which(nerdata[,"case_barcode"] %in% sig29_case_barcode),]
sig29_muts <- sig29_muts[order(sig29_muts[,"gene_symbol"]),]		#Nothing obvious here

sig15_case_barcode <- sapply(strsplit(sig15,"-"),function(x)paste(x[1:3],collapse="-"))
sig15_muts <- mmrdata[which(mmrdata[,"case_barcode"] %in% sig15_case_barcode),]
sig15_muts <- sig15_muts[order(sig15_muts[,"gene_symbol"]),]		#Nothing obvious here

sig15r_case_barcode <- sapply(strsplit(sig15r,"-"),function(x)paste(x[1:3],collapse="-"))
sig15r_muts <- mmrdata[which(mmrdata[,"case_barcode"] %in% sig15r_case_barcode),]
sig15r_muts <- sig15r_muts[order(sig15r_muts[,"gene_symbol"]),]		

sig3_case_barcode <- sapply(strsplit(sig3,"-"),function(x)paste(x[1:3],collapse="-"))
sig3_muts <- hrdata[which(hrdata[,"case_barcode"] %in% sig3_case_barcode),]
sig3_muts <- sig3_muts[order(sig3_muts[,"gene_symbol"]),]		#Nothing obvious here

tmp <- mutsig_matrix[order(mutsig_matrix[,"Signature_3_R"],decreasing=T),]

#Copy number using Samir's DNA repair genes
#----------------------
mutclass <- read.delim("/projects/varnf/PubDat/Gene_annotations/dna_repair_genes_mdacc_chae_2016.tsv")
cndata <- dbGetQuery(con, read_file("/projects/verhaak-lab/GLASS-analysis/sql/fred_cnv.sql"))
cndata <- cndata[-which(cndata[,"cnv_class"]=="Heterozygous"),]

mmr <- as.character(mutclass[grep("MMR",mutclass[,"dna_repair_pathways"]),"gene_symbol"])
ner <- as.character(mutclass[grep("NER",mutclass[,"dna_repair_pathways"]),"gene_symbol"])

mmrdata <- cndata[which(cndata[,"gene_symbol"] %in% mmr),]
nerdata <- cndata[which(cndata[,"gene_symbol"] %in% ner),]

sig29_case_barcode <- sapply(strsplit(sig29,"-"),function(x)paste(x[1:3],collapse="-"))
sig29_muts <- nerdata[which(nerdata[,"case_barcode"] %in% sig29_case_barcode),]
sig29_muts <- sig29_muts[order(sig29_muts[,"gene_symbol"]),]		#Nothing obvious here

sig15_case_barcode <- sapply(strsplit(sig15,"-"),function(x)paste(x[1:3],collapse="-"))
sig15_muts <- mmrdata[which(mmrdata[,"case_barcode"] %in% sig15_case_barcode),]
sig15_muts <- sig15_muts[order(sig15_muts[,"gene_symbol"]),]		#Nothing obvious here

sig15r_case_barcode <- sapply(strsplit(sig15r,"-"),function(x)paste(x[1:3],collapse="-"))
sig15r_muts <- mmrdata[which(mmrdata[,"case_barcode"] %in% sig15r_case_barcode),]
sig15r_muts <- sig15r_muts[order(sig15r_muts[,"gene_symbol"]),]		


#Driver mutations in each cluster: 
#----------------------
mutdata <- dbGetQuery(con, read_file("/projects/verhaak-lab/GLASS-analysis/sql/build_heatmap_data_mutation.sql"))

sig29_case_barcode <- sapply(strsplit(sig29,"-"),function(x)paste(x[1:3],collapse="-"))
sig29_muts <- mutdata[which(mutdata[,"case_barcode"] %in% sig29_case_barcode),]
sig29_muts <- sig29_muts[order(sig29_muts[,"gene_symbol"]),]		#Nothing obvious here

sig15_case_barcode <- sapply(strsplit(sig15,"-"),function(x)paste(x[1:3],collapse="-"))
sig15_muts <- mutdata[which(mutdata[,"case_barcode"] %in% sig15_case_barcode),]
sig15_muts <- sig15_muts[order(sig15_muts[,"gene_symbol"]),]		#Nothing obvious here

sig15r_case_barcode <- sapply(strsplit(sig15r,"-"),function(x)paste(x[1:3],collapse="-"))
sig15r_muts <- mutdata[which(mutdata[,"case_barcode"] %in% sig15r_case_barcode),]
sig15r_muts <- sig15r_muts[order(sig15r_muts[,"gene_symbol"]),]		
#Of slight interest: 2/4 samples had PTEN mutations that were shed at recurrence

#Driver copy number in each cluster:
#----------------------
cndata <- dbGetQuery(con, read_file("/projects/verhaak-lab/GLASS-analysis/sql/build_heatmap_data_cnv.sql"))

sig29_case_barcode <- sapply(strsplit(sig29,"-"),function(x)paste(x[1:3],collapse="-"))
sig29_cnas <- cndata[which(cndata[,"case_barcode"] %in% sig29_case_barcode),]
sig29_cnas <- sig29_cnas[order(sig29_cnas[,"gene_symbol"]),]		#Nothing obvious here

sig15_case_barcode <- sapply(strsplit(sig15,"-"),function(x)paste(x[1:3],collapse="-"))
sig15_cnas <- cndata[which(cndata[,"case_barcode"] %in% sig15_case_barcode),]
sig15_cnas <- sig15_cnas[order(sig15_cnas[,"gene_symbol"]),]		#Nothing obvious here

sig15r_case_barcode <- sapply(strsplit(sig15r,"-"),function(x)paste(x[1:3],collapse="-"))
sig15r_cnas <- cndata[which(cndata[,"case_barcode"] %in% sig15r_case_barcode),]
sig15r_cnas <- sig15r_cnas[order(sig15r_cnas[,"gene_symbol"]),]		
#Nothing obvious

#--------------------------------------------------

#Sanity check: Make a heatmap that is only primaries and R1 (eliminates false clusters due to multiple recurrences)

#Subset the matrix further 
TP_R1_big_annotation_table <- big_annotation_table[which(big_annotation_table[,"sample_type_a"]=="TP" & big_annotation_table[,"sample_type_b"]=="R1"),]
TP_R1_mutsig_matrix <- sub_mutsig_matrix[TP_R1_big_annotation_table[,"tumor_pair_barcode"],]

#Cluster differences
TP_R1_annotation_table <- TP_R1_big_annotation_table[,c("grade","idh_codel_subtype.x","driver_status","surgical_interval_mo.x","prior_tmz","prior_radiation","hypermutator_pair","cohort")]
colnames(TP_R1_annotation_table) <- c("grade","idh_codel_subtype","driver_status","surgical_interval_mo","prior_tmz","prior_radiation","hypermutator_pair","cohort")

#Set colors for sidebars
myset = brewer.pal(12,"Set3")
sidebar_labels <- apply(TP_R1_annotation_table,2,function(x)unique(x[!is.na(x)])[order(unique(x[!is.na(x)]))])
#sample_type <- c("black","white"); names(sample_type) = sidebar_labels[[1]]
grade = myset[1:3]; names(grade) = sidebar_labels[[1]]
idh_codel_subtype = myset[4:6]; names(idh_codel_subtype) = sidebar_labels[[2]]
driver_status = myset[7:10]; names(driver_status) = sidebar_labels[[3]]
prior_tmz = c("white","black"); names(prior_tmz) = sidebar_labels[[5]]
prior_radiation = c("white","black"); names(prior_radiation) = sidebar_labels[[6]]
hypermutator_pair = c("white","black"); names(hypermutator_pair) = sidebar_labels[[7]]
cohort = myset[1:12]; names(cohort) = sidebar_labels[[8]]
annotation_colors = list(grade=grade, idh_codel_subtype=idh_codel_subtype,
					driver_status=driver_status,
					prior_tmz=prior_tmz, prior_radiation=prior_radiation, 
					hypermutator_pair=hypermutator_pair, cohort=cohort)
ha = HeatmapAnnotation(df = TP_R1_annotation_table,which="row",
	 col=annotation_colors)
	 
pdf("/projects/varnf/GLASS/Figures/signatures/private_shared_sigs_heatmap_TP_R1.pdf",width=8,height=8)
hm = Heatmap(TP_R1_mutsig_matrix,
		col = brewer.pal(9,"YlGnBu"),
		cluster_rows=TRUE,
		show_row_names=FALSE,
		cluster_columns=FALSE,
		column_names_gp = gpar(fontsize = 8),
		) + ha
hm
dev.off()