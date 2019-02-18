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
q = "SELECT mf.tumor_pair_barcode, mf.fraction, mf.signature, mf.signature_score, mf.mut_count, ctp.*, tp.idh_codel_subtype, cs.grade
FROM analysis.rel_mutsig_fraction mf 
INNER JOIN analysis.silver_set ss ON ss.tumor_pair_barcode = mf.tumor_pair_barcode 
LEFT JOIN clinical.clinical_by_tumor_pair ctp ON mf.tumor_pair_barcode = ctp.tumor_pair_barcode
LEFT JOIN analysis.tumor_mut_comparison_anno tp ON mf.tumor_pair_barcode = tp.tumor_pair_barcode 
LEFT JOIN clinical.surgeries cs ON substring(mf.tumor_pair_barcode, 1, 15) = cs.sample_barcode
WHERE ss.priority = 1"

mutsig_table <- dbGetQuery(con,q)

#Clean table
#--------------------------------------------------

#Remove WGS
mutsig_table <- mutsig_table[grep("-WXS",mutsig_table[,1]),]

#Tabularize data for heatmap
#--------------------------------------------------

mutsig_table_mod <- mutsig_table
mutsig_table_mod <- mutsig_table_mod[order(nchar(mutsig_table_mod[,"signature"]),mutsig_table_mod[,"signature"]),]
mutsig_table_mod[,"signature"] <- paste("Signature",mutsig_table_mod[,"signature"],mutsig_table_mod[,"fraction"],sep="_")
mutsig_table_mod <- mutsig_table_mod[order(mutsig_table_mod[,"signature"]),]

#Lazy for loop
mutsig_matrix <- matrix(0,nrow=length(unique(mutsig_table_mod[,"tumor_pair_barcode"])),ncol=length(unique(mutsig_table_mod[,"signature"])))
rownames(mutsig_matrix) <- unique(mutsig_table_mod[,"tumor_pair_barcode"])
colnames(mutsig_matrix) <- unique(mutsig_table_mod[,"signature"])
for(i in 1:nrow(mutsig_table_mod))
{
	mutsig_matrix[mutsig_table_mod[i,"tumor_pair_barcode"],mutsig_table_mod[i,"signature"]] <- mutsig_table_mod[i,"signature_score"]
}


pdf("/projects/varnf/GLASS/Figures/signatures/SuppFig_heatmap_full_no_tracks.pdf",width=7,height=5)
Heatmap(mutsig_matrix,
		col = brewer.pal(9,"YlGnBu"),
		cluster_rows=TRUE,
		show_row_names=FALSE,
		cluster_columns=FALSE,
		column_names_gp = gpar(fontsize = 7),
		)
dev.off()

#Subset the matrix
keep_sigs <- c(1,3,6,11,12,13,15,16,18,21,24,26)
keep_sigs <- paste("Signature_",keep_sigs,"_",sep="")
keep_sigs <- paste(keep_sigs,collapse="|")
sub_mutsig_matrix <- mutsig_matrix[,grep(keep_sigs,colnames(mutsig_matrix),fixed=FALSE)]
sub_mutsig_matrix <- sub_mutsig_matrix[,order(nchar(colnames(sub_mutsig_matrix)), colnames(sub_mutsig_matrix))]

#Annotate matrix and add heatmap annotations

big_annotation_table <- mutsig_table[match(rownames(sub_mutsig_matrix),mutsig_table[,"tumor_pair_barcode"]),]
big_annotation_table[,"cohort"] <- sapply(strsplit(big_annotation_table[,"tumor_barcode_a"],"-"),function(x)paste(x[1:2],collapse="-"))

#Cluster differences
annotation_table <- big_annotation_table[,c("idh_codel_subtype","received_rt","received_tmz","hypermutator_status","cohort")]

#Set colors for sidebars
myset = brewer.pal(12,"Set3")
sidebar_labels <- apply(annotation_table,2,function(x)unique(x[!is.na(x)])[order(unique(x[!is.na(x)]))])
#sample_type <- c("black","white"); names(sample_type) = sidebar_labels[[1]]
#grade = myset[1:3]; names(grade) = sidebar_labels[[1]]
idh_codel_subtype = c("#F8766D","#00BA38","#619CFF"); names(idh_codel_subtype) = sidebar_labels[[1]]
received_tmz = c("#FFFFFF","#377eb8"); names(received_tmz) = sidebar_labels[[2]]
received_rt = c("#FFFFFF","#377eb8"); names(received_rt) = sidebar_labels[[3]]
hypermutator_status = c("#FFFFFF","#377eb8"); names(hypermutator_status) = sidebar_labels[[4]]
cohort = myset[1:12]; names(cohort) = sidebar_labels[[5]]
annotation_colors = list(idh_codel_subtype=idh_codel_subtype,
					received_rt=received_rt, received_tmz=received_tmz,  
					hypermutator_status=hypermutator_status, cohort=cohort)
ha = HeatmapAnnotation(df = annotation_table,which="row",
	 col=annotation_colors,show_annotation_name=TRUE,show_legend=FALSE,na_col="#E5E5E5")

pdf("/projects/varnf/GLASS/Figures/signatures/SuppFig_heatmap_with_tracks.pdf",width=5.2,height=4.9)
hm = Heatmap(sub_mutsig_matrix,
		col = brewer.pal(9,"YlGnBu"),
		cluster_rows=TRUE,
		show_row_names=FALSE,
		cluster_columns=FALSE,
		column_names_gp = gpar(fontsize = 8),
		show_heatmap_legend=FALSE) + ha
hm
dev.off()

ha = HeatmapAnnotation(df = annotation_table,which="row",
	 col=annotation_colors,show_annotation_name=TRUE,show_legend=TRUE,na_col="#E5E5E5")
pdf("/projects/varnf/GLASS/Figures/signatures/SuppFig_heatmap_legend.pdf",width=7,height=7)
hm = Heatmap(sub_mutsig_matrix,
		col = brewer.pal(9,"YlGnBu"),
		cluster_rows=TRUE,
		show_row_names=FALSE,
		cluster_columns=FALSE,
		column_names_gp = gpar(fontsize = 8),
		show_heatmap_legend=TRUE) + ha
hm
dev.off()