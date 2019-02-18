#This code makes a heatmap where the columns are private/shared signatures and the rows are samples.
#This visualization is used to identify clusters of patients with unique signature patterns.
#Version 1
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

q = "SELECT * \
FROM analysis.tumor_mut_comparison_anno"

pair_info <- dbGetQuery(con,q)

#Clean table
#--------------------------------------------------

mutsig_table <- merge(mutsig,pair_info,by="tumor_pair_barcode")
mutsig_table <- mutsig_table[order(mutsig_table[,"tumor_pair_barcode"],mutsig_table[,"signature"],mutsig_table[,"mutation_status"]),]

#Remove WGS
mutsig_table <- mutsig_table[grep("-WXS",mutsig_table[,1]),]

#Identify hypermutators
hypermutator_a <- as.numeric(mutsig_table[,"mf_a"] > 10)
hypermutator_b <- as.numeric(mutsig_table[,"mf_b"] > 10)
hypermutator_pair <- as.numeric(hypermutator_a | hypermutator_b)

mutsig_table <- cbind(mutsig_table,hypermutator_a,hypermutator_b,hypermutator_pair)

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

colnames(mutsig_matrix) <- gsub("primary","P",colnames(mutsig_matrix))
colnames(mutsig_matrix) <- gsub("recurrent","R",colnames(mutsig_matrix))
colnames(mutsig_matrix) <- gsub("shared","S",colnames(mutsig_matrix))

pdf("/projects/varnf/GLASS/Figures/signatures/private_shared_sigs_heatmap_full.pdf",width=7,height=5)
Heatmap(mutsig_matrix,
		col = brewer.pal(9,"YlGnBu"),
		cluster_rows=TRUE,
		show_row_names=FALSE,
		cluster_columns=FALSE,
		column_names_gp = gpar(fontsize = 7),
		)
dev.off()

#Subset the matrix
keep_sigs <- c(1,11,12,13,15,16,18,21,23,24,29)
keep_sigs <- paste("Signature.",keep_sigs,"_",sep="")
keep_sigs <- paste(keep_sigs,collapse="|")
sub_mutsig_matrix <- mutsig_matrix[,grep(keep_sigs,colnames(mutsig_matrix),fixed=FALSE)]


#Annotate matrix and add heatmap annotations

#Obtain clinical information
#Get annotations
q = "SELECT su.*, \
    CASE WHEN su.surgery_number > 1 \
    THEN \
        (SELECT treatment_tmz FROM clinical.surgeries WHERE case_barcode = su.case_barcode AND surgery_number = su.surgery_number - 1) \
    ELSE FALSE \
    END AS prior_tmz, \
    CASE WHEN su.surgery_number > 1 \
    THEN \
        (SELECT treatment_radiotherapy FROM clinical.surgeries WHERE case_barcode = su.case_barcode AND surgery_number = su.surgery_number - 1) \
    ELSE FALSE \
    END AS prior_radiation \
FROM clinical.surgeries su \
ORDER BY case_barcode, surgery_number"

clinical_info <- dbGetQuery(con,q)

mutsig_table[,"sample_barcode"] <- sapply(strsplit(mutsig_table[,"tumor_barcode_b"],"-"),function(x)paste(x[1:4],collapse="-"))
mutsig_table_full <- merge(mutsig_table, clinical_info, by = "sample_barcode")

q = "SELECT * FROM analysis.driver_status"
driver_status <- dbGetQuery(con,q)
driver_status <- driver_status[,c("tumor_pair_barcode","driver_status")]
driver_status[,"driver_status"] <- gsub(" ","_",driver_status[,"driver_status"])

mutsig_table_full[,"driver_status"] <- driver_status[match(mutsig_table_full[,"tumor_pair_barcode"],driver_status[,"tumor_pair_barcode"]),"driver_status"]

big_annotation_table <- mutsig_table_full[match(rownames(sub_mutsig_matrix),mutsig_table_full[,"tumor_pair_barcode"]),]
big_annotation_table[,"cohort"] <- sapply(strsplit(big_annotation_table[,"tumor_barcode_a"],"-"),function(x)paste(x[1:2],collapse="-"))

#Age and surgical interval correlation station
q = "SELECT * FROM clinical.cases"
clinical_cases <- dbGetQuery(con,q)
colnames(big_annotation_table)[which(colnames(big_annotation_table)=="case_barcode.x")] <- "case_barcode"
age_table <- merge(big_annotation_table,clinical_cases,by="case_barcode")

surg_cors <- apply(mutsig_matrix,2,function(x)cor(x,big_annotation_table[,"surgical_interval_mo.x"],method="s","pairwise.complete.obs"))
surg_cor_pval <- apply(mutsig_matrix,2,function(x)cor.test(x,big_annotation_table[,"surgical_interval_mo.x"],method="s")$p.value)
surg_cor_qval <- p.adjust(surg_cor_pval,"BH")	
#surg_cor_qval[which(surg_cor_qval < 0.05)]
# Signature.1_S Signature.12_S Signature.16_P Signature.16_S Signature.20_P 
#   0.007965045    0.006337910    0.014478115    0.006337910    0.044119521 
#Signature.21_S 
#   0.014478115 

age_cors <- apply(mutsig_matrix,2,function(x)cor(x,age_table[,"case_age_diagnosis_years"],method="s","pairwise.complete.obs"))
age_cor_pval <- apply(mutsig_matrix,2,function(x)cor.test(x,age_table[,"case_age_diagnosis_years"],method="s")$p.value)
age_cor_qval <- p.adjust(age_cor_pval,"BH")	#Only signature 6 recurrence was significantly negatively correlated with age after correction
#age_cor_qval[which(age_cor_qval < 0.05)]
#Signature.6_R 
#   0.01268544 


#Cluster differences
annotation_table <- big_annotation_table[,c("grade","idh_codel_subtype.x","driver_status","surgical_interval_mo.x","prior_tmz","prior_radiation","hypermutator_pair","cohort")]
colnames(annotation_table) <- c("grade","idh_codel_subtype","driver_status","surgical_interval_mo","prior_tmz","prior_radiation","hypermutator_pair","cohort")

#Set colors for sidebars
myset = brewer.pal(12,"Set3")
sidebar_labels <- apply(annotation_table,2,function(x)unique(x[!is.na(x)])[order(unique(x[!is.na(x)]))])
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
ha = HeatmapAnnotation(df = annotation_table,which="row",
	 col=annotation_colors)

pdf("/projects/varnf/GLASS/Figures/signatures/private_shared_sigs_heatmap_subset.pdf",width=8,height=8)
hm = Heatmap(sub_mutsig_matrix,
		col = brewer.pal(9,"YlGnBu"),
		cluster_rows=TRUE,
		show_row_names=FALSE,
		cluster_columns=FALSE,
		column_names_gp = gpar(fontsize = 8),
		) + ha
hm
dev.off()


#Pull out interesting samples:
draw_hm <- draw(hm)
row_order_hm <- row_order(draw_hm)

clustered_samps <- rownames(sub_mutsig_matrix)[row_order_hm[[1]]]

clustered_big_annotation <- big_annotation_table[match(clustered_samps,big_annotation_table[,"tumor_pair_barcode"]),]
clustered_mutsig_matrix <- sub_mutsig_matrix[clustered_samps,]

#Signature 29-high samples
sig29 <- rownames(clustered_mutsig_matrix)[which(clustered_mutsig_matrix[,"Signature.29_P" ]>0.6)]
clustered_big_annotation[match(sig29,clustered_big_annotation[,"tumor_pair_barcode"]),]
#No obvious similarities

#Signature 15-high samples
sig15 <- rownames(clustered_mutsig_matrix)[which(clustered_mutsig_matrix[,"Signature.15_P" ]>0.4 | clustered_mutsig_matrix[,"Signature.15_R" ]>0.4)]
clustered_big_annotation[match(sig15,clustered_big_annotation[,"tumor_pair_barcode"]),]

sig15r <- rownames(clustered_mutsig_matrix)[which(clustered_mutsig_matrix[,"Signature.15_R" ]>0.4 & clustered_mutsig_matrix[,"Signature.15_P" ]<0.4)]
clustered_big_annotation[match(sig15r,clustered_big_annotation[,"tumor_pair_barcode"]),]
#No obvious similarities

#Mutations in each cluster: 
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