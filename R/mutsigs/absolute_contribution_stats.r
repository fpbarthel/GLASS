library(DBI)
library(odbc)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(tidyverse)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

myinf1 <- "/projects/varnf/GLASS/analysis/signatures/absolute_contributions_by_sample.txt"

absolute_mutsig <- t(read.delim(myinf1,row.names=1))
rownames(absolute_mutsig) <- gsub("\\.","-",rownames(absolute_mutsig))
absolute_mutsig <- absolute_mutsig[-grep("-WGS-",rownames(absolute_mutsig)),]

#Age and surgical interval correlation station
q = "SELECT * FROM clinical.cases"
clinical_cases <- dbGetQuery(con,q)

#Primary age analysis
primary_mutsig_matrix <- absolute_mutsig[grep("-TP-",rownames(absolute_mutsig)),]
primary_barcode <- sapply(strsplit(rownames(primary_mutsig_matrix),"-"),function(x)paste(x[1:3],collapse="-"))
ordered_age <- clinical_cases[match(primary_barcode,clinical_cases[,"case_barcode"]),"case_age_diagnosis_years"]

age_cor_coef <- apply(primary_mutsig_matrix,2,function(x)cor(x,ordered_age,method="s"))
age_cor_pval <- apply(primary_mutsig_matrix,2,function(x)cor.test(x,ordered_age,method="s")$p.value)
age_cor_qval <- p.adjust(age_cor_pval,"BH")

age_sig_res <- data.frame(age_cor_coef,age_cor_pval,age_cor_qval)
age_sig_res <- age_sig_res[which(age_cor_qval<0.05),]
colnames(age_sig_res) <- c("coef","pval","qval")

#-----------------------------------------------

myinf1 <- "/projects/varnf/GLASS/analysis/signatures/absolute_contributions_by_private_shared.txt"

absolute_mutsig <- t(read.delim(myinf1,row.names=1))
rownames(absolute_mutsig) <- gsub("\\.","-",rownames(absolute_mutsig))
absolute_mutsig <- absolute_mutsig[-grep("-WGS-",rownames(absolute_mutsig)),]

#Age and surgical interval correlation station
q = "SELECT * FROM analysis.clinical_by_tumor_pair"
pair_info <- dbGetQuery(con,q)

recurrence_mutsig_matrix <- absolute_mutsig[grep("-recurrent",rownames(absolute_mutsig)),]
rownames(recurrence_mutsig_matrix) <- sapply(strsplit(rownames(recurrence_mutsig_matrix),"-"),function(x)paste(x[1:8],collapse="-"))

#Pull out R1 vs TP samples only
#sample_type_a <- sapply(strsplit(pair_info[,"tumor_barcode_a"],"-"),function(x)x[4])
#sample_type_b <- sapply(strsplit(pair_info[,"tumor_barcode_b"],"-"),function(x)x[4])
#mypairs <- pair_info[which(sample_type_a=="TP" & sample_type_b=="R1"),"tumor_pair_barcode"]
#mypairs <- intersect(mypairs,rownames(recurrence_mutsig_matrix))
#recurrence_mutsig_matrix <- recurrence_mutsig_matrix[mypairs,]

#Continuous correlations
#------------------------------------
#Recurrence surgical interval analysis
ordered_surg <- pair_info[match(rownames(recurrence_mutsig_matrix),pair_info[,"tumor_pair_barcode"]),"surgical_interval"]

surg_cor_coef <- apply(recurrence_mutsig_matrix,2,function(x)cor(x,ordered_surg,method="s","pairwise.complete.obs"))
surg_cor_pval <- apply(recurrence_mutsig_matrix,2,function(x)cor.test(x,ordered_surg,method="s")$p.value)
surg_cor_qval <- p.adjust(surg_cor_pval,"BH")

surg_sig_res <- data.frame(surg_cor_coef,surg_cor_pval,surg_cor_qval)
surg_sig_res <- surg_sig_res[which(surg_cor_qval<0.05),]
colnames(surg_sig_res) <- c("coef","pval","qval")					#0 correlated with surgical interval

#Recurrence radiation interval analysis
ordered_rad <- pair_info[match(rownames(recurrence_mutsig_matrix),pair_info[,"tumor_pair_barcode"]),"received_rt_sum_gy"]

rad_cor_coef <- apply(recurrence_mutsig_matrix,2,function(x)cor(x,ordered_rad,method="s","pairwise.complete.obs"))
rad_cor_pval <- apply(recurrence_mutsig_matrix,2,function(x)cor.test(x,ordered_rad,method="s")$p.value)
rad_cor_qval <- p.adjust(rad_cor_pval,"BH")

rad_sig_res <- data.frame(rad_cor_coef,rad_cor_pval,rad_cor_qval)
rad_sig_res <- rad_sig_res[which(rad_cor_qval<0.05),]
colnames(rad_sig_res) <- c("coef","pval","qval")					#0 correlated with surgical interval

#Clinical variable significance testing
#------------------------------------

vars <- c("received_tmz","received_rt","hypermutator_status")
var_res_list <- list()
for(i in 1:length(vars))
{
	ordered_var <- pair_info[match(rownames(recurrence_mutsig_matrix),pair_info[,"tumor_pair_barcode"]),vars[i]]
	
	g1 <- pair_info[which(pair_info[,vars[i]]==1),"tumor_pair_barcode"]
	g2 <- pair_info[which(pair_info[,vars[i]]==0),"tumor_pair_barcode"]
	
	eff <- apply(recurrence_mutsig_matrix, 2, function(x)median(x[g1],na.rm=T) - median(x[g2],na.rm=T))
	pval <- apply(recurrence_mutsig_matrix, 2, function(x)wilcox.test(x[g1], x[g2])$p.val)
	qval <- p.adjust(pval,"BH")
	
	res <- data.frame(eff,pval,qval)
	res <- res[which(res[,"qval"] < 0.05),]
	var_res_list[[i]] <- res
}
names(var_res_list) <- vars											#only hypermutators exhibited enrichment for a signature (3, 11, 13, 16, 23)

#TMZ-associated hypermutation:
ordered_var <- pair_info[match(rownames(recurrence_mutsig_matrix),pair_info[,"tumor_pair_barcode"]),]
g1 <- pair_info[which(pair_info[,"received_tmz"]==1 & pair_info[,"hypermutator_status"]==1),"tumor_pair_barcode"]
g2 <- pair_info[which(pair_info[,"received_tmz"]==0 | pair_info[,"hypermutator_status"]==0),"tumor_pair_barcode"]

eff <- apply(recurrence_mutsig_matrix, 2, function(x)mean(x[g1],na.rm=T) - mean(x[g2],na.rm=T))
pval <- apply(recurrence_mutsig_matrix, 2, function(x)wilcox.test(x[g1], x[g2])$p.val)
qval <- p.adjust(pval,"BH")

#Grade change comparisons removing hypermutators
#------------------------------------
#pair_info <- pair_info[-which(pair_info[,"hypermutator_status"]==1),]
ordered_var <- pair_info[match(rownames(recurrence_mutsig_matrix), pair_info[,"tumor_pair_barcode"]),"grade_change"]
g1 <- pair_info[which(pair_info[,"grade_change"]=="Grade up"),"tumor_pair_barcode"]
g2 <- pair_info[which(pair_info[,"grade_change"]=="Grade stable"),"tumor_pair_barcode"]

g1 <- intersect(g1,rownames(recurrence_mutsig_matrix))
g2 <- intersect(g2,rownames(recurrence_mutsig_matrix))

grade_eff <- apply(recurrence_mutsig_matrix, 2, function(x)mean(x[g1],na.rm=T) - mean(x[g2],na.rm=T))
grade_pval <- apply(recurrence_mutsig_matrix, 2, function(x)wilcox.test(x[g1], x[g2])$p.val)
grade_qval <- p.adjust(grade_pval,"BH")

grade_sig_res <- data.frame(grade_eff,grade_pval,grade_qval)
grade_sig_res <- grade_sig_res[which(grade_qval<0.05),]




#-----------------------------------------------------------------

#Molecular absolute contribution analyses

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
FROM analysis.abs_mutsig_fraction mf \
INNER JOIN analysis.silver_set ss ON ss.tumor_pair_barcode = mf.tumor_pair_barcode \
WHERE ss.priority = 1"

mutsig <- dbGetQuery(con,q)

mutdata <- dbGetQuery(con, read_file("/projects/verhaak-lab/GLASS-analysis/sql/fred_mutation.sql"))

mutsig_table <- mutsig
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

mutsig_matrix <- mutsig_matrix[order(mutsig_matrix[,"Signature_3_R"],decreasing=T),]

#Ordered signature 29 recurrence and primary only: Likely artifact?
#-----------------------
#> tmp <- mutsig_matrix[order(mutsig_matrix[,"Signature_29_R"],decreasing=T),]
#> head(tmp[,grep("29",colnames(tmp))])
#                              Signature_29_P Signature_29_R Signature_29_S
#GLSS-AT-00P3-TP-01-R1-01D-WXS      137.29669       257.8471      0.8509632
#GLSS-AT-00P1-TP-01-R1-01D-WXS      134.57414       179.5338      2.1236601
#GLSS-AT-00P7-TP-01-R1-01D-WXS      163.28447       154.4761      0.0000000
#GLSS-AT-00P5-TP-01-R1-01D-WXS      140.49809       141.4013      0.0000000
#GLSS-AT-00P6-TP-01-R1-01D-WXS       67.22821       130.3385      0.0000000
#GLSS-AT-00P2-TP-01-R1-01D-WXS      124.54295       110.1184      0.0000000
#> tmp[1:10,grep("29",colnames(tmp))]
#                              Signature_29_P Signature_29_R Signature_29_S
#GLSS-AT-00P3-TP-01-R1-01D-WXS     137.296690      257.84710      0.8509632
#GLSS-AT-00P1-TP-01-R1-01D-WXS     134.574140      179.53383      2.1236601
#GLSS-AT-00P7-TP-01-R1-01D-WXS     163.284473      154.47613      0.0000000
#GLSS-AT-00P5-TP-01-R1-01D-WXS     140.498093      141.40132      0.0000000
#GLSS-AT-00P6-TP-01-R1-01D-WXS      67.228214      130.33853      0.0000000
#GLSS-AT-00P2-TP-01-R1-01D-WXS     124.542945      110.11839      0.0000000
#GLSS-AT-00P4-TP-01-R1-01D-WXS      83.465488      104.34940      0.0000000
#GLSS-SM-R072-TP-01-R1-01D-WXS      36.744584       49.14541      3.4344817
#GLSS-MD-0094-TP-01-R1-01D-WXS       0.000000       38.43419      1.3214988
#GLSS-SF-0163-TP-01-R1-01D-WXS       3.178507       37.91446      0.0000000
#> tmp <- mutsig_matrix[order(mutsig_matrix[,"Signature_29_P"],decreasing=T),]
#> tmp[1:10,grep("29",colnames(tmp))]
#                              Signature_29_P Signature_29_R Signature_29_S
#GLSS-AT-00P7-TP-01-R1-01D-WXS      163.28447      154.47613      0.0000000
#GLSS-AT-00P5-TP-01-R1-01D-WXS      140.49809      141.40132      0.0000000
#GLSS-AT-00P3-TP-01-R1-01D-WXS      137.29669      257.84710      0.8509632
#GLSS-AT-00P1-TP-01-R1-01D-WXS      134.57414      179.53383      2.1236601
#GLSS-AT-00P2-TP-01-R1-01D-WXS      124.54295      110.11839      0.0000000
#GLSS-AT-00P4-TP-01-R1-01D-WXS       83.46549      104.34940      0.0000000
#GLSS-AT-00P6-TP-01-R1-01D-WXS       67.22821      130.33853      0.0000000
#GLSS-SF-0007-TP-01-R1-01D-WXS       36.92031       20.40166      0.2219343
#GLSS-SM-R072-TP-01-R1-01D-WXS       36.74458       49.14541      3.4344817
#GLSS-CU-R007-TP-01-R1-01D-WXS       33.07459        0.00000      0.0000000

#Clinical information
#----------------------
#Signature 29-high samples
sig29 <- rownames(mutsig_matrix)[which(mutsig_matrix[,"Signature_29_P" ]>0.6)]
sig15 <- rownames(mutsig_matrix)[which(mutsig_matrix[,"Signature_15_P" ]>0.4 | mutsig_matrix[,"Signature_15_R" ]>0.4)]
sig15r <- rownames(mutsig_matrix)[which(mutsig_matrix[,"Signature_15_R" ]>0.4 & mutsig_matrix[,"Signature_15_P" ]<0.4)]
sig3 <- rownames(mutsig_matrix)[which(mutsig_matrix[,"Signature_3_R" ]>0)]

#Mutations using Samir's DNA repair genes
#----------------------
mutclass <- read.delim("/projects/varnf/PubDat/Gene_annotations/dna_repair_genes_mdacc_chae_2016.tsv")

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

mmrdata <- cndata[which(cndata[,"gene_symbol"] %in% mmr),]
nerdata <- cndata[which(cndata[,"gene_symbol"] %in% ner),]
hrdata <- cndata[which(cndata[,"gene_symbol"] %in% hr),]
hrdata <- hrdata[-which(is.na(hrdata[,"cnv_retention"])),]

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
