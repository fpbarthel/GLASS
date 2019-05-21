library(DBI)
library(odbc)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(tidyverse)
library(reshape)
library(gridExtra)

rm(list=ls())

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

#Download
#--------------------------------------------------
neo <- dbGetQuery(con,read_file("/projects/verhaak-lab/GLASS-analysis/sql/neoag_with_norm_contam_pair.sql"))
neonew <- neo

pairs <- unique(neo[,"tumor_pair_barcode"])

n_primary <- n_recurrent <- n_shared <- med_primary <- med_recurrent <- med_shared <- surg_int <- norm_contam_a <- norm_contam_b <- hypermutator_recur <- tmz <- rt <- rep(0,length(pairs))
wilcox_ps <- wilcox_pr <- wilcox_rs <- rep(NA,length(pairs))
tumor_barcode_a <- tumor_barcode_b <- subtype <- immunotherapy <- rep("",length(pairs))
for(i in 1:length(pairs))
{
	sub_neo <- neo[which(neo[,"tumor_pair_barcode"]==pairs[i]),]
	
	primary_score <- sub_neo[which(sub_neo[,"fraction"]=="P"),"netmhcpan_mt_score"]
	recurrent_score <- sub_neo[which(sub_neo[,"fraction"]=="R"),"netmhcpan_mt_score"]
	shared_score <- sub_neo[which(sub_neo[,"fraction"]=="S"),"netmhcpan_mt_score"]
	
	n_primary[i] <- length(primary_score)
	n_recurrent[i] <- length(recurrent_score)
	n_shared[i] <- length(shared_score)
	
	med_primary[i] <- median(primary_score)
	med_recurrent[i] <- median(recurrent_score)
	med_shared[i] <- median(shared_score)
		
	if(length(primary_score)>0 & length(shared_score)>0){
		wilcox_ps[i] <- wilcox.test(primary_score,shared_score)$p.value}
	if(length(primary_score)>0 & length(recurrent_score)>0){
		wilcox_pr[i] <- wilcox.test(primary_score,recurrent_score)$p.value}
	if(length(recurrent_score)>0 & length(shared_score)>0){
		wilcox_rs[i] <- wilcox.test(recurrent_score,shared_score)$p.value}
	
	subtype[i] <- sub_neo[1,"idh_codel_subtype"]
	surg_int[i] <- sub_neo[1,"surgical_interval_mo"]
	norm_contam_a[i] <- sub_neo[1,"normal_contamination_a"]
	norm_contam_b[i] <- sub_neo[1,"normal_contamination_b"]
	hypermutator_recur[i] <- as.numeric(sub_neo[1,"hypermutator_status"])
	tmz[i] <- sub_neo[1,"received_tmz"]
	rt[i] <- sub_neo[1,"received_rt"]
	immunotherapy[i] <- sub_neo[1,"treatment_chemotherapy_other"]
	tumor_barcode_a[i] <- sapply(strsplit(sub_neo[1,"tumor_barcode_a"],"-"),function(x)paste(x[1:4],collapse="-"))
	tumor_barcode_b[i] <- sapply(strsplit(sub_neo[1,"tumor_barcode_b"],"-"),function(x)paste(x[1:4],collapse="-"))
}
prop_ps <- n_primary/(n_shared+n_primary)
res <- data.frame(pairs,tumor_barcode_a,tumor_barcode_b,n_primary,n_recurrent,n_shared,med_primary,med_recurrent,med_shared,prop_ps,
				  subtype,surg_int,tmz,rt,hypermutator_recur,immunotherapy,norm_contam_a,norm_contam_b)
res[which(is.na(res[,"surg_int"])),"surg_int"] = 0
res[,"immunotherapy"] <- as.character(res[,"immunotherapy"])
res[-grep("Pembrolizumab|dcVax",res[,"immunotherapy"]),"immunotherapy"] <- "None"
res[grep("Pembrolizumab",res[,"immunotherapy"]),"immunotherapy"] <- "Pembrolizumab"
res[,"tmz"] <- as.factor(res[,"tmz"])
res[,"rt"] <- as.factor(res[,"rt"])
res[,"hypermutator_recur"] <- as.factor(res[,"hypermutator_recur"])
res <- res[order(res[,"prop_ps"]),]

plot_res <- res
plot_res[,"ps_dif"] <- plot_res[,"med_primary"] - plot_res[,"med_shared"]
plot_res[,"pr_dif"] <- plot_res[,"med_primary"] - plot_res[,"med_recurrent"]
plot_res <- plot_res[order(plot_res[,"ps_dif"]),]

colnames(plot_res) <- c("tumor_pair_barcode","tumor_barcode_a","tumor_barcode_b",
								"n_primary","n_recurrent","n_shared","median_primary",
								"median_recurrent","median_shared","primary_shared_proportion",
								"subtype","surgical_int","tmz","rt","hypermutator","immunotherapy",
								"norm_contam_a","norm_contam_b","primary_shared","primary_recurrent")

dbWriteTable(con, Id(schema="analysis",table="immunoediting_scores"), plot_res, overwrite=TRUE, row.names=FALSE)

