#This code makes faceted barplots that compare the relative contribution of different signatures to private
#and shared mutations in samples stratified by IDH/codel subtype and hypermutator status
#Author: Fred Varn
#--------------------------------------------------

library(DBI)
library(odbc)
library(ggplot2)
library(tidyverse)

rm(list=ls())

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

#Download
#--------------------------------------------------
#Get mutational signature info
q = "SELECT mf.tumor_pair_barcode, tp.tumor_barcode_a, tp.tumor_barcode_b, mf.fraction, mf.signature, mf.abs_score, mf.rel_score, mf.mut_n, tp.idh_codel_subtype,cc.received_alk,cc.received_rt,cc.hypermutator_status
FROM analysis.mut_sig_fraction mf 
INNER JOIN analysis.silver_set ss ON ss.tumor_pair_barcode = mf.tumor_pair_barcode 
LEFT JOIN analysis.tumor_mut_comparison_anno tp ON mf.tumor_pair_barcode = tp.tumor_pair_barcode
LEFT JOIN analysis.tumor_clinical_comparison cc ON mf.tumor_pair_barcode = cc.tumor_pair_barcode"

mutsig_table <- dbGetQuery(con,q)

#Remove WGS
mutsig_table <- mutsig_table[grep("-WXS",mutsig_table[,1]),]
mutsig_table[,"fraction"] <- as.factor(mutsig_table[,"fraction"])

#Plot non-hypermutators
#--------------------------------------------------

nhm_mutsig_table <- mutsig_table[-which(mutsig_table[,"hypermutator_status"]==1),]
subtypes <- unique(nhm_mutsig_table[,"idh_codel_subtype"])

length(unique(nhm_mutsig_table[,"tumor_pair_barcode"]))

#-------------------------

#Kruskal wallis table
nhm_kruskal_table <- nhm_eff_table <- matrix(0,ncol=length(subtypes),nrow=30)
colnames(nhm_kruskal_table) = colnames(nhm_eff_table) = subtypes
rownames(nhm_kruskal_table) <- rownames(nhm_eff_table) <- as.character(1:30)
n = rep(0,3)
for(i in 1:length(subtypes))
{
	for(j in 1:30)
	{
		tmp_subtype <- subtypes[i]
		tmp_sig <- as.character(j)
		nhm_mutsig_table_sub <- nhm_mutsig_table[which(nhm_mutsig_table[,"idh_codel_subtype"]==tmp_subtype & nhm_mutsig_table[,"signature"]==tmp_sig),]
		nhm_kruskal_table[j,i] <- kruskal.test(rel_score~fraction,data=nhm_mutsig_table_sub)$p.value
		nhm_eff_table[j,i] <- mean(nhm_mutsig_table_sub[,"rel_score"])
		n[i] = nrow(nhm_mutsig_table_sub)
	}
}
#Adjusted p-value:
nhm_kruskal_table <- apply(nhm_kruskal_table,2,function(x)p.adjust(x,"BH"))

#Effects:
nhm_means <- apply(nhm_eff_table,1,mean)
nhm_means[order(nhm_means)]

#Plot all signatures where there was a significant difference in at least one subtype:
#keep <- apply(nhm_kruskal_table,1,function(x)sum(x<0.05))
keep <- apply(nhm_eff_table>0.04,1,sum)
keep <- names(keep[which(keep>0)])

nhm_mutsig_table <- nhm_mutsig_table[which(nhm_mutsig_table[,"signature"] %in% keep),]
nhm_mutsig_table[,"signature"] <- factor(nhm_mutsig_table[,"signature"],levels=as.character(keep))
nhm_mutsig_table[,"fraction"] <- factor(nhm_mutsig_table[,"fraction"],levels=c("S","P","R"))
nhm_mutsig_table[which(nhm_mutsig_table[,"idh_codel_subtype"]=="IDHwt_noncodel"),"idh_codel_subtype"] <- "IDHwt"
nhm_mutsig_table[which(nhm_mutsig_table[,"idh_codel_subtype"]=="IDHmut_codel"),"idh_codel_subtype"] <- "IDHmut codel"
nhm_mutsig_table[which(nhm_mutsig_table[,"idh_codel_subtype"]=="IDHmut_noncodel"),"idh_codel_subtype"] <- "IDHmut noncodel"
nhm_mutsig_table[,"idh_codel_subtype"] <- factor(nhm_mutsig_table[,"idh_codel_subtype"])

#Stats on signature 1 in manuscript
g1 <- mean(nhm_mutsig_table[which(nhm_mutsig_table[,"fraction"]=="S" & nhm_mutsig_table[,"signature"]=="1"),"rel_score"])		#0.59
g2 <- mean(nhm_mutsig_table[which(nhm_mutsig_table[,"fraction"]=="P" & nhm_mutsig_table[,"signature"]=="1"),"rel_score"])		#0.34
g3 <- mean(nhm_mutsig_table[which(nhm_mutsig_table[,"fraction"]=="R" & nhm_mutsig_table[,"signature"]=="1"),"rel_score"])		#0.23


source("/projects/varnf/SofWar/R/geom_uperrorbar.r")
#Plot contribution barplot
pdf("/projects/varnf/GLASS/Figures/signatures/finals/barplots_nhm_resub.pdf",width=7,height=1.8)
ggplot(nhm_mutsig_table, aes(y=rel_score, x=signature, fill=fraction)) +
stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="uperrorbar", position=position_dodge(width=0.90),color="gray", width=0.5) +
geom_bar(stat="summary", fun.y = "mean", position="dodge", colour="black") +
facet_grid(cols=vars(idh_codel_subtype)) +
scale_fill_manual(values=c("#CA932F","#CA2F66","#2FB3CA")) +
theme_bw()+
theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), 
	panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),
	axis.title.y = element_text(size=7), axis.text.y = element_text(size = 7), 
	axis.title.x = element_blank(), axis.text.x = element_text(size = 7, 
	vjust = 0.4), strip.text.x = element_text(size = 7), 
    strip.text.y = element_text(size = 7), 
    panel.spacing.x = unit(0, "lines"),
    legend.position="none") +
coord_cartesian(ylim=c(0,1.0))
dev.off()


#Plot hypermutators
#--------------------------------------------------

hm_mutsig_table <- mutsig_table[which(mutsig_table[,"hypermutator_status"]==1 & mutsig_table[,"received_alk"] ==1),]
length(unique(hm_mutsig_table[,"tumor_pair_barcode"]))	#n = 25


#Kruskal wallis table
hm_kruskal_table <- matrix(0,ncol=length(subtypes),nrow=30)
colnames(hm_kruskal_table) = subtypes
rownames(hm_kruskal_table) <- as.character(1:30)
n = rep(0,3)
for(i in 1:length(subtypes))
{
	for(j in 1:30)
	{
		tmp_subtype <- subtypes[i]
		tmp_sig <- as.character(j)
		hm_mutsig_table_sub <- hm_mutsig_table[which(hm_mutsig_table[,"idh_codel_subtype"]==tmp_subtype & hm_mutsig_table[,"signature"]==tmp_sig),]
		hm_kruskal_table[j,i] <- kruskal.test(rel_score~fraction,data=hm_mutsig_table_sub)$p.value
		n[i] = nrow(hm_mutsig_table_sub)
	}
}
#Adjusted p-value:
hm_kruskal_table <- apply(hm_kruskal_table,2,function(x)p.adjust(x,"BH"))

#Plot all signatures where there was a significant difference in at least one subtype:
keep <- apply(hm_kruskal_table,1,function(x)sum(x<0.05))
keep <- names(keep[which(keep>0)])

hm_mutsig_table <- hm_mutsig_table[which(hm_mutsig_table[,"signature"] %in% keep),]
hm_mutsig_table[,"signature"] <- factor(hm_mutsig_table[,"signature"],levels=as.character(keep))
hm_mutsig_table[,"fraction"] <- factor(hm_mutsig_table[,"fraction"],levels=c("S","P","R"))
hm_mutsig_table[which(hm_mutsig_table[,"idh_codel_subtype"]=="IDHwt_noncodel"),"idh_codel_subtype"] <- "IDHwt"
hm_mutsig_table[which(hm_mutsig_table[,"idh_codel_subtype"]=="IDHmut_codel"),"idh_codel_subtype"] <- "IDHmut codel"
hm_mutsig_table[which(hm_mutsig_table[,"idh_codel_subtype"]=="IDHmut_noncodel"),"idh_codel_subtype"] <- "IDHmut noncodel"
hm_mutsig_table[,"idh_codel_subtype"] <- factor(hm_mutsig_table[,"idh_codel_subtype"])

#Plot 1: Stacked barplots for all samples
#Plot contribution barplot
pdf("/projects/varnf/GLASS/Figures/signatures/finals/barplots_hm_resub.pdf",width=7,height=1.8)
ggplot(hm_mutsig_table, aes(y=rel_score, x=signature, fill=fraction)) +
stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="uperrorbar", position=position_dodge(width=0.90),color="gray", width=0.5) +
geom_bar(stat="summary", fun.y = "mean", position="dodge", colour="black") +
facet_grid(cols=vars(idh_codel_subtype)) +
scale_fill_manual(values=c("#CA932F","#CA2F66","#2FB3CA")) +
theme_bw()+
theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), 
	panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),
	axis.title.y = element_text(size=7), axis.text.y = element_text(size = 7), 
	axis.title.x = element_blank(), axis.text.x = element_text(size = 7, 
	vjust = 0.4), strip.text.x = element_text(size = 7), 
    strip.text.y = element_text(size = 7), 
    panel.spacing.x = unit(0, "lines"),
    legend.position="none") +
coord_cartesian(ylim=c(0,1.0))
dev.off()