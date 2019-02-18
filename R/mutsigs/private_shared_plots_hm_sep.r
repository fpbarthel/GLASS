#This code makes faceted barplots that compare the relative contribution of different signatures to private
#and shared mutations in samples stratified by IDH/codel subtype and hypermutator status
#Author: Fred Varn
#--------------------------------------------------

library(DBI)
library(odbc)
library(ggplot2)

rm(list=ls())

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

#Download
#--------------------------------------------------
#Get mutational signature info
q = "SELECT mf.tumor_pair_barcode, mf.fraction, mf.signature, mf.signature_score, mf.mut_count \
FROM analysis.mutsig_fraction mf \
INNER JOIN analysis.silver_set ss ON ss.tumor_pair_barcode = mf.tumor_pair_barcode \
WHERE ss.priority = 1"

mutsig <- dbGetQuery(con,q)

#Get pair information 

q = "SELECT ctp.*, tp.idh_codel_subtype \
FROM analysis.clinical_by_tumor_pair ctp \
INNER JOIN analysis.tumor_mut_comparison_anno tp ON tp.tumor_pair_barcode = ctp.tumor_pair_barcode"

pair_info <- dbGetQuery(con,q)

#Clean table
#--------------------------------------------------

mutsig_table <- merge(mutsig,pair_info,by="tumor_pair_barcode")
mutsig_table <- mutsig_table[order(mutsig_table[,"tumor_pair_barcode"],mutsig_table[,"signature"],mutsig_table[,"fraction"]),]

#Remove WGS
mutsig_table <- mutsig_table[grep("-WXS",mutsig_table[,1]),]

#Plot non-hypermutators
#--------------------------------------------------

nhm_mutsig_table <- mutsig_table[-which(mutsig_table[,"hypermutator_status"]==1),]
subtypes <- unique(nhm_mutsig_table[,"idh_codel_subtype"])

#Remove signatures with tiny contributions
#sig_avg <- aggregate(nhm_mutsig_table[,"signature_score"],by=list(nhm_mutsig_table[,"signature"]),mean)
#sig_pull <- sig_avg[which(sig_avg[,2]>0.05),1]
#nhm_mutsig_table_sub <- nhm_mutsig_table[which(nhm_mutsig_table[,"signature"]%in%sig_pull),]

#-------------------------
#Compare primary to recurrent
for(i in 1:length(subtypes))
{
	for(j in 1:30)
	{
		tmp_subtype <- subtypes[i]
		tmp_sig <- as.character(j)
		nhm_mutsig_table_sub <- nhm_mutsig_table[which(nhm_mutsig_table[,"idh_codel_subtype"]==tmp_subtype & nhm_mutsig_table[,"signature"]==tmp_sig),]
		g1 = nhm_mutsig_table_sub[which(nhm_mutsig_table_sub[,"fraction"]=="P"),"signature_score"]
		g2 = nhm_mutsig_table_sub[which(nhm_mutsig_table_sub[,"fraction"]=="R"),"signature_score"]
		nhm_wilcox_table[j,i] <- wilcox.test(g1,g2)$p.value
		nhm_eff_table[j,i] <- mean(g1) - mean(g2)
	}
}
nhm_wilcox_table <- apply(nhm_wilcox_table,2,function(x)p.adjust(x,"BH"))
#-------------------------

#Compare shared to private (use this for plotting)
nhm_wilcox_table <- nhm_eff_table <- matrix(0,ncol=length(subtypes),nrow=30)
colnames(nhm_wilcox_table) = colnames(nhm_eff_table) = subtypes
rownames(nhm_wilcox_table) = rownames(nhm_eff_table) =1:30
n = rep(0,3)
for(i in 1:length(subtypes))
{
	for(j in 1:30)
	{
		tmp_subtype <- subtypes[i]
		tmp_sig <- as.character(j)
		nhm_mutsig_table_sub <- nhm_mutsig_table[which(nhm_mutsig_table[,"idh_codel_subtype"]==tmp_subtype & nhm_mutsig_table[,"signature"]==tmp_sig),]
		g1 = nhm_mutsig_table_sub[which(nhm_mutsig_table_sub[,"fraction"]=="S"),"signature_score"]
		g2 = nhm_mutsig_table_sub[which(nhm_mutsig_table_sub[,"fraction"]!="S"),"signature_score"]
		nhm_wilcox_table[j,i] <- wilcox.test(g1,g2)$p.value
		nhm_eff_table[j,i] <- mean(g1) - mean(g2)
		n[i] = length(g1)
	}
}
#Adjusted p-value:
nhm_wilcox_table <- apply(nhm_wilcox_table,2,function(x)p.adjust(x,"BH"))

#Plot all signatures where there was a significant difference in at least one subtype:
keep <- apply(nhm_wilcox_table,1,function(x)sum(x<0.1))
keep <- names(keep[which(keep>0)])

nhm_mutsig_table <- nhm_mutsig_table[which(nhm_mutsig_table[,"signature"] %in% keep),]
nhm_mutsig_table[,"signature"] <- factor(nhm_mutsig_table[,"signature"],levels=as.character(c(1,3,12,15,18,21,23,24,29)))

#Plot 1: Stacked barplots for all samples
#Plot contribution barplot
pdf("/projects/varnf/GLASS/Figures/signatures/nonhypermut_private_shared_wxs_barplots.pdf",width=7,height=1.8)
ggplot(nhm_mutsig_table, aes(y=signature_score, x=signature, fill=fraction)) +
geom_bar(stat="summary", fun.y = "mean", position="dodge", colour="black") +
facet_grid(cols=vars(idh_codel_subtype)) +
theme_bw()+
theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), 
	panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),
	axis.title.y = element_text(size=7), axis.text.y = element_text(size = 7), 
	axis.title.x = element_blank(), axis.text.x = element_text(size = 7, 
	vjust = 0.4), strip.text.x = element_text(size = 7), 
    strip.text.y = element_text(size = 7), 
    legend.title = element_blank(),
    legend.text = element_text(size=7),
    legend.key.size = unit(0.5,unit="line"),
    panel.spacing.x = unit(0, "lines"))
dev.off()


#Plot hypermutators
#--------------------------------------------------

hm_mutsig_table <- mutsig_table[which(mutsig_table[,"hypermutator_status"]==1 & mutsig_table[,"received_tmz"] ==1),]

#Average signature 11 in recurrence
sig11_recur <- hm_mutsig_table[grep("11",hm_mutsig_table[,"signature"]),]
sig11_recur <- sig11_recur[grep("R",sig11_recur[,"fraction"]),]

subtypes <- unique(hm_mutsig_table[,"idh_codel_subtype"])

#Remove signatures with tiny contributions
#sig_avg <- aggregate(hm_mutsig_table[,"signature_score"],by=list(hm_mutsig_table[,"signature"]),mean)
#sig_pull <- sig_avg[which(sig_avg[,2]>0.05),1]
#hm_mutsig_table_sub <- hm_mutsig_table[which(hm_mutsig_table[,"signature"]%in%sig_pull),]

hm_wilcox_table <- hm_eff_table <- matrix(0,ncol=length(subtypes),nrow=30)
colnames(hm_wilcox_table) = colnames(hm_eff_table) = subtypes
rownames(hm_wilcox_table) = rownames(hm_eff_table) = 1:30
for(i in 1:length(subtypes))
{
	for(j in 1:30)
	{
		tmp_subtype <- subtypes[i]
		tmp_sig <- j
		hm_mutsig_table_sub <- hm_mutsig_table[which(hm_mutsig_table[,"idh_codel_subtype"]==tmp_subtype & hm_mutsig_table[,"signature"]==tmp_sig),]
		g1 = hm_mutsig_table_sub[which(hm_mutsig_table_sub[,"fraction"]=="S"),"signature_score"]
		g2 = hm_mutsig_table_sub[which(hm_mutsig_table_sub[,"fraction"]!="S"),"signature_score"]
		hm_wilcox_table[j,i] <- wilcox.test(g1,g2)$p.value
	}
}
hm_wilcox_table <- apply(hm_wilcox_table,2,function(x)p.adjust(x,"BH"))

#Primary vs recurrent (plotting)
hm_wilcox_table <- hm_eff_table <- matrix(0,ncol=length(subtypes),nrow=30)
colnames(hm_wilcox_table) = colnames(hm_eff_table) = subtypes
rownames(hm_wilcox_table) = rownames(hm_eff_table) = 1:30
n = rep(0,3)
for(i in 1:length(subtypes))
{
	for(j in 1:30)
	{
		tmp_subtype <- subtypes[i]
		tmp_sig <- j
		hm_mutsig_table_sub <- hm_mutsig_table[which(hm_mutsig_table[,"idh_codel_subtype"]==tmp_subtype & hm_mutsig_table[,"signature"]==tmp_sig),]
		g1 = hm_mutsig_table_sub[which(hm_mutsig_table_sub[,"fraction"]=="R"),"signature_score"]
		g2 = hm_mutsig_table_sub[which(hm_mutsig_table_sub[,"fraction"]=="P"),"signature_score"]
		hm_wilcox_table[j,i] <- wilcox.test(g1,g2)$p.value
		hm_eff_table[j,i] <- mean(g1) - mean(g2)
		n[i] = length(g1)
	}
}
hm_wilcox_table <- apply(hm_wilcox_table,2,function(x)p.adjust(x,"BH"))


#Plot all signatures where there was a significant difference in at least one subtype:
keep <- apply(hm_wilcox_table,1,function(x)sum(x<0.1,na.rm=T))
keep <- names(keep[which(keep>0)])

hm_mutsig_table <- hm_mutsig_table[which(hm_mutsig_table[,"signature"] %in% keep),]
hm_mutsig_table[,"signature"] <- as.character(hm_mutsig_table[,"signature"])

#Plot 1: Stacked barplots for all samples
#Plot contribution barplot
pdf("/projects/varnf/GLASS/Figures/signatures/hypermut_private_shared_wxs_barplots.pdf",width=7,height=1.8)
ggplot(hm_mutsig_table, aes(y=signature_score, x=signature, fill=fraction)) +
geom_bar(stat="summary", fun.y = "mean", position="dodge", colour="black") +
facet_grid(cols=vars(idh_codel_subtype)) +
theme_bw()+
theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), 
	panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),
	axis.title.y = element_text(size=7), axis.text.y = element_text(size = 7), 
	axis.title.x = element_blank(), axis.text.x = element_text(size = 7, 
	vjust = 0.4), strip.text.x = element_text(size = 7), 
    strip.text.y = element_text(size = 7), 
    legend.title = element_blank(),
    legend.text = element_text(size=7),
    legend.key.size = unit(0.5,unit="line"),
    panel.spacing.x = unit(0, "lines"))
dev.off()