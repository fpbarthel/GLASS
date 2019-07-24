#Code for Figure 4A: Compares neoantigens/nonsynonymous mutation rates
#Also provides correlation between mutation and neoantigen load (used in manuscript)
#Query at the top joins analysis.neoag_freq materialized view to other tables to make figure that shows subtype and timepoint-specific fractions
#-----------------------------------------------------
library(DBI)
library(odbc)
library(ggplot2)
library(reshape)

rm(list=ls())

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

q = "SELECT gs.tumor_pair_barcode, 
COALESCE(neo1.neoag_count::numeric,0) AS neoag_count_a, 
COALESCE(neo1.mt_count::numeric,0) AS mt_count_a, 
COALESCE(neo1.prop_immunogenic::numeric,0) AS prop_immunogenic_a, 
COALESCE(neo2.neoag_count::numeric,0) AS neoag_count_b, 
COALESCE(neo2.mt_count::numeric,0) AS mt_count_b, 
COALESCE(neo2.prop_immunogenic::numeric,0) AS prop_immunogenic_b, 
clin.idh_codel_subtype AS subtype, 
CASE WHEN mf.coverage_adj_mut_freq >= 10 THEN 1 WHEN mf.coverage_adj_mut_freq < 10 THEN 0 END AS hypermutator
FROM analysis.gold_set gs
LEFT JOIN analysis.neoag_freq neo1 ON neo1.aliquot_barcode = gs.tumor_barcode_a 
LEFT JOIN analysis.neoag_freq neo2 ON neo2.aliquot_barcode = gs.tumor_barcode_b
LEFT JOIN clinical.subtypes clin ON  clin.case_barcode = gs.case_barcode
INNER JOIN analysis.mut_freq mf ON mf.aliquot_barcode = gs.tumor_barcode_b
ORDER BY gs.tumor_pair_barcode"

neo <- dbGetQuery(con,q)


#get the n's for this plot
tmp <- neo[,"subtype"]
tmp <- tmp[order(tmp)]
rle(tmp)

#get average proportion of mutations giving rise to neoantigens
mean(c(neo[,"prop_immunogenic_a"],neo[,"prop_immunogenic_b"]))		#42%

#get correlation between mutation and neoantigen load
cor(c(neo[,"neoag_count_a"],neo[,"neoag_count_b"]),c(neo[,"mt_count_a"],neo[,"mt_count_b"]),method="s")	#0.89


prop <- c(neo[,"prop_immunogenic_a"],neo[,"prop_immunogenic_b"])
pair <- rep(neo[,"tumor_pair_barcode"],2)
subtype <- rep(neo[,"subtype"],2)
timepoint <- c(rep("Initial",nrow(neo)),rep("Recurrent",nrow(neo)))
res <- data.frame(pair,prop,subtype,timepoint)

subtypes <- unique(neo[,"subtype"])
prop_pri <- prop_rec <- neoag_sd_pri <- neoag_sd_rec <- rep(0,length(subtypes))
for(i in 1:length(subtypes))
{
	sub_res <- res[which(res[,"subtype"]==subtypes[i]),]
	prop_pri[i] <- mean(sub_res[which(sub_res[,"timepoint"]=="Initial"),"prop"])
	prop_rec[i] <- mean(sub_res[which(sub_res[,"timepoint"]=="Recurrent"),"prop"])
	neoag_sd_pri[i] <- sd(sub_res[which(sub_res[,"timepoint"]=="Initial"),"prop"])
	neoag_sd_rec[i] <- sd(sub_res[which(sub_res[,"timepoint"]=="Recurrent"),"prop"])

}

#Significance tests
g1 <- res[which(res[,"timepoint"]=="Initial" & res[,"subtype"] == "IDHmut-codel"),"prop"]
g2 <- res[which(res[,"timepoint"]=="Initial" & res[,"subtype"] == "IDHmut-noncodel"),"prop"]
g3 <- res[which(res[,"timepoint"]=="Initial" & res[,"subtype"] == "IDHwt"),"prop"]

g4 <- res[which(res[,"timepoint"]=="Recurrent" & res[,"subtype"] == "IDHmut-codel"),"prop"]
g5 <- res[which(res[,"timepoint"]=="Recurrent" & res[,"subtype"] == "IDHmut-noncodel"),"prop"]
g6 <- res[which(res[,"timepoint"]=="Recurrent" & res[,"subtype"] == "IDHwt"),"prop"]

wilcox.test(g1,g4)		#0.35
wilcox.test(g2,g5)		#0.82
wilcox.test(g3,g6)		#0.90

wilcox.test(g1,g2)		#0.19
wilcox.test(g1,g3)		#0.38
wilcox.test(g2,g3)		#0.61

wilcox.test(g4,g5)		#0.66
wilcox.test(g4,g6)		#0.92
wilcox.test(g5,g6)		#0.26

#Plot A
plot_res <- data.frame(c(prop_pri,prop_rec),
			c(neoag_sd_pri,neoag_sd_rec),
			c(rep("Initial",3),rep("Recurrent",3)),
			c(rep("prop_immunogenic",6)),
			rep(subtypes,2))
colnames(plot_res) <- c("fraction","sd","status","neoag","subtypes")

pdf("/projects/varnf/GLASS/Figures/resubmission/final/Figure4A.pdf",width=2,height=2)
p1 <- ggplot(plot_res,aes(y = fraction, x = status, fill = subtypes)) +
	geom_bar(position="dodge",stat="identity",color="black") +
	geom_errorbar(aes(ymin=fraction-sd,ymax=fraction+sd),size=0.4,width=0.35,position=position_dodge(.9))+
	theme_bw() +
	labs(x = "", y = "Neoantigens/nonsynonymous") +
	theme(axis.text.x=element_text(size=7),axis.text.y = element_text(size=7),
	axis.title.x = element_blank(),axis.title.y = element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	legend.position="none") +
	coord_cartesian(ylim=c(0,0.60))
p1
dev.off()


#Other analyses
#-------------------------------------

#Mut and neoag load comparisons
neo_box <- neo
neo_box[,"neoag_count"] <- as.numeric(neo_box[,"neoag_count"])
neo_box[,"mt_count"] <- as.numeric(neo_box[,"mt_count"])
neo_box <- neo_box[which(neo_box[,"hypermutator"]==0),]
pdf("/projects/varnf/GLASS/Figures/resubmission/neaog_load_boxplot.pdf",width=4,height=2)
p1 <- ggplot(neo_box,aes(y = neoag_count, x = sample_type,fill=subtype)) +
	geom_boxplot(outlier.shape=NA) +
	geom_line(aes(group=case_barcode),linetype=2,colour="gray50",size=0.1) +
	geom_point(size=0.5) +
	theme_bw() +
	facet_grid(.~subtype) +
	labs(x = "", y = "Mutation count") +
	theme(axis.text.x=element_text(size=7),axis.text.y = element_text(size=7),
	axis.title.x = element_blank(),axis.title.y = element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	legend.position="none")
p1
dev.off()

pdf("/projects/varnf/GLASS/Figures/resubmission/mut_load_boxplot.pdf",width=4,height=2)
p1 <- ggplot(neo_box,aes(y = mt_count, x = sample_type,fill=subtype)) +
	geom_boxplot(outlier.shape=NA) +
	geom_line(aes(group=case_barcode),linetype=2,colour="gray50",size=0.1) +
	geom_point(size=0.5) +
	theme_bw() +
	facet_grid(.~subtype) +
	labs(x = "", y = "Neoantigen count") +
	theme(axis.text.x=element_text(size=7),axis.text.y = element_text(size=7),
	axis.title.x = element_blank(),axis.title.y = element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	legend.position="none")
p1
dev.off()


#Hypermutator plot
g1 <- neo[which(neo[,"hypermutator"]==1),"prop_immunogenic"]
g2 <- neo[which(neo[,"hypermutator"]==0 & neo[,"sample_type"]=="R"),"prop_immunogenic"]


hm <- neo[which(neo[,"hypermutator"]==1),]
subtypes <- unique(neo[,"subtype"])
prop_immunogenic_pri <- prop_immunogenic_rec <- neoag_sd_pri <- neoag_sd_rec <- rep(0,length(subtypes))
for(i in 1:length(subtypes))
{
	sub_neo <- neo[which(neo[,"subtype"]==subtypes[i]),]
	prop_immunogenic_pri[i] <- mean(sub_neo[which(sub_neo[,"sample_type"]=="P"),"prop_immunogenic"])
	prop_immunogenic_rec[i] <- mean(sub_neo[which(sub_neo[,"sample_type"]=="R"),"prop_immunogenic"])
	neoag_sd_pri[i] <- sd(sub_neo[which(sub_neo[,"sample_type"]=="P"),"prop_immunogenic"])
	neoag_sd_rec[i] <- sd(sub_neo[which(sub_neo[,"sample_type"]=="R"),"prop_immunogenic"])

}

plot_res <- data.frame(c(prop_immunogenic_pri,prop_immunogenic_rec),
			c(neoag_sd_pri,neoag_sd_rec),
			c(rep("Primary",3),rep("Recurrent",3)),
			c(rep("prop_immunogenic",6)),
			rep(subtypes,2))
colnames(plot_res) <- c("fraction","sd","status","neoag","subtypes")
