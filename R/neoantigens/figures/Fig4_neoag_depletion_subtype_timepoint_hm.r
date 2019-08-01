#Code for Figure 4B,4C: Compares neoantigen depletion rates between subtypes and timepoints
#Code for Extended Data Figure 12B: Hypermutator analyses
#Also performs other exploratory analyses and statistical tests between groups
#Code making Figures is labelled clearly
#Query at the top joins analysis.neoantigen_depletion materialized view with other tables to get timepoint, hypermutator, and subtype information
#All samples in these analyses were required to have at least 3 missense mutations in their initial and recurrent tumors (avoid overfitting ratio)
#-----------------------------------------------------

library(odbc)
library(DBI)
library(ggplot2)
library(ggpubr)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")  

q <- "SELECT gs.* , nd1.rneo AS nd_a, nd2.rneo AS nd_b, clin.idh_codel_subtype AS subtype,
CASE WHEN mf1.coverage_adj_mut_freq >= 10 THEN 1 WHEN mf1.coverage_adj_mut_freq < 10 THEN 0 END AS hm_a,
CASE WHEN mf2.coverage_adj_mut_freq >= 10 THEN 1 WHEN mf2.coverage_adj_mut_freq < 10 THEN 0 END AS hm_b
FROM analysis.gold_set gs
LEFT JOIN analysis.neoantigen_depletion nd1 ON gs.tumor_barcode_a = nd1.aliquot_barcode
LEFT JOIN analysis.neoantigen_depletion nd2 ON gs.tumor_barcode_b = nd2.aliquot_barcode
LEFT JOIN analysis.mut_freq mf1 ON mf1.aliquot_barcode = gs.tumor_barcode_a
LEFT JOIN analysis.mut_freq mf2 ON mf2.aliquot_barcode = gs.tumor_barcode_b
LEFT JOIN clinical.subtypes clin ON clin.case_barcode = gs.case_barcode
WHERE nd1.rneo IS NOT NULL AND nd2.rneo IS NOT NULL AND (nd1.nobs >= 3 AND nd2.nobs >= 3)
ORDER BY nd1.rneo"

res <- dbGetQuery(con, q)

tmp <- res[,"subtype"]
tmp <- tmp[order(tmp)]
rle(tmp)
#Scatterplot of primary vs recurrent correlations with depletion
cor.test(res[,"nd_a"],res[,"nd_b"])
cor.test(res[,"nd_a"],res[,"nd_b"])$p.value

summary(aov(c(res[,"nd_a"],res[,"nd_b"])~rep(res[,"subtype"],2))) #P = 0.09

#******************************
#Figure 4C
#******************************
gtsize = 7/(14/5)
pdf("/projects/varnf/GLASS/Figures/resubmission/Figure4C.pdf",width=2,height=2)
ggplot(res, aes(x = nd_a, y = nd_b,colour=subtype)) + 
	geom_point(size=1.5,alpha=0.4) +
	geom_abline(intercept=0, slope=1, color = "gray50", size=0.5) +
 	geom_text(aes(x=1.75, y=0, label = "R = 0.73"), size=gtsize, group = NA, colour = "black") +
	labs(x="Observed/expected neoantigens\n(initial)", y="Observed/expected neoantigens\n(recurrent)") +
	theme_classic() +
	theme(axis.title=element_text(size=7),
	axis.text=element_text(colour="black",size=7),
	legend.position="none") +
	coord_cartesian(ylim=c(0,2),xlim=c(0,2))
dev.off()

gtsize = 7/(14/5)
pdf("/projects/varnf/GLASS/Figures/resubmission/Figure4C_bigger.pdf",width=2,height=2)
ggplot(res, aes(x = nd_a, y = nd_b,colour=subtype)) + 
	geom_point(size=1.5,alpha=0.4) +
	geom_abline(intercept=0, slope=1, color = "gray50", size=0.5) +
 	geom_text(aes(x=1.75, y=0, label = "R = 0.73"), size=gtsize, group = NA, colour = "black") +
	labs(x="Observed/expected neoantigens (initial)", y="Observed/expected neoantigens (recurrent)") +
	theme_classic() +
	theme(axis.title=element_text(size=7),
	axis.text=element_text(colour="black",size=7),
	legend.position="none") +
	coord_cartesian(ylim=c(0,2),xlim=c(0,2))
dev.off()

#Examine differences in the residuals
difs <- res[,"nd_b"] - res[,"nd_a"]
wilcox.test(difs[which(res[,"subtype"]=="IDHwt")],difs[which(res[,"subtype"]=="IDHmut-noncodel")])
wilcox.test(difs[which(res[,"subtype"]=="IDHwt")],difs[which(res[,"subtype"]=="IDHmut-codel")])
wilcox.test(difs[which(res[,"subtype"]=="IDHmut-noncodel")],difs[which(res[,"subtype"]=="IDHmut-codel")])
wilcox.test(difs[which(res[,"subtype"]=="IDHwt")],difs[which(res[,"subtype"]=="IDHmut-codel")])
wilcox.test(difs[which(res[,"subtype"]!="IDHmut-codel")],difs[which(res[,"subtype"]=="IDHmut-codel")])
wilcox.test(difs[which(res[,"subtype"]!="IDHwt")],difs[which(res[,"subtype"]!="IDHwt")])


#Compare primary vs recurrent irrespective of subtype:
t.test(res[,"nd_a"],res[,"nd_b"],paired=TRUE)		#Initial vs recurrent:	0.77
t.test(res[,"nd_a"], mu=1)							#Initial: 				0.47
t.test(res[,"nd_b"], mu=1)							#Recurrent:				0.59

#Compare subtypes irrespective of timepoint:
idhwt <- res[which(res[,"subtype"]=="IDHwt"),]
idhmutnoncodel <- res[which(res[,"subtype"]=="IDHmut-noncodel"),]
idhmutcodel <- res[which(res[,"subtype"]=="IDHmut-codel"),]
idhmut <- res[which(res[,"subtype"]!="IDHwt"),]

t.test(c(idhwt[,"nd_a"],idhwt[,"nd_b"]),c(idhmutnoncodel[,"nd_a"],idhmutnoncodel[,"nd_b"]))						#0.04
t.test(c(idhwt[,"nd_a"],idhwt[,"nd_b"]),c(idhmutcodel[,"nd_a"],idhmutcodel[,"nd_b"]))							#0.39
t.test(c(idhmutnoncodel[,"nd_a"],idhmutnoncodel[,"nd_b"]),c(idhmutcodel[,"nd_a"],idhmutcodel[,"nd_b"]))			#0.61
t.test(c(idhwt[,"nd_a"],idhwt[,"nd_b"]),c(idhmut[,"nd_a"],idhmut[,"nd_b"]))										#0.04

#Compare primary depletion in each subtype:
g1 <- res[which(res[,"subtype"]=="IDHwt"),"nd_a"]
g2 <- res[which(res[,"subtype"]=="IDHmut-noncodel"),"nd_a"]
g3 <- res[which(res[,"subtype"]=="IDHmut-codel"),"nd_a"]


t.test(g1, mu=1)					#IDHwt:									0.66
t.test(g2, mu=1)					#IDHmut-noncodel:						0.22
t.test(g3, mu=1)					#IDHmut-codel:							0.32
t.test(c(g2,g3), mu=1)				#IDHmut:								0.12

t.test(g1,g2)						#IDHwt vs IDHmut-noncodel:				0.20
t.test(g1,g3)						#IDHwt vs IDHmut-codel:					0.28
t.test(g2,g3)						#IDHmut-noncodel vs IDHmut-codel:		0.89
t.test(g1,c(g2,g3))				#IDHwt vs IDHmut:							0.13


#Compare recurrent depletion in each subtype:
g4 <- res[which(res[,"subtype"]=="IDHwt"),"nd_b"]
g5 <- res[which(res[,"subtype"]=="IDHmut-noncodel"),"nd_b"]
g6 <- res[which(res[,"subtype"]=="IDHmut-codel"),"nd_b"]

t.test(g4, mu=1)					#IDHwt:									0.62
t.test(g5, mu=1)					#IDHmut-noncodel:						0.09
t.test(g6, mu=1)					#IDHmut-codel:							0.98
t.test(c(g5,g6), mu=1)				#IDHmut:								0.17

t.test(g4,g5)						#IDHwt vs IDHmut-noncodel:				0.09
t.test(g4,g6)						#IDHwt vs IDHmut-codel:					0.88
t.test(g5,g6)						#IDHmut-noncodel vs IDHmut-codel:		0.41
t.test(g4,c(g5,g6))				#IDHwt vs IDHmut:							0.16


#Compare primary vs recurrent depletion in each subtype
t.test(g1,g4,paired=TRUE)				#IDHwt:								0.92
t.test(g2,g5,paired=TRUE)				#IDHmut-noncodel:					0.75
t.test(g3,g6,paired=TRUE)				#IDHmut-codel:						0.09
t.test(c(g2,g3),c(g5,g6),paired=TRUE)	#IDHmut:							0.77

#Compare hypermutators at recurrence vs non-hypermutators at recurrence
t.test(res[which(res[,"hm_b"]==1),"nd_b"],res[which(res[,"hm_b"]==0),"nd_b"])	#0.54


#Subtype and timepoint figures
#-----------------------------------------------------

pair <- c(rep(res[,"tumor_pair_barcode"],2))
samp <- c(res[,"tumor_barcode_a"],res[,"tumor_barcode_b"])
rneo <- c(res[,"nd_a"],res[,"nd_b"])
subtype <- rep(res[,"subtype"],2)
timepoint <- c(rep("Initial",nrow(res)),rep("Recurrent",nrow(res)))
plot_res <- data.frame(pair,samp,rneo,subtype,timepoint)

#Old figure don't use
pdf("/projects/varnf/GLASS/Figures/resubmission/depletion_subtype_timepoint_box.pdf",width=2,height=2)
ggplot(plot_res, aes(y = rneo, x = timepoint,fill=subtype)) + 
	geom_hline(yintercept=1, color = "gray50", size=1) +
	geom_boxplot(lwd=0.5,outlier.size=0.1,fatten=2) +
	labs(y="observed/expected neoantigens") +
	theme_classic() +
	theme(strip.text.x = element_text(size=7),strip.background = element_blank(),
	axis.title.y=element_text(size=7),axis.title.x=element_blank(),
	axis.text=element_text(colour="black",size=7), axis.text.x=element_text(hjust=0.5,size=7),
	legend.position="none")
dev.off()

#******************************
#Figure 4B
#******************************
t.test(plot_res[which(plot_res[,"subtype"]=="IDHwt"),"rneo"],plot_res[which(plot_res[,"subtype"]=="IDHmut-noncodel"),"rneo"]) #0.04
pdf("/projects/varnf/GLASS/Figures/resubmission/Figure4B.pdf",width=3.75,height=2)
ggplot(plot_res, aes(y = rneo, x = subtype,colour=subtype)) + 
	geom_boxplot(lwd=0.5,outlier.size=0.1,fatten=2,colour="black",fill="white") +
	geom_jitter(size=1,width=0.05) +	
	labs(y="Observed/expected neoantigens") +
 	geom_text(aes(x=2.5, y=1.8, label = "P = 0.04"), size=gtsize, group = NA, color = "black") +
	theme_classic() + 
	theme(strip.text.x = element_text(size=7),strip.background = element_blank(),
	axis.title.x=element_text(size=7),axis.title.y=element_blank(),
	axis.text=element_text(colour="black",size=7), axis.text.x=element_text(hjust=0.5,size=7),
	legend.position="none") + 
	coord_flip(ylim=c(0,2.5))
dev.off()

#Old figure, replaced with correlation plot
pval <- rep("",nrow(plot_res))
pval[which(plot_res[,"subtype"]=="IDHwt")] <- "P = 0.70"
pval[which(plot_res[,"subtype"]=="IDHmut-noncodel")] <- "P = 0.77"
pval[which(plot_res[,"subtype"]=="IDHmut-codel")] <- "P = 0.14"
p_text <- cbind(plot_res,pval)
pdf("/projects/varnf/GLASS/Figures/resubmission/depletion_subtype_timepoint_ladder.pdf",width=4,height=2)
ggplot(plot_res, aes(x = timepoint, y = rneo, group = pair,colour=subtype)) + 
	geom_line(size=0.45) +
	geom_point(size=1,colour="black") +
	labs(y="Observed/expected neoantigens") +
 	geom_text(data=p_text, aes(x=1.5, y=2, label = pval), size=gtsize, color = "black") +
	facet_grid(.~subtype) +
	theme_classic() +
	theme(strip.text.x = element_text(size=7),strip.background = element_blank(),
	axis.title.y=element_text(size=7),axis.title.x=element_blank(),
	axis.text=element_text(colour="black",size=7), axis.text.x=element_text(hjust=0.5,size=7),
	legend.position="none")
dev.off()

#**********************************************************
#Extended Data Figure 12B: Hypermutator figures (two panels)
#**********************************************************

samp <- res[,"tumor_barcode_b"]
rneo <- res[,"nd_b"]
subtype <- res[,"subtype"]
hm <- as.factor(res[,"hm_b"])
plot_res_hm <- data.frame(samp,rneo,subtype,hm)
t.test(plot_res_hm[which(plot_res_hm[,"hm"]==1),"rneo"],plot_res_hm[which(plot_res_hm[,"hm"]==0),"rneo"]) #0.54


pdf("/projects/varnf/GLASS/Figures/resubmission/EDF12_depletion_hypermutator_recur_box.pdf",width=2,height=2)
ggplot(plot_res_hm, aes(y = rneo, x = hm, colour=hm)) + 
	geom_boxplot(lwd=0.5,outlier.size=0.1,fatten=2,colour="black",fill="white") +
	geom_jitter(size=1,width=0.05) +	
	labs(y="Observed/expected neoantigens") +
 	geom_text(aes(x=1.5, y=2, label = "P = 0.54"), size=gtsize, group = NA, color = "black") +
	scale_colour_manual(values=c("royalblue4","tomato3")) +
	theme_classic() +
	theme(strip.text.x = element_text(size=7),strip.background = element_blank(),
	axis.title.y=element_text(size=7),axis.title.x=element_blank(),
	axis.text=element_text(colour="black",size=7), axis.text.x=element_text(hjust=0.5,size=7),
	legend.position="none")
dev.off()


tmp <- res[which(res[,"hm_b"]==1),]
pair <- c(rep(tmp[,"tumor_pair_barcode"],2))
samp <- c(tmp[,"tumor_barcode_a"],tmp[,"tumor_barcode_b"])
rneo <- c(tmp[,"nd_a"],tmp[,"nd_b"])
subtype <- rep(tmp[,"subtype"],2)
hm <-  c(tmp[,"hm_a"],tmp[,"hm_b"])
timepoint <- c(rep("Initial",nrow(tmp)),rep("Recurrent",nrow(tmp)))
plot_res_hmtp <- data.frame(pair,samp,rneo,hm,subtype,timepoint)

t.test(plot_res_hmtp[which(plot_res_hmtp[,"timepoint"]=="Initial"),"rneo"],plot_res_hmtp[which(plot_res_hmtp[,"timepoint"]=="Recurrent"),"rneo"],paired=TRUE) #P = 0.32
#P = 0.24

pdf("/projects/varnf/GLASS/Figures/resubmission/EDF12_depletion_hypermutator_timepoint_box.pdf",width=2,height=2)
ggplot(plot_res_hmtp, aes(y = rneo, x = timepoint, group=pair, colour=subtype)) + 
	geom_line(size=0.45) +
	geom_point(size=1,colour="black") +
	labs(y="Observed/expected neoantigens") +
 	geom_text(aes(x=1.5, y=2, label = "P = 0.32"), size=gtsize, group = NA, color = "black") +
	theme_classic() +
	theme(strip.text.x = element_text(size=7),strip.background = element_blank(),
	axis.title.y=element_text(size=7),axis.title.x=element_blank(),
	axis.text=element_text(colour="black",size=7), axis.text.x=element_text(hjust=0.5,size=7),
	legend.position="none")
dev.off()