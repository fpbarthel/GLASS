#Code for Figure 4D: Compares neoantigen depletion ratio between subclonal and clonal neoantigens
#Also performs statistical tests between groups (no differences)
#Query at the top joins analysis.neoantigen_depletion_clonality materialized view with other tables to get timepoint, hypermutator, and subtype information
#-----------------------------------------------------

library(DBI)
library(ggplot2)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")  

q <- "SELECT gs.* , nd1c.rneo AS nd_a_c, nd1s.rneo AS nd_a_s, nd2c.rneo AS nd_b_c, nd2s.rneo AS nd_b_s, clin.idh_codel_subtype AS subtype,
CASE WHEN mf1.coverage_adj_mut_freq >= 10 THEN 1 WHEN mf1.coverage_adj_mut_freq < 10 THEN 0 END AS hm_a,
CASE WHEN mf2.coverage_adj_mut_freq >= 10 THEN 1 WHEN mf2.coverage_adj_mut_freq < 10 THEN 0 END AS hm_b
FROM analysis.gold_set gs
LEFT JOIN analysis.neoantigen_depletion_clonality nd1c ON gs.tumor_barcode_a = nd1c.aliquot_barcode AND nd1c.clonality= 'C'
LEFT JOIN analysis.neoantigen_depletion_clonality nd1s ON gs.tumor_barcode_a = nd1s.aliquot_barcode AND nd1s.clonality= 'S'
LEFT JOIN analysis.neoantigen_depletion_clonality nd2c ON gs.tumor_barcode_b = nd2c.aliquot_barcode AND nd2c.clonality= 'C'
LEFT JOIN analysis.neoantigen_depletion_clonality nd2s ON gs.tumor_barcode_b = nd2s.aliquot_barcode AND nd2s.clonality= 'S'
LEFT JOIN analysis.mut_freq mf1 ON mf1.aliquot_barcode = gs.tumor_barcode_a
LEFT JOIN analysis.mut_freq mf2 ON mf2.aliquot_barcode = gs.tumor_barcode_b
LEFT JOIN clinical.subtypes clin ON clin.case_barcode = gs.case_barcode
WHERE nd1c.rneo IS NOT NULL AND nd1s.rneo IS NOT NULL AND nd2c.rneo IS NOT NULL AND nd2s.rneo IS NOT NULL AND
(nd1c.nobs >= 3 AND nd1s.nobs >= 3 AND nd2c.nobs >= 3 AND nd2s.nobs >= 3)
ORDER BY gs.tumor_pair_barcode"

res <- dbGetQuery(con, q)

#Get sample sizes
tmp <- res[,"subtype"]
tmp <- tmp[order(tmp)]
rle(tmp)

#Clonal vs subclonal in each tumor
wilcox.test(res[,"nd_a_c"],res[,"nd_a_s"],paired=TRUE)
wilcox.test(res[,"nd_b_c"],res[,"nd_b_s"],paired=TRUE)

#Clonality status by timepoint
wilcox.test(res[,"nd_a_c"],res[,"nd_b_c"],paired=TRUE)
wilcox.test(res[,"nd_a_s"],res[,"nd_b_s"],paired=TRUE)

subtypes <- unique(res[,"subtype"])
pvals <- matrix(0, ncol=4,nrow=3)
colnames(pvals) <- c("Primary_CS","Recurrent_CS","Clonal_PR","Subclonal_PR")
rownames(pvals) = subtypes
for(i in 1:length(subtypes))
{
	sub_res <- res[which(res[,"subtype"]==subtypes[i]),]
	
	#Clonal vs subclonal in each tumor
	pvals[i,1] <- wilcox.test(sub_res[,"nd_a_c"],sub_res[,"nd_a_s"],paired=TRUE)$p.value
	pvals[i,2] <- wilcox.test(sub_res[,"nd_b_c"],sub_res[,"nd_b_s"],paired=TRUE)$p.value

	#Clonality status by timepoint
	pvals[i,3] <- wilcox.test(sub_res[,"nd_a_c"],sub_res[,"nd_b_c"],paired=TRUE)$p.value
	pvals[i,4] <- wilcox.test(sub_res[,"nd_a_s"],sub_res[,"nd_b_s"],paired=TRUE)$p.value	
	
}

pair <- c(rep(res[,"tumor_pair_barcode"],4))
samp <- rep(c(res[,"tumor_barcode_a"],res[,"tumor_barcode_b"]),2)
rneo <- c(res[,"nd_a_c"],res[,"nd_b_c"],res[,"nd_a_s"],res[,"nd_b_s"])
subtype <- rep(res[,"subtype"],4)
timepoint <- rep(c(rep("Initial",nrow(res)),rep("Recurrent",nrow(res))),2)
clonality <- c(rep("Clonal",nrow(res)*2),rep("Subclonal",nrow(res)*2))
plot_res <- data.frame(pair,samp,rneo,subtype,timepoint,clonality)

pval <- rep("",nrow(plot_res))
pval[which(plot_res[,"timepoint"]=="Initial")] <- "P = 0.44"
pval[which(plot_res[,"timepoint"]=="Recurrent")] <- "P = 0.81"
p_text <- cbind(plot_res,pval)

gtsize = 7/(14/5)

pdf("/projects/varnf/GLASS/Figures/resubmission/clonality_ladderplot.pdf",width=3.25,height=2)
ggplot(plot_res, aes(y = rneo, x = clonality, group=pair, colour=subtype)) + 
	geom_line(size=0.45,alpha=0.4) +
	geom_point(size=1,colour="black") +
	facet_grid(.~timepoint) +
	labs(y="Observed/expected neoantigens") +
 	geom_text(data=p_text, aes(x=1.5, y=2, label = pval), size=gtsize, color = "black") +
	theme_classic() +
	theme(strip.text.x = element_text(size=7),strip.background = element_blank(),
	axis.title.y=element_text(size=7),axis.title.x=element_blank(),
	axis.text=element_text(colour="black",size=7), axis.text.x=element_text(hjust=0.5,size=7),
	legend.position="none")
dev.off()