#Code that makes boxplots of observed-to-expected neoantigen ratio by fraction and subtype
#Did not make it into the original manuscript, have not thoroughly reviewed this code
#May be useful as an Extended Data Figure down the road
#-----------------------------------------------------

library(odbc)
library(DBI)
library(ggplot2)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")  

q <- "SELECT gs.tumor_pair_barcode, nd.fraction, nd.rneo, clin.idh_codel_subtype AS subtype,
CASE WHEN mf.coverage_adj_mut_freq >= 10 THEN 1 WHEN mf.coverage_adj_mut_freq < 10 THEN 0 END AS hm_rec
FROM analysis.gold_set gs
INNER JOIN analysis.neoantigen_depletion_fraction nd ON gs.tumor_pair_barcode = nd.tumor_pair_barcode
LEFT JOIN analysis.mut_freq mf ON mf.aliquot_barcode = gs.tumor_barcode_b
LEFT JOIN clinical.subtypes clin ON clin.case_barcode = gs.case_barcode
WHERE rneo IS NOT NULL
ORDER BY nd.rneo DESC"

res <- dbGetQuery(con, q)

myp <- mye <- rep(0,12)

g1 <- res[which(res[,"fraction"]=="P"),"rneo"]
g2 <- res[which(res[,"fraction"]=="R"),"rneo"]
g3 <- res[which(res[,"fraction"]=="S"),"rneo"]

myp[1] <- wilcox.test(g1, mu=1)$p.value
myp[2] <- wilcox.test(g2, mu=1)$p.value
myp[3] <- wilcox.test(g3, mu=1)$p.value

mye[1] <- median(g1)
mye[2] <- median(g2)
mye[3] <- median(g3)

subtypes <- unique(res[,"subtype"])
mysub <- c(rep("all",3),rep(subtypes,each=3))
myfrac <- rep(c("P","R","S"),4)
for(i in 1:length(subtypes))
{
	subres <- res[which(res[,"subtype"]==subtypes[i]),]
	g1 <- subres[which(subres[,"fraction"]=="P"),"rneo"]
	g2 <- subres[which(subres[,"fraction"]=="R"),"rneo"]
	g3 <- subres[which(subres[,"fraction"]=="S"),"rneo"]

	myp[(i*3)+1] <- wilcox.test(g1, mu=1)$p.value
	myp[(i*3)+2] <- wilcox.test(g2, mu=1)$p.value
	myp[(i*3)+3] <- wilcox.test(g3, mu=1)$p.value
	
	mye[(i*3)+1] <- median(g1)
	mye[(i*3)+2] <- median(g2)
	mye[(i*3)+3] <- median(g3)
}
substats <- data.frame(mysub,myfrac,myp,mye)

pdf("/projects/varnf/GLASS/Figures/resubmission/depletion_subtype_fraction.pdf",width=3.5,height=2)
ggplot(res, aes(y = rneo, x = fraction, fill=fraction)) + 
	geom_hline(yintercept=1, color = "gray50", size=1) +
	geom_boxplot(lwd=0.5,outlier.size=0.1,fatten=2) +
	labs(y="observed/expected neoantigens") +
	scale_fill_manual(values=c("#CA2F66","#2FB3CA","#CA932F")) +
	facet_grid(~subtype, scales = "free") + 
	theme_classic() +
	theme(strip.text.x = element_text(size=8),strip.background = element_blank(),
	axis.title.y=element_text(size=8),axis.title.x=element_blank(),
	axis.text=element_text(colour="black",size=8), axis.text.x=element_text(hjust=0.5,size=8),
	legend.position="none")
dev.off()

hm <- unique(res[,"hm_rec"])
myhm <- c(rep(0,3),rep(1,3))
myfrac <- rep(c("P","R","S"),2)
myp <- mye <- rep(0,6)
for(i in 1:length(hm))
{
	hmres <- res[which(res[,"hm_rec"]==hm[i]),]
	g1 <- res[which(res[,"fraction"]=="P"),"rneo"]
	g2 <- res[which(res[,"fraction"]=="R"),"rneo"]
	g3 <- res[which(res[,"fraction"]=="S"),"rneo"]

	myp[((i-1)*3)+1] <- wilcox.test(g1, mu=1)$p.value
	myp[((i-1)*3)+2] <- wilcox.test(g2, mu=1)$p.value
	myp[((i-1)*3)+3] <- wilcox.test(g3, mu=1)$p.value

	mye[((i-1)*3)+1] <- median(g1)
	mye[((i-1)*3)+2] <- median(g2)
	mye[((i-1)*3)+3] <- median(g3)
}
hmstats <- data.frame(myhm,myfrac,myp,mye)


pdf("/projects/varnf/GLASS/Figures/resubmission/depletion_hypermutation_fraction.pdf",width=2.45,height=2)
ggplot(res, aes(y = rneo, x = fraction, fill=fraction)) + 
	geom_hline(yintercept=1, color = "gray50", size=1) +
	geom_boxplot(lwd=0.5,outlier.size=0.1,fatten=2) +
	labs(y="observed/expected neoantigens") +
	scale_fill_manual(values=c("#CA2F66","#2FB3CA","#CA932F")) +
	facet_grid(~hm_rec, scales = "free") + 
	theme_classic() +
	theme(strip.text.x = element_text(size=8),strip.background = element_blank(),
	axis.title.y=element_text(size=8),axis.title.x=element_blank(),
	axis.text=element_text(colour="black",size=8), axis.text.x=element_text(hjust=0.5,size=8),
	legend.position="none")
dev.off()