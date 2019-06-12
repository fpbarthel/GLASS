library(odbc)
library(DBI)
library(ggplot2)
library(tidyverse)
library(reshape)
library(plyr)

rm(list=ls())

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")  

res <- dbGetQuery(con, read_file("/projects/varnf/GLASS/GLASS/sql/neoag/cibersort_depletion.sql"))
colnames(res)[6:16] <- c("B cells","PCs","CD8","CD4","gd T cells","NKs","Mono/MFs","DCs","MCs","Eos","PMN cells")

#Code to select paired samples only
myrle <- rle(substr(res[,"aliquot_barcode"],1,12))
pull <- myrle[["values"]][which(myrle[["lengths"]]==2)]
res <- res[which(substr(res[,"aliquot_barcode"],1,12)%in%pull),]

ini_res <- res[grep("-TP-",res[,"aliquot_barcode"]),]
rec_res <- res[grep('-R1-|-R2-|-R3-',res[,"aliquot_barcode"]),]

g1 <- which(ini_res[,"nd"]<=1)
g2 <- which(ini_res[,"nd"]>1)

apply(ini_res[,6:16],2,function(x)wilcox.test(x[g1],x[g2])$p.value)
apply(ini_res[,6:16],2,function(x)(median(x[g1])-median(x[g2])))

g1 <- which(rec_res[,"nd"]<=1)
g2 <- which(rec_res[,"nd"]>1)

apply(rec_res[,6:16],2,function(x)wilcox.test(x[g1],x[g2])$p.value)
apply(rec_res[,6:16],2,function(x)(median(x[g1])-median(x[g2])))

ini1 <- res[which(res[,"timepoint"]=="Initial" & res[,"nd"]<1),]
ini2 <- res[which(res[,"timepoint"]=="Initial" & res[,"nd"]>=1),]

rec1 <- res[which(res[,"timepoint"]=="Recurrent" & res[,"nd"]<1),]
rec2 <- res[which(res[,"timepoint"]=="Recurrent" & res[,"nd"]>=1),]

ini1_mean <- apply(ini1[,6:16],2,mean)
ini1_mean <- ini1_mean/(sum(ini1_mean))
ini2_mean <- apply(ini2[,6:16],2,mean)
ini2_mean <- ini2_mean/(sum(ini2_mean))

rec1_mean <- apply(rec1[,6:16],2,mean)
rec1_mean <- rec1_mean/(sum(rec1_mean))
rec2_mean <- apply(rec2[,6:16],2,mean)
rec2_mean <- rec2_mean/(sum(rec2_mean))

plot_res <- rbind(ini1_mean,ini2_mean,rec1_mean,rec2_mean)
plot_res <- melt(plot_res)

rec <- rep(rep(c("Initial","Recurrent"),each=2),nrow(plot_res)/4)
dep <- rep(c("< 1","> 1"),nrow(plot_res)/2)

plot_res <- cbind(plot_res,rec,dep)
colnames(plot_res) <- c("group","cell","proportion","recurrence","depletion")
plot_res[,"cell"] <- factor(plot_res[,"cell"],levels=rev(unique(plot_res[,"cell"])))
plot_res[,"sig"] <- rep("",nrow(plot_res))
plot_res[which(plot_res[,"recurrence"]=="Recurrent" & plot_res[,"cell"]=="Mono/MFs"),"sig"] <- "*"
plot_res[which(plot_res[,"recurrence"]=="Recurrent" & plot_res[,"cell"]=="PMN cells"),"sig"] <- "*"
plot_res <- ddply(plot_res,c("group","recurrence","depletion"),transform,ypos=cumsum(proportion))

#CIBERSORT stacked barplot
colors <- rev(c("#435D71","#52B473","#31C2EA","#EBC646","#D74729","#886B9D","#EBDBC4","#A9D4E9","#5EA6C2","#938DBE","#BDBAD1"))

gtsize = 7/(14/5)
pdf("/projects/varnf/GLASS/Figures/resubmission/CIBERSORT_stacked_barplot.pdf",width=2.2,height=2)
p1 <- ggplot(plot_res, aes(y = proportion, x = depletion, fill = cell)) +
geom_bar(stat="identity",width=1,colour="black",size=0.15) +
geom_text(aes(y=ypos,label=sig),vjust=1.6,colour="black",size=gtsize) +
facet_wrap(~recurrence) +
scale_fill_manual(values=colors) +
theme(
strip.text.x = element_text(size=7),
strip.background = element_blank(),
legend.title=element_blank(),
legend.text=element_text(size=7),
legend.key.size = unit(0.35, "cm"),
panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(),
plot.background=element_blank(),
axis.text.x = element_text(hjust=0.5),
axis.title=element_blank()) +
scale_x_discrete(expand=c(0,0)) +
scale_y_continuous(limits=c(0,1.0001),expand=c(0,0))
p1
dev.off()
