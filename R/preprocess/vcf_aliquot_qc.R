library(VariantAnnotation)
library(DBI)
library(odbc)

rm(list=ls())
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")
q <- "SELECT * FROM biospecimen.aliquots"
aliquots <- dbGetQuery(con,q)
aliquots[,"case_barcode"] <- sapply(strsplit(aliquots[,"sample_barcode"],"-"),function(x)paste(x[1:3],collapse="-"))

myDir1 <- "/projects/verhaak-lab/GLASS-analysis/results/mutect2/m2filter"

mytag <- dir(myDir1)
mytag <- mytag[grep("filtered.vcf.gz$",mytag)]
vcff <- paste(myDir1,mytag,sep="/")
mytag <- gsub(".filtered.vcf.gz","",mytag)

check <- matrix(NA,nrow=length(vcff),ncol=15)
rownames(check) <- mytag
aliquot_match <- rep(0,length(vcff))
for(i in 1:length(vcff))
{
	cat("\r",i)
	vcf = readVcf(vcff[i], "hg19")
	samp_names <- rownames(colData(vcf))
	case_names <- sapply(strsplit(samp_names,"-"),function(x)paste(x[1:3],collapse="-"))
	
	samp_boo <- as.numeric(case_names == mytag[i])
	nsamp <- length(samp_boo)
	
	check[i,1:nsamp] <- samp_boo
	check[i,ncol(check)] <- nsamp
	
	sub_aliquots <- aliquots[which(aliquots[,"case_barcode"]==mytag[i]),]
	aliquot_match[i] <- sum(samp_names %in% sub_aliquots[,"aliquot_barcode"])/nrow(sub_aliquots)
}

sums <- apply(check[,1:(ncol(check)-1)],1,function(x)sum(x,na.rm=T))
sums == check[,ncol(check)]
sum(sums == check[,ncol(check)]) == nrow(check)
aliquot_match

