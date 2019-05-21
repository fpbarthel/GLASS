library(DBI)
library(odbc)

rm(list=ls())

cibersort <- read.delim("/projects/varnf/GLASS/data/CIBERSORT/CIBERSORT_cancer_cell.txt",sep="\t",header=T,stringsAsFactor=F)
mapping_table <- read.delim("/projects/varnf/GLASS/data/CIBERSORT/cancer_cell_RNAseq_mapping.txt",sep="\t",header=T,stringsAsFactor=F)
myoutf <- "/projects/varnf/GLASS/data/CIBERSORT/CIBERSORT_GLASS_format.txt"

#Add information for TCGA-14-1402
mapping_table[which(mapping_table[,"SampleId"]=="TCGA.14.1402.01"),"GLSS_barcodeTP"] <- "TCGA-14-1402-TP"
mapping_table[which(mapping_table[,"SampleId2"]=="TCGA.14.1402.02A"),"GLSS_barcodeR1"] <- "TCGA-14-1402-R1"

cibersort[,2] <- gsub("-",".",cibersort[,2])

mapping <- c(mapping_table[,"GLSS_barcodeTP"],mapping_table[,"GLSS_barcodeR1"])
names(mapping) <- c(mapping_table[,"SampleId"],mapping_table[,"SampleId2"])

ordered_names <- mapping[cibersort[,"SampleId"]]
cibersort[,"sample_barcode"] <- ordered_names

cibersort <- cibersort[-which(is.na(cibersort[,"sample_barcode"])),]
cibersort <- cibersort[-which(cibersort[,"sample_barcode"]=="not in data freeze"),]

rownames(cibersort) <- cibersort[,"sample_barcode"]
cibersort <- cibersort[,3:24]
cibersort <- cbind(rownames(cibersort),cibersort)
colnames(cibersort)[1] <- "sample_barcode"

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")
dbWriteTable(con, Id(schema="analysis",table="cibersort"), cibersort, overwrite=TRUE, row.names=FALSE)

write.table(cibersort, myoutf,sep="\t",quote=F,row.names=F)
