#This script saves Supplementary Table 6 (generated using the neoantigen_peptide_counts.sql query) to a text file
#-----------------------------------------------------

library(DBI)
library(odbc)
library(ggplot2)
library(reshape)

rm(list=ls())

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

res <- dbGetQuery(con, read_file("sql/neoag/neoantigen_peptide_counts.sql"))

write.table(res,"/projects/varnf/GLASS/Figures/resubmission/final/SuppTableS4.txt",sep="\t",quote=F,row.names=F)