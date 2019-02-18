#This code calculates the fraction of mutations that give rise to neoantigens by looking
#dividing the neoantigen rate by the mutation rate
#Version 1
#Author: Fred Varn
#--------------------------------
library(DBI)
library(odbc)
library(ggplot2)
library(reshape)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

q = "SELECT sig.*, clin.surgical_interval,LOG(clin.surgical_interval) AS log_surg_int, sub.idh_codel_subtype
FROM analysis.abs_mutsig_fraction sig
INNER JOIN analysis.silver_set ss ON sig.tumor_pair_barcode = ss.tumor_pair_barcode
LEFT JOIN clinical.clinical_by_tumor_pair clin ON sig.tumor_pair_barcode = clin.tumor_pair_barcode
LEFT JOIN clinical.subtypes sub ON substring(sig.tumor_pair_barcode,1,12) = sub.case_barcode
WHERE sig.fraction = 'R' AND sig.signature = 1 AND clin.surgical_interval IS NOT NULL"

sig <- dbGetQuery(con,q)
sig [,"idh_codel_subtype"] <- factor(sig [,"idh_codel_subtype"], levels=c("IDHwt_noncodel","IDHmut_noncodel","IDHmut_codel"))

summary(lm(signature_score~log_surg_int+idh_codel_subtype,data=sig))
cor.test(sig[,"surgical_interval"],sig[,"signature_score"])