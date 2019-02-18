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

q = "SELECT sig.*,sub.case_age_diagnosis_years
FROM analysis.abs_mutsig_sample sig
INNER JOIN analysis.silver_set ss ON sig.aliquot_barcode = ss.tumor_barcode_a
LEFT JOIN clinical.cases sub ON substring(sig.aliquot_barcode,1,12) = sub.case_barcode
WHERE sig.signature = 1"

sig <- dbGetQuery(con,q)

cor.test(sig[,"signature_score"],sig[,"case_age_diagnosis_years"],method="p")