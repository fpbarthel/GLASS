#Version 1
#Author: Fred Varn
#--------------------------------
library(DBI)
library(odbc)
library(ggplot2)
library(reshape)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

q = "SELECT mf.*,sub.case_age_diagnosis_years
FROM analysis.mut_sig_aliquot mf 
INNER JOIN analysis.silver_set ss ON mf.aliquot_barcode = ss.tumor_barcode_a
LEFT JOIN clinical.cases sub ON substring(mf.aliquot_barcode,1,12) = sub.case_barcode"

sig <- dbGetQuery(con,q)

sigs <- unique(sig[,"signature"])
pearson <- p.val <- rep(0,length(sigs))
for(i in 1:length(sigs))
{
	sub_sigs <- sig[which(sig[,"signature"]==sigs[i]),]
	pearson[i] <- cor(sub_sigs[,"abs_score"],sub_sigs[,"case_age_diagnosis_years"],method="p",use="pairwise.complete.obs")
	p.val[i] <- cor.test(sub_sigs[,"abs_score"],sub_sigs[,"case_age_diagnosis_years"],method="p")$p.value
}
q.val <- p.adjust(p.val)
res <- data.frame(pearson,p.val,q.val)