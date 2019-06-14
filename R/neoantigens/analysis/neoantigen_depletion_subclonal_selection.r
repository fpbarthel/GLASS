#Code to that compares each the observed-to-expected neoantigen ratios between samples marked as "selected" and "neutral" evolution using the subclonalSelection method
#Query at the top joins analysis.neoantigen_depletion table to analysis.subclonalselection table (and others)
#Comparisons are made for initial tumors and for recurrent tumors
#No significant associations
#This analysis is reported in the manuscript
#------------------------------------------------------------------------------

library(odbc)
library(DBI)
library(ggplot2)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")  

q <- "SELECT gs.* , nd1.rneo AS nd_a, nd2.rneo AS nd_b, clin.idh_codel_subtype AS subtype,sc1.most_probable_classification AS neut1, sc2.most_probable_classification AS neut2, sc1.probability_neutral AS prob1, sc2.probability_neutral AS prob2
FROM analysis.gold_set gs
LEFT JOIN analysis.neoantigen_depletion nd1 ON nd1.aliquot_barcode = gs.tumor_barcode_a
LEFT JOIN analysis.neoantigen_depletion nd2 ON nd2.aliquot_barcode = gs.tumor_barcode_b
LEFT JOIN analysis.mut_freq mf1 ON mf1.aliquot_barcode = gs.tumor_barcode_a
LEFT JOIN analysis.mut_freq mf2 ON mf2.aliquot_barcode = gs.tumor_barcode_b
LEFT JOIN clinical.subtypes clin ON clin.case_barcode = gs.case_barcode
LEFT JOIN analysis.subclonalselection sc1 ON sc1.aliquot_barcode = gs.tumor_barcode_a
LEFT JOIN analysis.subclonalselection sc2 ON sc2.aliquot_barcode = gs.tumor_barcode_b
WHERE nd1.rneo IS NOT NULL AND nd2.rneo IS NOT NULL AND (nd1.nobs >= 3 AND nd2.nobs >= 3) --AND
--sc1.most_probable_classification IS NOT NULL AND sc2.most_probable_classification IS NOT NULL
ORDER BY nd1.rneo"

res <- dbGetQuery(con, q)

subtypes <- unique(res[,"subtype"])
for(i in 1:length(subtypes))
{
	sub_res <- res[which(res[,"subtype"]==subtypes[i]),]
	pri_s <- sub_res[which(sub_res[,"neut1"]=="S"),"nd_a"]
	pri_n <- sub_res[which(sub_res[,"neut1"]=="N"),"nd_a"]
	
	rec_s <- sub_res[which(sub_res[,"neut2"]=="S"),"nd_b"]
	rec_n <- sub_res[which(sub_res[,"neut2"]=="N"),"nd_b"]
	
	wilcox.test(pri_s,pri_n)
	wilcox.test(rec_s,rec_n)
	
	s1 <- c(pri_s,rec_s)
	n1 <- c(pri_n,rec_n)
	wilcox.test(s1,n1)	
}
