library(DBI)
library(odbc)
library(ggplot2)
library(reshape)

rm(list=ls())

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

q = "WITH hla_table AS
(
	SELECT aliquot_barcode, hla_allele
	FROM analysis.neoantigens_by_aliquot
	GROUP BY aliquot_barcode, hla_allele
),
hla_tot AS
(
	SELECT aliquot_barcode, COUNT(*) AS hla_num
	FROM hla_table
	GROUP BY aliquot_barcode
)

SELECT nd.*,hla.hla_num
FROM analysis.neoantigen_depletion nd
INNER JOIN hla_tot hla ON hla.aliquot_barcode = nd.aliquot_barcode
ORDER BY rneo"

res <- dbGetQuery(con,q)

res[,"hla_num"] <- as.numeric(res[,"hla_num"])

cor.test(res[,"rneo"],res[,"hla_num"],method="s")		#R = 0.29	P = 5.1e-10

pri <- res[grep("-TP-",res[,1]),]
rec <- res[grep("-R1-|-R2-|-R3-|-R4-",res[,1]),]

cor.test(pri[,"rneo"],pri[,"hla_num"],method="s")		#R = 0.24	P = 4e-4
cor.test(rec[,"rneo"],rec[,"hla_num"],method="s")		#R = 0.34	P = 1.3e-7
