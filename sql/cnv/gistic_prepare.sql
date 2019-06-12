/*
Prepare a set of primaries and recurrences from the gold set (these are good quality CNV data)
As input for running GISTIC
==
Ensures one sample per patient
*/
WITH t1 AS
(
	SELECT
		aliquot_barcode,
		(CASE WHEN chrom = 23 THEN 'X' ELSE chrom::varchar(2) END) AS chrom,
		lower(pos) AS "start",
		upper(pos)-1 AS "end",
		num_points AS num_snps,
		log2_copy_ratio
	FROM variants.gatk_seg gs
)
SELECT t1.*,'P' AS sample_type FROM t1
INNER JOIN analysis.gold_set gs ON gs.tumor_barcode_a = t1.aliquot_barcode

UNION

SELECT t1.*,'R' AS sample_type FROM t1
INNER JOIN analysis.gold_set gs ON gs.tumor_barcode_b = t1.aliquot_barcode
