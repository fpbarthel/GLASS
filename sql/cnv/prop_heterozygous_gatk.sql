/*
- For each tumor pair in the tumor pairs table, compute the proportion of the original genome changed
*/
WITH
cnv AS
(
	SELECT
		gs.tumor_pair_barcode, 
		gs.case_barcode, 
		gs.tumor_barcode_a, 
		gs.tumor_barcode_b,
		sum(upper(pos) - lower(pos) -1) AS seg_size,
		sum(CASE WHEN gs.cnv_call = 0 THEN (upper(pos) - lower(pos) -1) ELSE 0 END) AS het_size
	FROM analysis.gatk_seg_diff_call gs
	WHERE chrom < 23
	GROUP BY 1,2,3,4
)
SELECT
	tumor_pair_barcode, 
	case_barcode, 
	tumor_barcode_a, 
	tumor_barcode_b,
	round(1.0 - het_size::decimal/seg_size,4) AS prop_change
FROM cnv
ORDER BY 2