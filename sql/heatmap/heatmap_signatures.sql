/*SELECT gs.case_barcode, idh_codel_subtype, fraction, signature, mut_n, abs_score, rel_score
FROM analysis.mut_sig_fraction ms
INNER JOIN analysis.gold_set gs ON gs.tumor_pair_barcode = ms.tumor_pair_barcode
INNER JOIN clinical.subtypes st ON st.case_barcode = gs.case_barcode
WHERE signature IN (1,3,11,15,26)*/
WITH t1 AS (
	SELECT gs.case_barcode, idh_codel_subtype, fraction, signature, mut_n, abs_score, rel_score, RANK() OVER (PARTITION BY gs.case_barcode,fraction ORDER BY rel_score DESC) AS rnk
	FROM analysis.mut_sig_fraction_limited ms
	INNER JOIN analysis.gold_set gs ON gs.tumor_pair_barcode = ms.tumor_pair_barcode
	INNER JOIN clinical.subtypes st ON st.case_barcode = gs.case_barcode
	)
	SELECT * FROM t1 --WHERE rnk = 1