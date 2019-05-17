WITH t1 AS
(
	SELECT gs.case_barcode, idh_codel_subtype, hypermutator_status::integer, fraction, signature, mut_n, abs_score, rel_score, RANK() OVER (PARTITION BY gs.case_barcode,fraction ORDER BY rel_score DESC) AS rnk
	FROM analysis.mut_sig_fraction_limited ms
	INNER JOIN analysis.gold_set gs ON gs.tumor_pair_barcode = ms.tumor_pair_barcode
	INNER JOIN analysis.tumor_clinical_comparison tcc ON tcc.tumor_pair_barcode = ms.tumor_pair_barcode
	INNER JOIN clinical.subtypes st ON st.case_barcode = gs.case_barcode
)
SELECT * --signature,fraction,hypermutator_status,sum(rel_score) / count(rel_score) AS avg_rel_score, stddev(rel_score) AS sd_rel_score
FROM t1
--GROUP BY 1,2,3
--ORDER BY 4 DESC