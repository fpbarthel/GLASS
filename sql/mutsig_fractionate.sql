WITH paired_fractions AS
(
	SELECT tp.tumor_pair_barcode, fraction
	FROM analysis.tumor_pairs tp, fractions
),
variant_contexts AS
(
	SELECT DISTINCT ref_context AS trinucleotide_context, alt
	FROM ref.signature_proba sp
),
variant_contexts_pairs AS
(
	SELECT *
	FROM paired_fractions, variant_contexts
),
variant_contexts_counts AS
(	
	SELECT tumor_pair_barcode, 
	(CASE WHEN mutect2_call_a AND mutect2_call_b THEN 'S' WHEN mutect2_call_a AND NOT mutect2_call_b THEN 'P' WHEN mutect2_call_b AND NOT mutect2_call_a THEN 'R' END) AS fraction,
	trinucleotide_context, pa.alt, COUNT(*) AS mut_n
	FROM variants.pgeno pg
	INNER JOIN variants.passanno pa ON pa.variant_id = pg.variant_id
	WHERE pg.comparison_type = 'longitudinal' AND
		pg.variant_type = 'SNP' AND
		(pg.mutect2_call_a  AND (pg.ref_count_a + pg.alt_count_a) > 14) OR  (pg.mutect2_call_b AND (pg.ref_count_b + pg.alt_count_b) > 14)
	GROUP BY 1,2,3,4
),
variant_contexts_counts_pairs AS
(
	SELECT vcp.*, COALESCE(mut_n,0) AS mut_n, SUM(COALESCE(mut_n,0)) OVER (PARTITION BY vcp.tumor_pair_barcode, vcp.fraction) AS mut_n_total
	FROM variant_contexts_pairs vcp
	LEFT JOIN variant_contexts_counts vcc ON vcc.tumor_pair_barcode = vcp.tumor_pair_barcode AND vcc.fraction = vcp.fraction AND vcc.trinucleotide_context = vcp.trinucleotide_context AND vcc.alt = vcp.alt
),
ref_context_array AS
(
	SELECT array_agg(a ORDER BY signature) AS ref_context_arr
	FROM (SELECT signature, array_agg(proba ORDER BY ref_context,alt) a FROM ref.signature_proba sp GROUP BY 1) t
),
context_reconstruction AS
(
	SELECT tumor_pair_barcode,fraction,ref_context_arr,sum(mut_n) AS mut_n, array_agg(mut_n ORDER BY trinucleotide_context,alt), lsqnonneg(ref_context_arr, array_agg(mut_n ORDER BY trinucleotide_context,alt)) AS mut_sigs
	FROM variant_contexts_counts_pairs, ref_context_array
	WHERE mut_n_total > 1
	GROUP BY 1,2,3
)
SELECT tumor_pair_barcode, fraction, generate_series(1,30) AS signature, mut_n, unnest(mut_sigs) AS abs_score, UNNEST(mut_sigs) / (SELECT SUM(s) FROM UNNEST(mut_sigs) s) AS rel_score
FROM context_reconstruction
ORDER BY 1, 2, 3																																	   