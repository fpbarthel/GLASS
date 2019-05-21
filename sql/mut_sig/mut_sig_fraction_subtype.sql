/*
Calculate mutational signatures
By subtype, fraction and hypermutation status
*/
WITH
selected_tumor_pairs AS
(
	SELECT ss.tumor_pair_barcode, ss.tumor_barcode_a, ss.tumor_barcode_b, idh_codel_subtype, hypermutator_status
	FROM analysis.gold_set ss
	INNER JOIN clinical.subtypes st ON st.case_barcode = ss.case_barcode
	INNER JOIN analysis.tumor_clinical_comparison tcc ON tcc.tumor_barcode_a = ss.tumor_barcode_a AND tcc.tumor_barcode_b = ss.tumor_barcode_b
),
selected_aliquots AS
(
	SELECT tumor_barcode_a AS aliquot_barcode, idh_codel_subtype, 'P' AS sample_type FROM selected_tumor_pairs
	UNION
	SELECT tumor_barcode_b AS aliquot_barcode, idh_codel_subtype, 'R' AS sample_type FROM selected_tumor_pairs
),
selected_subtypes AS
(
	SELECT DISTINCT idh_codel_subtype FROM selected_aliquots
),
variant_contexts AS
(
	SELECT DISTINCT ref_context AS trinucleotide_context, alt
	FROM ref.signature_proba sp
),
variant_contexts_aliquots AS
(
	SELECT *
	FROM variant_contexts
	CROSS JOIN selected_subtypes
	CROSS JOIN fractions
),
variant_context_counts AS
(	
	SELECT
		CASE
		WHEN pg.mutect2_call_a AND pg.mutect2_call_b THEN 'S'
		WHEN pg.mutect2_call_a AND NOT pg.mutect2_call_b THEN 'P'
		WHEN pg.mutect2_call_b AND NOT pg.mutect2_call_a THEN 'R'
		ELSE NULL
		END AS fraction,
		idh_codel_subtype,
		trinucleotide_context, pa.alt, COUNT(*) AS mut_n
	FROM variants.pgeno pg
	INNER JOIN selected_tumor_pairs sa ON sa.tumor_pair_barcode = pg.tumor_pair_barcode
	INNER JOIN variants.passanno pa ON pa.variant_id = pg.variant_id
	WHERE pg.comparison_type = 'longitudinal' AND pg.variant_type = 'SNP' AND 
		((pg.mutect2_call_a AND NOT pg.mutect2_call_b AND (pg.ref_count_a + pg.alt_count_a) > 14) OR 
		 (pg.mutect2_call_b AND NOT pg.mutect2_call_a AND (pg.ref_count_b + pg.alt_count_b) > 14) OR
         (pg.mutect2_call_a AND pg.mutect2_call_b AND (pg.ref_count_a + pg.alt_count_a) > 14 AND (pg.ref_count_b + pg.alt_count_b) > 14))
	GROUP BY 1,2,3,4
),
variant_context_counts_aliquots AS
(
	SELECT vca.*, COALESCE(mut_n,0) AS mut_n, SUM(COALESCE(mut_n,0)) OVER (PARTITION BY vca.fraction, vca.idh_codel_subtype) AS mut_n_total
	FROM variant_contexts_aliquots vca
	LEFT JOIN variant_context_counts vcc ON vcc.fraction = vca.fraction AND vcc.idh_codel_subtype = vca.idh_codel_subtype AND vcc.trinucleotide_context = vca.trinucleotide_context AND vcc.alt = vca.alt
),
ref_context_array AS
(
	SELECT array_agg(a ORDER BY signature) AS ref_context_arr
	FROM (SELECT signature, array_agg(proba ORDER BY ref_context,alt) a FROM ref.signature_proba sp GROUP BY 1) t
),
context_reconstruction AS
(
	SELECT fraction,idh_codel_subtype,ref_context_arr,sum(mut_n) AS mut_n, array_agg(mut_n ORDER BY trinucleotide_context,alt), lsqnonneg(ref_context_arr, array_agg(mut_n ORDER BY trinucleotide_context,alt)) AS mut_sigs
	FROM variant_context_counts_aliquots, ref_context_array
	WHERE mut_n_total > 1
	GROUP BY 1,2,3
)
--SELECT fraction, idh_codel_subtype, hypermutator_status, generate_series(1,30) AS signature, mut_n, unnest(mut_sigs) AS abs_score, UNNEST(mut_sigs) / (SELECT SUM(s) FROM UNNEST(mut_sigs) s) AS rel_score
--FROM context_reconstruction
SELECT * FROM variant_context_counts_aliquots