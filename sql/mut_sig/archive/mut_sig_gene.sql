/*
Calculate mutational signatures
By Gene
*/
WITH selected_genes AS
(
	SELECT DISTINCT gene_symbol
	FROM variants.anno
	ORDER BY 1
),
variant_contexts AS
(
	SELECT DISTINCT ref_context AS trinucleotide_context, alt
	FROM ref.signature_proba sp
),
variant_contexts_genes AS
(
	SELECT *
	FROM selected_genes, variant_contexts
),
variant_context_counts AS
(	
	SELECT gene_symbol, trinucleotide_context, pa.alt, COUNT(*) AS mut_n
	FROM variants.passgeno pg
	INNER JOIN variants.passanno pa ON pa.variant_id = pg.variant_id
	LEFT JOIN analysis.blocklist bl ON bl.aliquot_barcode = pg.aliquot_barcode
	WHERE ssm2_pass_call IS TRUE AND variant_type = 'SNP' AND ad_alt + ad_ref >= 15 AND fingerprint_exclusion = 'allow' AND coverage_exclusion = 'allow' AND variant_classification = 'MISSENSE'
	GROUP BY 1,2,3
),
variant_context_counts_genes AS
(
	SELECT vca.*, COALESCE(mut_n,0) AS mut_n, SUM(COALESCE(mut_n,0)) OVER (PARTITION BY vca.gene_symbol) AS mut_n_total
	FROM variant_contexts_genes vca
	LEFT JOIN variant_context_counts vcc ON vcc.gene_symbol = vca.gene_symbol AND vcc.trinucleotide_context = vca.trinucleotide_context AND vcc.alt = vca.alt
),
ref_context_array AS
(
	SELECT array_agg(a ORDER BY signature) AS ref_context_arr
	FROM (SELECT signature, array_agg(proba ORDER BY ref_context,alt) a FROM ref.signature_proba sp GROUP BY 1) t
),
context_reconstruction AS
(
	SELECT gene_symbol,ref_context_arr,sum(mut_n) AS mut_n, array_agg(mut_n ORDER BY trinucleotide_context,alt), lsqnonneg(ref_context_arr, array_agg(mut_n ORDER BY trinucleotide_context,alt)) AS mut_sigs
	FROM variant_context_counts_genes, ref_context_array
	WHERE mut_n_total > 9
	GROUP BY 1,2
)
SELECT gene_symbol, generate_series(1,30) AS signature, mut_n, unnest(mut_sigs) AS abs_score, UNNEST(mut_sigs) / (SELECT SUM(s) FROM UNNEST(mut_sigs) s) AS rel_score
FROM context_reconstruction