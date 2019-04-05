/*
Calculate mutational signatures
By Aliquot
*/
WITH selected_aliquots AS
(
	SELECT ba.aliquot_barcode
	FROM biospecimen.aliquots ba
	LEFT JOIN analysis.blocklist bl ON bl.aliquot_barcode = ba.aliquot_barcode
	LEFT JOIN biospecimen.samples bs ON bs.sample_barcode = ba.sample_barcode
	WHERE fingerprint_exclusion = 'allow' AND coverage_exclusion = 'allow' AND sample_type NOT IN ('NM','NB') --AND ba.aliquot_barcode IN ('GLSS-HK-0002-R1-01D-WGS-S3QETN') --,'GLSS-DK-0008-R1-01D-WXS-DDD4B8','TCGA-06-0190-R1-01D-WGS-P20F5P','TCGA-14-1402-R1-01D-WGS-2EHMQ2')--('GLSS-CU-R008-TP-01D-WXS-0238UJ','GLSS-HK-0004-R1-01D-WGS-RYFPEB')
),
variant_contexts AS
(
	SELECT DISTINCT ref_context AS trinucleotide_context, alt
	FROM ref.signature_proba sp
),
variant_contexts_aliquots AS
(
	SELECT *
	FROM selected_aliquots, variant_contexts
),
variant_context_counts AS
(	
	SELECT aliquot_barcode, trinucleotide_context, pa.alt, COUNT(*) AS mut_n
	FROM variants.passgeno pg
	INNER JOIN variants.passanno pa ON pa.variant_id = pg.variant_id
	WHERE ssm2_pass_call IS TRUE AND variant_type = 'SNP' AND ad_alt + ad_ref >= 15
	GROUP BY 1,2,3
),
variant_context_counts_aliquots AS
(
	SELECT vca.*, COALESCE(mut_n,0) AS mut_n, SUM(COALESCE(mut_n,0)) OVER (PARTITION BY vca.aliquot_barcode) AS mut_n_total
	FROM variant_contexts_aliquots vca
	LEFT JOIN variant_context_counts vcc ON vcc.aliquot_barcode = vca.aliquot_barcode AND vcc.trinucleotide_context = vca.trinucleotide_context AND vcc.alt = vca.alt
),
ref_context_array AS
(
	SELECT array_agg(a ORDER BY signature) AS ref_context_arr
	FROM (SELECT signature, array_agg(proba ORDER BY ref_context,alt) a FROM ref.signature_proba sp GROUP BY 1) t
),
context_reconstruction AS
(
	SELECT aliquot_barcode,ref_context_arr,sum(mut_n) AS mut_n, array_agg(mut_n ORDER BY trinucleotide_context,alt), lsqnonneg(ref_context_arr, array_agg(mut_n ORDER BY trinucleotide_context,alt)) AS mut_sigs
	FROM variant_context_counts_aliquots, ref_context_array
	WHERE mut_n_total > 1
	GROUP BY 1,2
)
SELECT aliquot_barcode, generate_series(1,30) AS signature, mut_n, unnest(mut_sigs) AS abs_score, UNNEST(mut_sigs) / (SELECT SUM(s) FROM UNNEST(mut_sigs) s) AS rel_score
FROM context_reconstruction