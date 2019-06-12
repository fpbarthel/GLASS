/*
Compute mutation, neoantigen, and strong neoantigen frequencies for each aliquot_barcode
- Frequencies are output in mutations/neoantigens per megabase (1e6 basepairs)
- Only mutations/neoantigens with >= 15x are counted
- Mutations/neoantigen counts are divided by the number of basepairs with at least 15x coverage
- Strong neoantigens represent neoantigens with a binding affinity score < 50
- COALESCE is used to prevent "divison by zero" problems
- JOIN to blocklist so we don't report any aliquots that were excluded based on fingerprinting or coverage
- Currently INNER JOIN-ing to silver_set to avoid samples used in multiple comparisons getting counted multiple times- not ideal and should be fixed later
- DO NOT USE THESE NUMBERS FOR WGS. Neoantigens are only found in exons so the coverage-adjusted neoantigen rate is SEVERELY underestimated
*/

WITH mutation_counts AS
(
	SELECT 
		m2.aliquot_barcode,
		cumulative_coverage,
		ssm2_call_count AS mutation_count,
		COALESCE(ROUND(ssm2_call_count::numeric / cumulative_coverage::numeric * 1e6, 4), 0::numeric) AS coverage_adj_mut_freq
	FROM variants.ssm2_count m2
	INNER JOIN analysis.coverage cov ON cov.aliquot_barcode = m2.aliquot_barcode
	WHERE m2.ad_depth = 14 AND cov.coverage = 15
	ORDER BY coverage_adj_mut_freq DESC
),
neoag_table AS
(
	SELECT 
		gt.tumor_pair_barcode,gt.chrom,gt.pos,gt.alt,gt.gene_symbol,gt.tumor_barcode_a,gt.tumor_barcode_b,alt_count_a, alt_count_b, ref_count_a, ref_count_b,netmhcpan_mt_score,row_number() OVER w AS neoag_number,
		(CASE WHEN mutect2_call_a AND mutect2_call_b THEN 'S' WHEN mutect2_call_a AND NOT mutect2_call_b THEN 'P' WHEN mutect2_call_b AND NOT mutect2_call_a THEN 'R' END) AS fraction,
		(CASE WHEN netmhcpan_mt_score IS NULL THEN 0 WHEN netmhcpan_mt_score IS NOT NULL THEN 1 END) AS neoantigen,
		(CASE WHEN netmhcpan_mt_score IS NULL OR netmhcpan_mt_score IS NOT NULL AND netmhcpan_mt_score > 50 THEN 0 WHEN netmhcpan_mt_score IS NOT NULL AND netmhcpan_mt_score <= 50 THEN 1 END) AS strong_neoantigen
	FROM variants.pgeno gt
	LEFT JOIN analysis.neoantigens_by_pair neo ON 
		neo.tumor_pair_barcode = gt.tumor_pair_barcode AND 
		neo.chrom = CAST(gt.chrom AS int) AND 
		neo.pos = gt.pos AND
		neo.alt = gt.alt
	INNER JOIN analysis.silver_set ss ON ss.tumor_pair_barcode = neo.tumor_pair_barcode
	WHERE ((mutect2_call_a  AND (ref_count_a + alt_count_a) > 14) OR  (mutect2_call_b AND (ref_count_b + alt_count_b) > 14)) AND netmhcpan_mt_score IS NOT NULL
	WINDOW w AS (PARTITION BY gt.tumor_pair_barcode,gt.chrom,gt.pos,gt.alt ORDER BY netmhcpan_mt_score)
),
neoag_counts AS
(
	SELECT tumor_barcode_a AS aliquot_barcode, SUM(neoantigen) AS neoag_count, SUM(strong_neoantigen) AS strong_neoag_count
	FROM neoag_table
	WHERE (fraction = 'P' OR fraction = 'S') AND neoag_number = 1
	GROUP BY tumor_barcode_a
	UNION
	SELECT tumor_barcode_b AS aliquot_barcode, SUM(neoantigen) AS neoag_count, SUM(strong_neoantigen) AS strong_neoag_count
	FROM neoag_table
	WHERE (fraction = 'R' OR fraction = 'S') AND neoag_number = 1
	GROUP BY tumor_barcode_b
)
SELECT
	cov.aliquot_barcode,
	cov.cumulative_coverage,
	mut.mutation_count,
	neo.neoag_count,
	neo.strong_neoag_count,
	COALESCE(ROUND(mut.mutation_count::numeric / cov.cumulative_coverage::numeric * 1e6, 4), 0::numeric) AS coverage_adj_mut_freq,
	COALESCE(ROUND(neo.neoag_count::numeric / cov.cumulative_coverage::numeric * 1e6, 4), 0::numeric) AS coverage_adj_neoag_freq,
	COALESCE(ROUND(neo.strong_neoag_count::numeric / cov.cumulative_coverage::numeric * 1e6, 4), 0::numeric) AS coverage_adj_strong_neoag_freq
INTO analysis.neoag_freq
FROM analysis.coverage cov
LEFT JOIN mutation_counts mut ON mut.aliquot_barcode = cov.aliquot_barcode
LEFT JOIN analysis.blocklist bl ON bl.aliquot_barcode = cov.aliquot_barcode
INNER JOIN neoag_counts neo ON neo.aliquot_barcode = cov.aliquot_barcode
WHERE
	cov.coverage = 15 AND
	bl.fingerprint_exclusion = 'allow' AND
	bl.coverage_exclusion = 'allow' AND
	mutation_count IS NOT NULL
ORDER BY 6 DESC

