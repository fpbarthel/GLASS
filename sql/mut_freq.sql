/*
Compute mutation frequencies for each aliquot_barcode
- Mutation frequencies is output in mutations per megabase (1e6 basepairs)
- Only mutations with >= 15x are counted
- Mutation counts are divided by the number of basepairs with at least 15x coverage
- COALESCE is used to prevent "divison by zero" problems
- JOIN to blocklist so we don't report any aliquots that were excluded based on fingerprinting or coverage
*/

WITH mutation_counts AS
(
	SELECT
		aliquot_barcode,
		COUNT(*) AS mutation_count
	FROM analysis.genotypes
	WHERE mutect2_call AND (alt_count + ref_count) > 14
	GROUP BY aliquot_barcode
)

SELECT
	cov.aliquot_barcode,
	cov.cumulative_coverage,
	mut.mutation_count,
	COALESCE(ROUND(mut.mutation_count::numeric / cov.cumulative_coverage::numeric * 1e6, 4), 0::numeric) AS coverage_adj_mut_freq
FROM analysis.coverage cov
LEFT JOIN mutation_counts mut ON mut.aliquot_barcode = cov.aliquot_barcode
LEFT JOIN analysis.blocklist bl ON bl.aliquot_barcode = cov.aliquot_barcode
WHERE
	cov.coverage = 15 AND
	bl.fingerprint_exclusion = 'allow' AND
	bl.coverage_exclusion = 'allow' AND
	mutation_count IS NOT NULL
ORDER BY 4 DESC