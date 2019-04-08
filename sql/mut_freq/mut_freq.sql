/*
Compute mutation frequencies for each aliquot_barcode
- Mutation frequencies is output in mutations per megabase (1e6 basepairs)
- Only mutations with >= 15x are counted
- Mutation counts are divided by the number of basepairs with at least 15x coverage
- COALESCE is used to prevent "divison by zero" problems
- JOIN to blocklist so we don't report any aliquots that were excluded based on fingerprinting or coverage

Note:
- for the ssm2_count table I counted mutation using greater than (>) --> 14 threshold
- for the coverage table I counted coverage using greater than or equal to (>=) --> 15 threshold
*/

SELECT 
	m2.aliquot_barcode,
	cumulative_coverage,
	ssm2_call_count AS mutation_count,
	COALESCE(ROUND(ssm2_call_count::numeric / cumulative_coverage::numeric * 1e6, 4), 0::numeric) AS coverage_adj_mut_freq
FROM variants.ssm2_count m2
INNER JOIN analysis.coverage cov ON cov.aliquot_barcode = m2.aliquot_barcode
WHERE m2.ad_depth = 14 AND cov.coverage = 15

