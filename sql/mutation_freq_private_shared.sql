WITH selected_tumor_pairs AS
(
	SELECT
 		tumor_pair_barcode,
		row_number() OVER (PARTITION BY case_barcode ORDER BY surgical_interval_mo DESC, portion_a ASC, portion_b ASC, substring(tumor_pair_barcode from 27 for 3) ASC) AS priority
	FROM analysis.tumor_pairs ps
	LEFT JOIN analysis.blocklist b1 ON b1.aliquot_barcode = ps.tumor_barcode_a
	LEFT JOIN analysis.blocklist b2 ON b2.aliquot_barcode = ps.tumor_barcode_b
	WHERE
		comparison_type = 'longitudinal' AND
		sample_type_b <> 'M1' AND
		b1.coverage_exclusion = 'allow' AND b2.coverage_exclusion = 'allow'
)
SELECT
	tmc.tumor_pair_barcode,
	tmc.case_barcode,
	tmc.surgical_interval_mo,
	tmc.count_a,
	tmc.count_b,
	tmc.union_ab,
	tmc.intersection_ab,
	tmc.setdiff_a,
	tmc.setdiff_b,
	mf1.cumulative_coverage AS cov_a,
	mf2.cumulative_coverage AS cov_b,
	LEAST(mf1.cumulative_coverage, mf2.cumulative_coverage) AS min_cov,
	ROUND(setdiff_a::decimal / mf1.cumulative_coverage * 1e6, 4) AS mf_private_a,
	ROUND(setdiff_b::decimal / mf2.cumulative_coverage * 1e6, 4) AS mf_private_b,
	ROUND(intersection_ab::decimal / LEAST(mf1.cumulative_coverage, mf2.cumulative_coverage) * 1e6, 4) AS mf_shared,
	mf1.coverage_adj_mut_freq AS mf_a,
	mf2.coverage_adj_mut_freq AS mf_b
FROM analysis.tumor_mut_comparison tmc
INNER JOIN selected_tumor_pairs stp ON tmc.tumor_pair_barcode = stp.tumor_pair_barcode AND priority = 1
LEFT JOIN analysis.mutation_freq mf1 ON mf1.aliquot_barcode = tmc.tumor_barcode_a 
LEFT JOIN analysis.mutation_freq mf2 ON mf2.aliquot_barcode = tmc.tumor_barcode_b 