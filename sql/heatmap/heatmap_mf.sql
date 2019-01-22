SELECT
	tmc.tumor_pair_barcode,
	tmc.case_barcode,
	idh_codel_subtype,
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
INNER JOIN analysis.silver_set stp ON tmc.tumor_pair_barcode = stp.tumor_pair_barcode
LEFT JOIN analysis.mutation_freq mf1 ON mf1.aliquot_barcode = tmc.tumor_barcode_a 
LEFT JOIN analysis.mutation_freq mf2 ON mf2.aliquot_barcode = tmc.tumor_barcode_b 
LEFT JOIN clinical.subtypes su ON su.case_barcode = stp.case_barcode