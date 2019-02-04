SELECT
	tmc.tumor_pair_barcode,
	tmc.case_barcode,
	idh_codel_subtype,
	received_alkylating_agent,
	hypermutator_status,
	0 AS time_birth,
	ca.case_age_diagnosis_years AS time_initial,
	ROUND(ca.case_age_diagnosis_years + (tmc.surgical_interval_mo / 12.0),2) AS time_recurrence,
	0 AS mf_birth,
	mf1.coverage_adj_mut_freq AS mf_initial,
	mf2.coverage_adj_mut_freq AS mf_recurrence
	/*tmc.count_a,
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
	ROUND(intersection_ab::decimal / LEAST(mf1.cumulative_coverage, mf2.cumulative_coverage) * 1e6, 4) AS mf_shared,*/
FROM analysis.tumor_mut_comparison tmc
INNER JOIN analysis.silver_set stp ON tmc.tumor_pair_barcode = stp.tumor_pair_barcode
LEFT JOIN clinical.clinical_by_tumor_pair ctp ON ctp.tumor_pair_barcode = stp.tumor_pair_barcode
LEFT JOIN analysis.mutation_freq mf1 ON mf1.aliquot_barcode = tmc.tumor_barcode_a 
LEFT JOIN analysis.mutation_freq mf2 ON mf2.aliquot_barcode = tmc.tumor_barcode_b 
LEFT JOIN clinical.subtypes su ON su.case_barcode = stp.case_barcode
LEFT JOIN clinical.cases ca ON ca.case_barcode = stp.case_barcode