SELECT tmc.case_barcode, tumor_barcode_a, tumor_barcode_b, aliquot_analysis_type, case_source, idh_codel_subtype, received_alk, hypermutator_status, tq1.length AS p_len, tq2.length AS r_len, tqn.length AS n_len, tq1.length / tq2.length AS pr_len_ratio, 
	a1.prop_aneuploidy AS aneuploidy_a,
	a2.prop_aneuploidy AS aneuploidy_b,
	a1.aneuploidy_amp_score::integer AS aneuploidy_amp_score_a,
	a2.aneuploidy_amp_score::integer AS aneuploidy_amp_score_b,
	a1.aneuploidy_del_score::integer AS aneuploidy_del_score_a,
	a2.aneuploidy_del_score::integer AS aneuploidy_del_score_b,
	a1.aneuploidy_score::integer AS aneuploidy_score_a,
	a2.aneuploidy_score::integer AS aneuploidy_score_b
FROM analysis.tumor_mut_comparison_anno tmc
LEFT JOIN analysis.gatk_aneuploidy a1 ON a1.aliquot_barcode = tmc.tumor_barcode_a
LEFT JOIN analysis.gatk_aneuploidy a2 ON a2.aliquot_barcode = tmc.tumor_barcode_b
LEFT JOIN analysis.pairs pa1 ON pa1.tumor_barcode = tmc.tumor_barcode_a
LEFT JOIN biospecimen.aliquots al1 ON al1.aliquot_barcode = tmc.tumor_barcode_a
LEFT JOIN biospecimen.samples sa1 ON sa1.sample_barcode = al1.sample_barcode
LEFT JOIN clinical.cases ca1 ON ca1.case_barcode = sa1.case_barcode
LEFT JOIN analysis.telseq tqn ON tqn.aliquot_barcode = pa1.normal_barcode
LEFT JOIN analysis.telseq tq1 ON tq1.aliquot_barcode = tmc.tumor_barcode_a
LEFT JOIN analysis.telseq tq2 ON tq2.aliquot_barcode = tmc.tumor_barcode_b