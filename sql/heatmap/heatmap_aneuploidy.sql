SELECT
	tumor_pair_barcode,
	tumor_barcode_a,
	tumor_barcode_b,
	ss.case_barcode,
	idh_codel_subtype,
	a1.aneuploidy AS aneuploidy_a,
	a2.aneuploidy AS aneuploidy_b,
	t1.aneuploidy_amp_score::integer AS aneuploidy_amp_score_a,
	t2.aneuploidy_amp_score::integer AS aneuploidy_amp_score_b,
	t1.aneuploidy_del_score::integer AS aneuploidy_del_score_a,
	t2.aneuploidy_del_score::integer AS aneuploidy_del_score_b,
	t1.aneuploidy_score::integer AS aneuploidy_score_a,
	t2.aneuploidy_score::integer AS aneuploidy_score_b,
	(CASE WHEN b1.cnv_exclusion <> 'allow' OR b2.cnv_exclusion <> 'allow' THEN 1 ELSE 0 END) qc_fail
FROM analysis.silver_set ss
LEFT JOIN analysis.aneuploidy a1 ON a1.aliquot_barcode = ss.tumor_barcode_a
LEFT JOIN analysis.aneuploidy a2 ON a2.aliquot_barcode = ss.tumor_barcode_b
LEFT JOIN analysis.taylor_aneuploidy t1 ON t1.aliquot_barcode = ss.tumor_barcode_a
LEFT JOIN analysis.taylor_aneuploidy t2 ON t2.aliquot_barcode = ss.tumor_barcode_b
LEFT JOIN clinical.subtypes su ON su.case_barcode = ss.case_barcode
LEFT JOIN analysis.blocklist b1 ON b1.aliquot_barcode = ss.tumor_barcode_a
LEFT JOIN analysis.blocklist b2 ON b2.aliquot_barcode = ss.tumor_barcode_b