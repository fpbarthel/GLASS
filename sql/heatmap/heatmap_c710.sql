SELECT
	tumor_pair_barcode,
	ss.case_barcode,
	idh_codel_subtype,
	(CASE
	 WHEN c1.c710 IS FALSE AND c2.c710 IS FALSE THEN 'WT'
	 WHEN c1.c710 IS TRUE AND c2.c710 IS TRUE THEN 'S'
	 WHEN c1.c710 IS TRUE AND c2.c710 IS FALSE THEN 'P'
	 WHEN c1.c710 IS FALSE AND c2.c710 IS TRUE THEN 'R'
	 ELSE NULL
	END) c710_status
FROM analysis.gold_set ss
LEFT JOIN analysis.gatk_c710_status c1 ON c1.aliquot_barcode = ss.tumor_barcode_a
LEFT JOIN analysis.gatk_c710_status c2 ON c2.aliquot_barcode = ss.tumor_barcode_b
LEFT JOIN clinical.subtypes su ON su.case_barcode = ss.case_barcode