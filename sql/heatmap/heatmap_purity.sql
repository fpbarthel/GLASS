SELECT
	ss.case_barcode,
	tp1.purity AS purity_a,
	tp2.purity AS purity_b,
	sp1.cellularity AS seqz_purity_a,
	sp1.cellularity AS seqz_purity_b,
	idh_codel_subtype
FROM analysis.gold_set ss
LEFT JOIN analysis.pairs p1 ON p1.tumor_barcode = ss.tumor_barcode_a
LEFT JOIN variants.titan_params tp1 ON tp1.pair_barcode = p1.pair_barcode
LEFT JOIN analysis.pairs p2 ON p2.tumor_barcode = ss.tumor_barcode_b
LEFT JOIN variants.titan_params tp2 ON tp2.pair_barcode = p2.pair_barcode
LEFT JOIN clinical.subtypes su ON su.case_barcode = ss.case_barcode
LEFT JOIN variants.seqz_params sp1 ON sp1.pair_barcode = p1.pair_barcode
LEFT JOIN variants.seqz_params sp2 ON sp2.pair_barcode = p2.pair_barcode