SELECT
	tumor_pair_barcode, ss.case_barcode, idh_codel_subtype, tumor_barcode_a, tumor_barcode_b, 
	s1.most_probable_classification AS evolution_a,
	s2.most_probable_classification AS evolution_b,
	(CASE WHEN s1.most_probable_classification IS NOT NULL AND s2.most_probable_classification IS NOT NULL THEN s1.most_probable_classification || '-' || s2.most_probable_classification END) as evolution_ab
FROM analysis.silver_set ss
LEFT JOIN clinical.subtypes st ON st.case_barcode = ss.case_barcode
LEFT JOIN analysis.subclonalselection s1 ON s1.aliquot_barcode = ss.tumor_barcode_a
LEFT JOIN analysis.subclonalselection s2 ON s2.aliquot_barcode = ss.tumor_barcode_b