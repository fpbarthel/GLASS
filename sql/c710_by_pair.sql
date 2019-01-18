SELECT tumor_pair_barcode, c1.c710 AS pri_c710, c2.c710 AS rec_c710
FROM analysis.silver_set ss
LEFT JOIN analysis.c710 c1 ON c1.aliquot_barcode = ss.tumor_barcode_a
LEFT JOIN analysis.c710 c2 ON c2.aliquot_barcode = ss.tumor_barcode_b