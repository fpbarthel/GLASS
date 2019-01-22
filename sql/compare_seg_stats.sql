/*
For each `tumor_pair_barcode` in the `tumor_pairs` table:
- Compare number of segments between (a) and (b)
- Compare proportion of genome that is heterozyous between (a) and (b)
*/
SELECT
	tumor_pair_barcode,
	tumor_barcode_a,
	tumor_barcode_b,
	s1.num_seg AS num_seg_a,
	s2.num_seg AS num_seg_b,
	s1.prop_het AS prop_het_a,
	s2.prop_het AS prop_het_b,
	s2.num_seg - s1.num_seg AS delta_num_seg,
	s2.prop_het - s1.prop_het AS delta_prop_het
FROM analysis.tumor_pairs pa
LEFT JOIN analysis.pairs p1 ON p1.tumor_barcode = pa.tumor_barcode_a
LEFT JOIN analysis.pairs p2 ON p2.tumor_barcode = pa.tumor_barcode_b
LEFT JOIN analysis.titan_seg_prop_het s1 ON s1.pair_barcode = p1.pair_barcode
LEFT JOIN analysis.titan_seg_prop_het s2 ON s2.pair_barcode = p2.pair_barcode
ORDER BY 9