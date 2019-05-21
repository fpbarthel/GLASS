WITH selected_samples AS
(
	SELECT tumor_barcode_a AS aliquot_barcode, pair_barcode FROM analysis.diamond_set ds
	INNER JOIN analysis.pairs pa ON pa.tumor_barcode = ds.tumor_barcode_a
	UNION
	SELECT tumor_barcode_b AS aliquot_barcode, pair_barcode FROM analysis.diamond_set ds
	INNER JOIN analysis.pairs pa ON pa.tumor_barcode = ds.tumor_barcode_b
)
SELECT ss.aliquot_barcode, cellularity, purity
FROM selected_samples ss
INNER JOIN variants.titan_params tp ON tp.pair_barcode = ss.pair_barcode
INNER JOIN variants.seqz_params sp ON sp.pair_barcode = ss.pair_barcode