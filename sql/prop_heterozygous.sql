/*
For each `pair_barcode` with segmentation data in the `titan_seg` table:
- Quantify the number of segments
- Quantify the sum of all segment sizes (this should roughly add up to the total genome size)
- Quantify the sum of all *heterozygous* segment sizes
- Quantify the proportion of heterozygous segments sizes vs non-heterozygous segments
*/
WITH
segquant AS
(
	SELECT
		pair_barcode,
		COUNT(*) AS num_seg,
		sum(upper(pos) - lower(pos)) AS seg_size,
		sum(CASE WHEN corrected_call IN ('HET','NEUT') THEN upper(pos) - lower(pos) ELSE 0 END) AS het_size
	FROM analysis.titan_seg
	GROUP BY pair_barcode
)
SELECT *, round(het_size::decimal/seg_size,4) AS prop_het FROM segquant ORDER BY 4