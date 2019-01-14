/*
For each `aliquot_barcode` with segmentation data in the `gatk_seg` table:
- Quantify the sum of all segment sizes (this should roughly add up to the total genome size, around 3b)
- Quantify the sum of all *heterozygous* segment sizes
- Quantify the proportion of heterozygous segments sizes vs non-heterozygous segments
*/
WITH
segquant AS
(
	SELECT
		aliquot_barcode,
		sum(upper(pos) - lower(pos) -1) AS seg_size,
		sum(CASE WHEN log2_copy_ratio > log(2.0, 1.0/2.0) AND log2_copy_ratio < log(2.0, 3.0/2.0) THEN (upper(pos) - lower(pos) -1) ELSE 0 END) AS het_size
	FROM analysis.gatk_seg
	GROUP BY aliquot_barcode
)
SELECT *, round(het_size::decimal/seg_size,4) AS prop_het FROM segquant ORDER BY 4