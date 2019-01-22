/*
- For each aliquot in the gatk_seg table, calculate sum of all segment sizes (this should roughly add up to the total genome size, around 3b)
- Define a heterozygous segment as as segment that has a CallCopyRatioSegments value of zero (hetozygous)
- Quantify the proportion of heterozygous segments sizes vs non-heterozygous segments
*/
WITH
cnv AS
(
	SELECT
		aliquot_barcode,
		sum(upper(pos) - lower(pos) -1) AS seg_size,
		sum(CASE WHEN gs.call = '0' THEN (upper(pos) - lower(pos) -1) ELSE 0 END) AS het_size
	FROM analysis.gatk_seg gs
	WHERE chrom NOT IN ('X','Y')
	GROUP BY 1
)
SELECT
	aliquot_barcode,
	round(1.0 - het_size::decimal/seg_size,4) AS aneuploidy
FROM cnv
ORDER BY 2