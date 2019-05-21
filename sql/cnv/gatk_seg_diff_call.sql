/*
For each `tumor_pair_barcode`, consisting of `tumor_barcode_a` (a) and `tumor_barcode_b` (b) in the `tumor_pairs` table:
- Intersect all segments from (a) and (b) --> `pos_intersect`
- Determine change in `minor_cn` between (a) and (b) as (b-a)
- Determine change in `major_cn` between (a) and (b) as (b-a)
	* For the two items above, a value of 0 indicates no change, positive numbers indicate (b > a) and negative numbers (a > b)
- Determine the difference in log2(tumor/normal) copy number between (a) and (b) as (b-a)
*/
/*
For each `tumor_pair_barcode`, consisting of `tumor_barcode_a` (a) and `tumor_barcode_b` (b) in the `tumor_pairs` table:
- Intersect all segments from (a) and (b) --> `pos_intersect`
- Determine change in `minor_cn` between (a) and (b) as (b-a)
- Determine change in `major_cn` between (a) and (b) as (b-a)
	* For the two items above, a value of 0 indicates no change, positive numbers indicate (b > a) and negative numbers (a > b)
- Determine the difference in log2(tumor/normal) copy number between (a) and (b) as (b-a)
*/
WITH
seg_diff AS
(
	SELECT
		pa.tumor_pair_barcode, 
		pa.case_barcode, 
		pa.tumor_barcode_a, 
		pa.tumor_barcode_b,
		s1.chrom,
		s1.pos * s2.pos AS pos, 
		s2.log2_copy_ratio - s1.log2_copy_ratio AS delta_log2_copy_ratio
	FROM analysis.tumor_pairs pa
	INNER JOIN variants.gatk_seg s1 ON s1.aliquot_barcode = pa.tumor_barcode_a
	INNER JOIN variants.gatk_seg s2 ON s2.aliquot_barcode = pa.tumor_barcode_b AND s1.chrom = s2.chrom AND s1.pos && s2.pos
	LIMIT 1000
),
unfiltered_seg_wmean_wsd AS
(
	SELECT
		tumor_pair_barcode, 
		case_barcode, 
		tumor_barcode_a, 
		tumor_barcode_b,
		COUNT(*) AS num_seg,
		(sum((upper(pos) - lower(pos) -1) * 2^delta_log2_copy_ratio) / sum(upper(pos) - lower(pos) -1) )::decimal AS wmean,
		(sqrt((sum((upper(pos) - lower(pos) -1) * (2^delta_log2_copy_ratio)^2) - (sum((upper(pos) - lower(pos) -1) * 2^delta_log2_copy_ratio)^2) / sum(upper(pos) - lower(pos) -1)) / (sum(upper(pos) - lower(pos) -1) -1)))::decimal AS wsd
	FROM seg_diff gs
	WHERE
		2^delta_log2_copy_ratio >= 0.9 AND
		2^delta_log2_copy_ratio <= 1.1
	GROUP BY 1,2,3,4
),
filtered_seg_wmean_wsd AS
(
	SELECT
		gs.tumor_pair_barcode, 
		gs.case_barcode, 
		gs.tumor_barcode_a, 
		gs.tumor_barcode_b,
		num_seg,
		wmean,
		wsd,
		COUNT(*) AS fnum_seg,
		(sum((upper(pos) - lower(pos) -1) * 2^delta_log2_copy_ratio) / sum(upper(pos) - lower(pos) -1))::decimal  AS fwmean,
		(sqrt((sum((upper(pos) - lower(pos) -1)*(2^delta_log2_copy_ratio)^2) -(sum((upper(pos) - lower(pos) -1)*2^delta_log2_copy_ratio)^2)/sum(upper(pos) - lower(pos) -1))/(sum(upper(pos) - lower(pos) -1) -1)))::decimal AS fwsd
	FROM seg_diff gs
	INNER JOIN unfiltered_seg_wmean_wsd us ON us.tumor_barcode_a = gs.tumor_barcode_a AND us.tumor_barcode_b = gs.tumor_barcode_b
	WHERE
		2^delta_log2_copy_ratio >= 0.9 AND
		2^delta_log2_copy_ratio <= 1.1 AND
		(2^delta_log2_copy_ratio - wmean) > -2.0 * wsd AND
		(2^delta_log2_copy_ratio - wmean) < 2.0 * wsd 
	GROUP BY 1,2,3,4,5,6,7
),
call_cnv AS
(
	SELECT
		gs.tumor_pair_barcode, 
		gs.case_barcode, 
		gs.tumor_barcode_a, 
		gs.tumor_barcode_b,
		gs.chrom,
		gs.pos,
		gs.delta_log2_copy_ratio,
		(CASE
		 WHEN 2^delta_log2_copy_ratio >= 0.9 AND 2^delta_log2_copy_ratio <= 1.1 THEN 0
		 WHEN (2^delta_log2_copy_ratio - fwmean) < -2.0 * fwsd THEN -1
		 WHEN (2^delta_log2_copy_ratio - fwmean) > 2.0 * fwsd THEN 1
		 ELSE 0
		 END) cnv_call
	FROM seg_diff gs
	INNER JOIN filtered_seg_wmean_wsd fis ON fis.tumor_barcode_a = gs.tumor_barcode_a AND fis.tumor_barcode_b = gs.tumor_barcode_b
)
SELECT * FROM call_cnv