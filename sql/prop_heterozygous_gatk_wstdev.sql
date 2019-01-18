/*
---
PostgreSQL implementation of the GATK CallCopyRatioSegments caller.
Method described here: https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/copynumber/caller/SimpleCopyRatioCaller.java#L18-L26
Formula for weighted stdev taken from Stackoverflow: https://stackoverflow.com/questions/31021895/is-there-a-way-to-use-postgresql-function-to-calculate-weighted-average-in-rails
Also here: sqrt( (sum(w*v^2) - (sum(w*v)^2)/sum(w)) /(sum(w)-1))
Where w = weights and v = value
---
1.) Use the non-log2 mean copy ratio to determine copy-neutral segments,
2.) Weight segments by length for determining the mean and standard deviation of the non-log2 copy ratio in copy-neutral segments,
3.) Filter outlier copy-neutral segments by non-log2 copy ratio z-score,
4.) Use the filtered copy-neutral segments to determine a length-weighted mean and standard deviation,
5.) Call segments using z-score based on this mean and standard deviation.
--
*/
-- 

WITH
unfiltered_seg_wmean_wsd AS
(
	SELECT
		aliquot_barcode,
		COUNT(*) AS num_seg,
		(sum((upper(pos) - lower(pos) -1) * 2^log2_copy_ratio::decimal) / sum(upper(pos) - lower(pos) -1) ) AS wmean,
		(sqrt((sum((upper(pos) - lower(pos) -1) * (2^log2_copy_ratio::decimal)^2) - (sum((upper(pos) - lower(pos) -1) * 2^log2_copy_ratio::decimal)^2) / sum(upper(pos) - lower(pos) -1)) / (sum(upper(pos) - lower(pos) -1) -1)))::decimal AS wsd
	FROM analysis.gatk_seg
	WHERE 2^log2_copy_ratio >= 0.9 AND 2^log2_copy_ratio <= 1.1
	GROUP BY 1
),
filtered_seg_wmean_wsd AS
(
	SELECT
		gs.aliquot_barcode,
		COUNT(*) AS num_seg,
		(sum((upper(pos) - lower(pos) -1) * 2^log2_copy_ratio::decimal) / sum(upper(pos) - lower(pos) -1))::decimal  AS fwmean,
		(sqrt((sum((upper(pos) - lower(pos) -1)*(2^log2_copy_ratio::decimal)^2) -(sum((upper(pos) - lower(pos) -1)*2^log2_copy_ratio::decimal)^2)/sum(upper(pos) - lower(pos) -1))/(sum(upper(pos) - lower(pos) -1) -1)))::decimal AS fwsd
	FROM analysis.gatk_seg gs
	INNER JOIN unfiltered_seg_wmean_wsd us ON us.aliquot_barcode = gs.aliquot_barcode
	WHERE 2^log2_copy_ratio >= -2 * wsd + wmean AND 2^log2_copy_ratio <= 2 * wsd + wmean
	GROUP BY 1
),
call_cnv AS
(
	SELECT
		gs.aliquot_barcode,
		gs.chrom,
		gs.pos,
		gs.log2_copy_ratio,
		gs.call,
		fwmean,
		fwsd,
		2^log2_copy_ratio,
		(CASE
		 WHEN 2^log2_copy_ratio >= 0.9 AND 2^log2_copy_ratio <= 1.1 THEN '0'
		 WHEN (2^log2_copy_ratio - fwmean) < -2.0 * fwsd THEN '-'
		 WHEN (2^log2_copy_ratio - fwmean) > 2.0 * fwsd THEN '+'
		 ELSE '0'
		 END) floris_call
	FROM analysis.gatk_seg gs
	INNER JOIN filtered_seg_wmean_wsd fis ON fis.aliquot_barcode = gs.aliquot_barcode
)
SELECT *, call <> floris_call FROM call_cnv WHERE aliquot_barcode = 'GLSS-CU-R017-R1-01D-WXS-C0CBCW'
