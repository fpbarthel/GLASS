/*
---
PostgreSQL implementation of the GATK CallCopyRatioSegments caller.
Method described here: https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/copynumber/caller/SimpleCopyRatioCaller.java#L18-L26
See also: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_copynumber_CallCopyRatioSegments.php
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
Parameters:
- copy neutral boundaries: 0.9 - 1.1
- outlier SD boundaries for filtering neutral segments: 2 standard deviations
- SD boundaries beyond which to call CNV amp/loss: 2 standard deviations
-- 
*/

WITH
unfiltered_seg_wmean_wsd AS
(
	SELECT
		aliquot_barcode,
		COUNT(*) AS num_seg,
		(sum((upper(pos) - lower(pos) -1) * 2^log2_copy_ratio) / sum(upper(pos) - lower(pos) -1) )::decimal AS wmean,
		(sqrt((sum((upper(pos) - lower(pos) -1) * (2^log2_copy_ratio)^2) - (sum((upper(pos) - lower(pos) -1) * 2^log2_copy_ratio)^2) / sum(upper(pos) - lower(pos) -1)) / (sum(upper(pos) - lower(pos) -1) -1)))::decimal AS wsd
	FROM analysis.gatk_seg
	WHERE
		2^log2_copy_ratio >= 0.9 AND
		2^log2_copy_ratio <= 1.1
	GROUP BY 1
),
filtered_seg_wmean_wsd AS
(
	SELECT
		gs.aliquot_barcode,
		num_seg,
		wmean,
		wsd,
		COUNT(*) AS fnum_seg,
		(sum((upper(pos) - lower(pos) -1) * 2^log2_copy_ratio) / sum(upper(pos) - lower(pos) -1))::decimal  AS fwmean,
		(sqrt((sum((upper(pos) - lower(pos) -1)*(2^log2_copy_ratio)^2) -(sum((upper(pos) - lower(pos) -1)*2^log2_copy_ratio)^2)/sum(upper(pos) - lower(pos) -1))/(sum(upper(pos) - lower(pos) -1) -1)))::decimal AS fwsd
	FROM analysis.gatk_seg gs
	INNER JOIN unfiltered_seg_wmean_wsd us ON us.aliquot_barcode = gs.aliquot_barcode
	WHERE
		2^log2_copy_ratio >= 0.9 AND
		2^log2_copy_ratio <= 1.1 AND
		(2^log2_copy_ratio - wmean) > -2.0 * wsd AND
		(2^log2_copy_ratio - wmean) < 2.0 * wsd 
	GROUP BY 1,2,3,4
),
call_cnv AS
(
	SELECT
		gs.aliquot_barcode,
		gs.chrom,
		gs.pos,
		gs.log2_copy_ratio,
		--gs.call,
		--fwmean,
		--fwsd,
		--2^log2_copy_ratio,
		(CASE
		 WHEN 2^log2_copy_ratio >= 0.9 AND 2^log2_copy_ratio <= 1.1 THEN 0
		 WHEN (2^log2_copy_ratio - fwmean) < -2.0 * fwsd THEN -1
		 WHEN (2^log2_copy_ratio - fwmean) > 2.0 * fwsd THEN 1
		 ELSE 0
		 END) cnv_call
	FROM analysis.gatk_seg gs
	INNER JOIN filtered_seg_wmean_wsd fis ON fis.aliquot_barcode = gs.aliquot_barcode
)
--SELECT * FROM call_cnv WHERE call <> floris_call
--SELECT * FROM filtered_seg_wmean_wsd WHERE aliquot_barcode = 'GLSS-HF-3081-R2-01D-WXS-4LBL2G'																																	
--SELECT call,floris_call,COUNT(*) FROM call_cnv GROUP BY 1,2
SELECT * FROM call_cnv																																					  