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
---
o Updated to also compute lenght weighted mean and standard deviation across all amplified and deleted segments
o To be used as alternative high level thresholds in case arm level thresholds are not informative
---
Parameters:
- copy neutral boundaries: 0.9 - 1.1
- outlier SD boundaries for filtering neutral segments: 2 standard deviations
- SD boundaries beyond which to call CNV amp/loss: 2 standard deviations
---
This code is used to power two tables:
- The `call_cnv` table described here powers the table `gatk_seg_call`
- The `seg_stats` table described here powers the table `gatk_seg_stats`
---
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
		(CASE
		 WHEN 2^log2_copy_ratio >= 0.9 AND 2^log2_copy_ratio <= 1.1 THEN 0
		 WHEN (2^log2_copy_ratio - fwmean) < -2.0 * fwsd THEN -1
		 WHEN (2^log2_copy_ratio - fwmean) > 2.0 * fwsd THEN 1
		 ELSE 0
		 END) cnv_call
	FROM analysis.gatk_seg gs
	INNER JOIN filtered_seg_wmean_wsd fis ON fis.aliquot_barcode = gs.aliquot_barcode
),
unfiltered_loss_wmean_wsd AS
(
	SELECT
		aliquot_barcode,
		COUNT(*) AS num_seg,
		(sum((upper(pos) - lower(pos) -1) * 2^log2_copy_ratio) / sum(upper(pos) - lower(pos) -1) )::decimal AS wmean,
		(CASE WHEN COUNT(*) > 1 THEN (sqrt((sum((upper(pos) - lower(pos) -1) * (2^log2_copy_ratio)^2) - (sum((upper(pos) - lower(pos) -1) * 2^log2_copy_ratio)^2) / sum(upper(pos) - lower(pos) -1)) / (sum(upper(pos) - lower(pos) -1) -1)))::decimal ELSE NULL END) AS wsd
	FROM call_cnv gs
	WHERE cnv_call = -1
	GROUP BY 1
),
filtered_loss_wmean_wsd AS
(
	SELECT
		gs.aliquot_barcode,
		num_seg,
		wmean,
		wsd,
		COUNT(*) AS fnum_seg,
		(sum((upper(pos) - lower(pos) -1) * 2^log2_copy_ratio) / sum(upper(pos) - lower(pos) -1))::decimal  AS fwmean,
		(CASE WHEN COUNT(*) > 1 THEN (sqrt((sum((upper(pos) - lower(pos) -1)*(2^log2_copy_ratio)^2) -(sum((upper(pos) - lower(pos) -1)*2^log2_copy_ratio)^2)/sum(upper(pos) - lower(pos) -1))/(sum(upper(pos) - lower(pos) -1) -1)))::decimal ELSE NULL END) AS fwsd
	FROM call_cnv gs
	INNER JOIN unfiltered_loss_wmean_wsd us ON us.aliquot_barcode = gs.aliquot_barcode
	WHERE
		cnv_call = -1 AND
		(2^log2_copy_ratio - wmean) > -2.0 * wsd AND
		(2^log2_copy_ratio - wmean) < 2.0 * wsd 
	GROUP BY 1,2,3,4
),
unfiltered_gain_wmean_wsd AS
(
	SELECT
		aliquot_barcode,
		COUNT(*) AS num_seg,
		(sum((upper(pos) - lower(pos) -1) * 2^log2_copy_ratio) / sum(upper(pos) - lower(pos) -1) )::decimal AS wmean,
		(CASE WHEN COUNT(*) > 1 THEN (sqrt((sum((upper(pos) - lower(pos) -1) * (2^log2_copy_ratio)^2) - (sum((upper(pos) - lower(pos) -1) * 2^log2_copy_ratio)^2) / sum(upper(pos) - lower(pos) -1)) / (sum(upper(pos) - lower(pos) -1) -1)))::decimal ELSE NULL END) AS wsd
	FROM call_cnv gs
	WHERE cnv_call = 1
	GROUP BY 1
),
filtered_gain_wmean_wsd AS
(
	SELECT
		gs.aliquot_barcode,
		num_seg,
		wmean,
		wsd,
		COUNT(*) AS fnum_seg,
		(sum((upper(pos) - lower(pos) -1) * 2^log2_copy_ratio) / sum(upper(pos) - lower(pos) -1))::decimal  AS fwmean,
		(CASE WHEN COUNT(*) > 1 THEN (sqrt((sum((upper(pos) - lower(pos) -1)*(2^log2_copy_ratio)^2) -(sum((upper(pos) - lower(pos) -1)*2^log2_copy_ratio)^2)/sum(upper(pos) - lower(pos) -1))/(sum(upper(pos) - lower(pos) -1) -1)))::decimal ELSE NULL END) AS fwsd
	FROM call_cnv gs
	INNER JOIN unfiltered_gain_wmean_wsd us ON us.aliquot_barcode = gs.aliquot_barcode
	WHERE
		cnv_call = 1 AND
		(2^log2_copy_ratio - wmean) > -2.0 * wsd AND
		(2^log2_copy_ratio - wmean) < 2.0 * wsd 
	GROUP BY 1,2,3,4
),
seg_stats AS
(
	SELECT
		al.aliquot_barcode,
		neu.fnum_seg AS neu_n,
		neu.fwmean AS neu_fwmean,
		neu.fwsd AS neu_fwsd,
		amp.fnum_seg AS amp_n,
		amp.fwmean AS amp_fwmean,
		amp.fwsd AS amp_fwsd,
		del.fnum_seg AS del_n,
		del.fwmean AS del_fwmean,
		del.fwsd AS del_fwsd
	FROM biospecimen.aliquots al
	LEFT JOIN filtered_seg_wmean_wsd neu ON al.aliquot_barcode = neu.aliquot_barcode
	LEFT JOIN filtered_loss_wmean_wsd del ON al.aliquot_barcode = del.aliquot_barcode
	LEFT JOIN filtered_gain_wmean_wsd amp ON al.aliquot_barcode = amp.aliquot_barcode
)
--SELECT * FROM call_cnv
SELECT * FROM seg_stats																																					  