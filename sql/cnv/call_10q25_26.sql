/*
Call CNV for 10q25-26 region
*/
WITH
selected_regions AS
(
	SELECT '10q25-26' AS region, * FROM ref.cytobands WHERE chrom = 10 AND substring(cytoband from 1 for 3) IN ('q25','q26')
),
gene_seg_intersect AS
(
    SELECT aliquot_barcode, gs.chrom, (upper(t0.pos * gs.pos) - lower(t0.pos * gs.pos) -1) AS w, 2^log2_copy_ratio::decimal As cr
    FROM variants.gatk_seg gs
    INNER JOIN selected_regions t0 ON t0.chrom = gs.chrom AND t0.pos && gs.pos
),
gene_sample_call AS
(
    SELECT aliquot_barcode, region, 
		sum(w * cr) / sum(w) AS wcr
    FROM gene_seg_intersect
    GROUP BY aliquot_barcode, region
),
seg_stats_optimized AS
(
	SELECT
		gs.aliquot_barcode,
		LEAST(0.9, neu_fwmean - 2 * neu_fwsd) AS del_thres,
		GREATEST(1.1, neu_fwmean + 2 * neu_fwsd) AS amp_thres,
		(CASE
		 WHEN max_loss_arm_wmean < 0.9 AND max_loss_arm_n >= 3 THEN GREATEST(0,max_loss_arm_wmean - 2 * max_loss_arm_wsd)
		 WHEN del_fwmean < 0.9 AND del_n >= 3 THEN GREATEST(0,del_fwmean - 2 * del_fwsd)
		 ELSE NULL
		END) AS hldel_thres,
		(CASE
		 WHEN max_gain_arm_wmean > 1.1 AND max_gain_arm_n >= 3 THEN max_gain_arm_wmean + 2 * max_gain_arm_wsd
		 WHEN amp_fwmean > 1.1 AND amp_n >= 3 THEN amp_fwmean + 2 * amp_fwsd
		 ELSE NULL
		END) AS hlamp_thres
	FROM analysis.gatk_seg_stats gs
	LEFT JOIN analysis.gatk_aneuploidy gsa ON gsa.aliquot_barcode = gs.aliquot_barcode
),
gene_cp AS
(
	SELECT ts.aliquot_barcode, region, ts.chrom, (upper(t0.pos * ts.pos) - lower(t0.pos * ts.pos) -1) AS w, cellular_prevalence As cp
	FROM variants.titan_seg ts
	INNER JOIN selected_regions t0 ON t0.chrom = ts.chrom AND t0.pos && ts.pos
),
gene_cp_agg AS
(
	SELECT aliquot_barcode, region, 
		COALESCE(sum(w * cp) / NULLIF(sum(w),0),NULL) AS wcp
    FROM gene_cp
    GROUP BY 1, 2
)
SELECT
	gc.aliquot_barcode,
	gc.region,
	(CASE
	 WHEN gc.wcr >= del_thres AND gc.wcr <= amp_thres THEN 0
	 WHEN gc.wcr < hldel_thres THEN -2
	 WHEN gc.wcr < del_thres THEN -1
	 WHEN gc.wcr > hlamp_thres THEN 2
	 WHEN gc.wcr > amp_thres THEN 1
	 ELSE NULL
	 END) hlvl_call,
	gc.wcr,
	wcp AS cellular_prevalence
FROM gene_sample_call gc
LEFT JOIN seg_stats_optimized ss ON ss.aliquot_barcode = gc.aliquot_barcode
LEFT JOIN gene_cp_agg cp ON cp.aliquot_barcode = gc.aliquot_barcode AND cp.region = gc.region
ORDER BY 3