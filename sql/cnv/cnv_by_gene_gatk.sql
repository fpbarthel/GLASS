/*
---
Perform gene level copy number calling using a modified GATK CallCopyRatioSegments / ReCapSeg caller combined with arm level thresholds inspired by the
Taylor aneuploidy paper and GISTIC2 methods
---
1. Define a list of genes on which to call copy number
2. Intersect gene start/stop with segments
3. Group results by gene and sample, and compute weighted average of non-log2 copy number
4. Perform simple calling using -1 and 1 thresholds
---
*/
WITH
selected_genes AS
(
	SELECT dr.gene_symbol,chrom,pos
	FROM ref.driver_genes dr
	LEFT JOIN ref.genes ge ON ge.gene_symbol = dr.gene_symbol
	WHERE has_cnv IS TRUE
),
gene_seg_intersect AS
(
    SELECT aliquot_barcode, gene_symbol, gs.chrom, (upper(t0.pos * gs.pos) - lower(t0.pos * gs.pos) -1) AS w, 2^log2_copy_ratio::decimal As cr
    FROM analysis.gatk_seg gs
    INNER JOIN selected_genes t0 ON t0.chrom = gs.chrom AND t0.pos && gs.pos
),
gene_sample_call AS
(
    SELECT aliquot_barcode, gene_symbol, 
		sum(w * cr) / sum(w) AS wcr
    FROM gene_seg_intersect
    GROUP BY aliquot_barcode, gene_symbol
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
	LEFT JOIN analysis.taylor_aneuploidy gsa ON gsa.aliquot_barcode = gs.aliquot_barcode
)
SELECT
	gc.aliquot_barcode,
	gc.gene_symbol,
	(CASE
	 WHEN gc.wcr >= del_thres AND gc.wcr <= amp_thres THEN 0
	 WHEN gc.wcr < hldel_thres THEN -2
	 WHEN gc.wcr < del_thres THEN -1
	 WHEN gc.wcr > hlamp_thres THEN 2
	 WHEN gc.wcr > amp_thres THEN 1
	 ELSE NULL
	 END) hlvl_call 
FROM gene_sample_call gc
LEFT JOIN seg_stats_optimized ss ON ss.aliquot_barcode = gc.aliquot_barcode
ORDER BY 3