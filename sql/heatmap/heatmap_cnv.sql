WITH
selected_tumor_pairs AS
(
	SELECT ss.tumor_pair_barcode, ss.tumor_barcode_a, ss.tumor_barcode_b, ss.case_barcode, (CASE WHEN gs.tumor_pair_barcode IS NULL THEN 'Silver set' ELSE 'Gold set' END) AS gold_set
	FROM analysis.gold_set ss
	LEFT JOIN analysis.gold_set gs ON gs.tumor_pair_barcode = ss.tumor_pair_barcode
),
selected_genes AS
(
	SELECT dr.gene_symbol,chrom,pos,cnv_direction,pathway
	FROM ref.driver_genes dr
	LEFT JOIN ref.genes ge ON ge.gene_symbol = dr.gene_symbol
	WHERE has_cnv IS TRUE
),
selected_genes_samples AS
(
	SELECT * FROM selected_tumor_pairs, selected_genes
),
seg_stats_optimized AS
(
	SELECT
		gs.aliquot_barcode,
		LEAST(0.9, neu_fwmean - 2 * neu_fwsd) AS del_thres,
		GREATEST(1.1, neu_fwmean + 2 * neu_fwsd) AS amp_thres,
		
		-- high level deletion threshold
		(CASE
		 WHEN max_loss_arm_wmean < 0.9 AND max_loss_arm_n >= 3 THEN GREATEST(0,max_loss_arm_wmean - 2 * max_loss_arm_wsd)
		 WHEN del_fwmean < 0.9 AND del_n >= 3 THEN GREATEST(0,del_fwmean - 2 * del_fwsd)
		 ELSE NULL
		END) AS hldel_thres,
		
		-- high level deletion filtered/weighted mean
		(CASE
		 WHEN max_loss_arm_wmean < 0.9 AND max_loss_arm_n >= 3 THEN max_loss_arm_wmean
		 WHEN del_fwmean < 0.9 AND del_n >= 3 THEN del_fwmean
		 ELSE NULL
		END) AS hldel_fwmean,
	
		-- high level deletion filtered/weighted standard deviation
		(CASE
		 WHEN max_loss_arm_wmean < 0.9 AND max_loss_arm_n >= 3 THEN max_loss_arm_wsd
		 WHEN del_fwmean < 0.9 AND del_n >= 3 THEN del_fwsd
		 ELSE NULL
		END) AS hldel_fwsd,
		
		-- high level amplification threshold
		(CASE
		 WHEN max_gain_arm_wmean > 1.1 AND max_gain_arm_n >= 3 THEN max_gain_arm_wmean + 2 * max_gain_arm_wsd
		 WHEN amp_fwmean > 1.1 AND amp_n >= 3 THEN amp_fwmean + 2 * amp_fwsd
		 ELSE NULL
		END) AS hlamp_thres,
	
		-- high level amplification filtered/weighted mean
		(CASE
		 WHEN max_gain_arm_wmean > 1.1 AND max_gain_arm_n >= 3 THEN max_gain_arm_wmean
		 WHEN amp_fwmean > 1.1 AND amp_n >= 3 THEN amp_fwmean
		 ELSE NULL
		END) AS hlamp_fwmean,
	
		-- high level amplification filtered/weighted standard deviation
		(CASE
		 WHEN max_gain_arm_wmean > 1.1 AND max_gain_arm_n >= 3 THEN max_gain_arm_wsd
		 WHEN amp_fwmean > 1.1 AND amp_n >= 3 THEN amp_fwsd
		 ELSE NULL
		END) AS hlamp_fwsd
	FROM analysis.gatk_seg_stats gs
	LEFT JOIN analysis.gatk_aneuploidy gsa ON gsa.aliquot_barcode = gs.aliquot_barcode
),
cnv_by_pair_gene AS
(
	SELECT
		sgs.tumor_pair_barcode,
		sgs.case_barcode,
		sgs.gene_symbol,
		gold_set,
		pathway,
		ROUND(c1.wcr,6) AS cr_a,
		ROUND(c2.wcr,6) AS cr_b,
		(CASE
		 WHEN cnv_direction = -1 AND (c1.hlvl_call = -2 OR c2.hlvl_call = -2) THEN 'HLDEL'
		 WHEN cnv_direction = 1 AND (c1.hlvl_call = 2 OR c2.hlvl_call = 2) THEN 'HLAMP'
		 ELSE NULL
		 END) cnv_state,
		(CASE
		 WHEN cnv_direction = -1 AND (c1.hlvl_call = -2 OR c2.hlvl_call = -2) AND c1.hlvl_call < c2.hlvl_call THEN 'P'
		 WHEN cnv_direction = -1 AND (c1.hlvl_call = -2 OR c2.hlvl_call = -2) AND c1.hlvl_call = c2.hlvl_call THEN 'S'
		 WHEN cnv_direction = -1 AND (c1.hlvl_call = -2 OR c2.hlvl_call = -2) AND c1.hlvl_call > c2.hlvl_call THEN 'R'
		 WHEN cnv_direction = 1 AND (c1.hlvl_call = 2 OR c2.hlvl_call = 2) AND c1.hlvl_call > c2.hlvl_call THEN 'P'
		 WHEN cnv_direction = 1 AND (c1.hlvl_call = 2 OR c2.hlvl_call = 2) AND c1.hlvl_call = c2.hlvl_call THEN 'S'
		 WHEN cnv_direction = 1 AND (c1.hlvl_call = 2 OR c2.hlvl_call = 2) AND c1.hlvl_call < c2.hlvl_call THEN 'R'
		 ELSE NULL
		 END) cnv_change,
		--ss1.hldel_fwmean,
		--ss1.hlamp_fwmean,
		--ss1.hldel_fwsd,
		--ss1.hlamp_fwsd,
		ROUND((c1.wcr - ss1.hldel_fwmean)/ss1.hldel_fwsd,6) AS hldel_zs_a,
		ROUND((c2.wcr - ss2.hldel_fwmean)/ss2.hldel_fwsd,6) AS hldel_zs_b,
		ROUND((c1.wcr - ss1.hlamp_fwmean)/ss1.hlamp_fwsd,6) AS hlamp_zs_a,
		ROUND((c2.wcr - ss2.hlamp_fwmean)/ss2.hlamp_fwsd,6) AS hlamp_zs_b,
		c1.hlvl_call AS cn_a,
		c2.hlvl_call AS cn_b,
		ss1.del_thres AS del_thres_a,
		ss1.amp_thres AS amp_thres_a,
		ss1.hldel_thres AS hldel_thres_a,
		ss1.hlamp_thres AS hlamp_thres_a,
		ss2.del_thres AS del_thres_b,
		ss2.amp_thres AS amp_thres_b,
		ss2.hldel_thres AS hldel_thres_b,
		ss2.hlamp_thres AS hlamp_thres_b
	FROM selected_genes_samples sgs
	LEFT JOIN analysis.gatk_cnv_by_gene c1 ON c1.aliquot_barcode = sgs.tumor_barcode_a AND c1.gene_symbol = sgs.gene_symbol
	LEFT JOIN analysis.gatk_cnv_by_gene c2 ON c2.aliquot_barcode = sgs.tumor_barcode_b AND c2.gene_symbol = sgs.gene_symbol
	LEFT JOIN seg_stats_optimized ss1 ON ss1.aliquot_barcode = sgs.tumor_barcode_a
	LEFT JOIN seg_stats_optimized ss2 ON ss2.aliquot_barcode = sgs.tumor_barcode_b
)
SELECT cn.*,idh_codel_subtype --cn_a,cn_b,COUNT(*) --case_barcode, gene_symbol, cnv_retention, cnv_class
FROM cnv_by_pair_gene cn --GROUP BY 1,2
--LEFT JOIN analysis.gatk_seg_retention_states gs ON gs.cn_a = cn.cn_a AND gs.cn_b = cn.cn_b
LEFT JOIN clinical.subtypes st ON st.case_barcode = cn.case_barcode