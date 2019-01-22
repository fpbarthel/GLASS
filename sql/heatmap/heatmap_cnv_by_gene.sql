/*
Fractionated CNV proportions per gene (and subtype)
*/
WITH
selected_tumor_pairs AS
(
	SELECT * FROM analysis.silver_set
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
	LEFT JOIN analysis.taylor_aneuploidy gsa ON gsa.aliquot_barcode = gs.aliquot_barcode
),
cnv_by_pair_gene AS
(
	SELECT
		sgs.tumor_pair_barcode,
		sgs.case_barcode,
		sgs.gene_symbol,
		pathway,
		ROUND(c1.wcr,6) AS cr_a,
		ROUND(c2.wcr,6) AS cr_b,
		(CASE
		 WHEN cnv_direction = -1 AND c1.hlvl_call = -2 THEN true
		 WHEN cnv_direction = 1 AND c1.hlvl_call = 2 THEN true
		 ELSE false
		 END) selected_call_a,
		(CASE
		 WHEN cnv_direction = -1 AND c2.hlvl_call = -2 THEN true
		 WHEN cnv_direction = 1 AND c2.hlvl_call = 2 THEN true
		 ELSE false
		 END) selected_call_b,
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
	LEFT JOIN analysis.cnv_by_gene_gatk c1 ON c1.aliquot_barcode = sgs.tumor_barcode_a AND c1.gene_symbol = sgs.gene_symbol
	LEFT JOIN analysis.cnv_by_gene_gatk c2 ON c2.aliquot_barcode = sgs.tumor_barcode_b AND c2.gene_symbol = sgs.gene_symbol
	LEFT JOIN seg_stats_optimized ss1 ON ss1.aliquot_barcode = sgs.tumor_barcode_a
	LEFT JOIN seg_stats_optimized ss2 ON ss2.aliquot_barcode = sgs.tumor_barcode_b
),
cnv_by_gene AS
(
	SELECT
		gene_symbol,
		SUM(CASE WHEN selected_call_a THEN 1 ELSE 0 END)::integer		AS count_a,											-- count of mutations in A
		SUM(CASE WHEN selected_call_b THEN 1 ELSE 0 END)::integer		AS count_b,											-- count of mutations in B
		SUM(CASE WHEN selected_call_a AND selected_call_b THEN 1 ELSE 0 END)::integer 		AS shared, 						-- if both A and B are true then a variant is shared
		SUM(CASE WHEN selected_call_a AND NOT selected_call_b THEN 1 ELSE 0 END)::integer 	AS private_a,					-- if A is true and B is not then it is unique to A
		SUM(CASE WHEN selected_call_b AND NOT selected_call_a THEN 1 ELSE 0 END)::integer 	AS private_b,					-- if B is true and A is not then it is unique to B,
		SUM(CASE WHEN selected_call_b OR selected_call_a THEN 1 ELSE 0 END)::integer 		AS total,						-- if a is true OR B is true, sums to total variants between A and B
		ROUND(SUM(CASE WHEN selected_call_a AND selected_call_b THEN 1 ELSE 0 END) / 
			  SUM(CASE WHEN selected_call_b OR selected_call_a THEN 1 ELSE 0 END)::decimal,2) AS prop_shared,		-- divides `shared` by total to get the proportion shared
		ROUND(SUM(CASE WHEN selected_call_a AND NOT selected_call_b THEN 1 ELSE 0 END) /
			  SUM(CASE WHEN selected_call_b OR selected_call_a THEN 1 ELSE 0 END)::decimal,2) AS prop_private_a,	-- divides `private_a` by total to get the proportion private to A
		ROUND(SUM(CASE WHEN selected_call_b AND NOT selected_call_a THEN 1 ELSE 0 END) /
			  SUM(CASE WHEN selected_call_b OR selected_call_a THEN 1 ELSE 0 END)::decimal,2) AS prop_private_b		-- divides `private_b` by total to get the proportion private to B
	FROM cnv_by_pair_gene
	GROUP BY 1	-- aggregate across genes
	ORDER BY 7 DESC -- order by column 7 (= total)
),
cnv_by_gene_subtype AS
(
	SELECT
		idh_codel_subtype,
		gene_symbol,	
		SUM(CASE WHEN selected_call_a THEN 1 ELSE 0 END)::integer		AS count_a,											-- count of mutations in A
		SUM(CASE WHEN selected_call_b THEN 1 ELSE 0 END)::integer		AS count_b,											-- count of mutations in B
		SUM(CASE WHEN selected_call_a AND selected_call_b THEN 1 ELSE 0 END)::integer 		AS shared, 						-- if both A and B are true then a variant is shared
		SUM(CASE WHEN selected_call_a AND NOT selected_call_b THEN 1 ELSE 0 END)::integer 	AS private_a,					-- if A is true and B is not then it is unique to A
		SUM(CASE WHEN selected_call_b AND NOT selected_call_a THEN 1 ELSE 0 END)::integer 	AS private_b,					-- if B is true and A is not then it is unique to B,
		SUM(CASE WHEN selected_call_b OR selected_call_a THEN 1 ELSE 0 END)::integer 		AS total,						-- if a is true OR B is true, sums to total variants between A and B
		ROUND(SUM(CASE WHEN selected_call_a AND selected_call_b THEN 1 ELSE 0 END) / 
			  NULLIF(SUM(CASE WHEN selected_call_b OR selected_call_a THEN 1 ELSE 0 END),0)::decimal,2) AS prop_shared,		-- divides `shared` by total to get the proportion shared
		ROUND(SUM(CASE WHEN selected_call_a AND NOT selected_call_b THEN 1 ELSE 0 END) /
			  NULLIF(SUM(CASE WHEN selected_call_b OR selected_call_a THEN 1 ELSE 0 END),0)::decimal,2) AS prop_private_a,	-- divides `private_a` by total to get the proportion private to A
		ROUND(SUM(CASE WHEN selected_call_b AND NOT selected_call_a THEN 1 ELSE 0 END) /
			  NULLIF(SUM(CASE WHEN selected_call_b OR selected_call_a THEN 1 ELSE 0 END),0)::decimal,2) AS prop_private_b		-- divides `private_b` by total to get the proportion private to B
	FROM cnv_by_pair_gene vg
	LEFT JOIN clinical.subtypes st ON st.case_barcode = vg.case_barcode
	GROUP BY 1,2	-- aggregate across genes
	ORDER BY 7 DESC -- order by column 7 (= total)
)
(SELECT 'all' AS idh_codel_subtype, vg.*
FROM cnv_by_gene vg
		 
UNION
		 
SELECT vs.*
FROM cnv_by_gene_subtype vs)
ORDER BY 1 ASC, 8 DESC