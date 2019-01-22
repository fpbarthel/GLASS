/*
Determine CNV driver changes between selected primary and recurrences
*/
WITH
selected_tumor_pairs AS
(
	SELECT
		tumor_pair_barcode,
		ps.case_barcode,
		idh_codel_subtype,
		tumor_barcode_a,
		tumor_barcode_b
	FROM analysis.tumor_pairs ps
	LEFT JOIN analysis.blocklist b1 ON b1.aliquot_barcode = ps.tumor_barcode_a
	LEFT JOIN analysis.blocklist b2 ON b2.aliquot_barcode = ps.tumor_barcode_b
	LEFT JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
	WHERE
		sample_type_b <> 'M1' AND
		b1.coverage_exclusion = 'allow' AND b2.coverage_exclusion = 'allow' 
),
selected_genes AS
(
	SELECT
		dg.gene_symbol,
		chrom,
		pos,
		cnv_direction,
		(CASE WHEN cnv_direction = 1 THEN 'amp' WHEN cnv_direction = -1 THEN 'del' ELSE NULL END) effect
	FROM ref.driver_genes dg
	LEFT JOIN ref.genes rg ON rg.gene_symbol = dg.gene_symbol
	WHERE has_cnv IS TRUE
),
selected_genes_pairs AS
(
	SELECT
		tumor_pair_barcode,
		tumor_barcode_a,
		tumor_barcode_b,
		case_barcode,
		idh_codel_subtype,
		gene_symbol,
		chrom,
		pos,
		cnv_direction,
		effect
	FROM selected_tumor_pairs, selected_genes
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
		sgs.idh_codel_subtype,
		sgs.tumor_barcode_a,
		sgs.tumor_barcode_b,
		sgs.gene_symbol,
		sgs.effect,
		--pathway,
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
		ROUND((c1.wcr - ss1.hldel_fwmean)/ss1.hldel_fwsd,6) AS hldel_zs_a,
		ROUND((c2.wcr - ss2.hldel_fwmean)/ss2.hldel_fwsd,6) AS hldel_zs_b,
		ROUND((c1.wcr - ss1.hlamp_fwmean)/ss1.hlamp_fwsd,6) AS hlamp_zs_a,
		ROUND((c2.wcr - ss2.hlamp_fwmean)/ss2.hlamp_fwsd,6) AS hlamp_zs_b,
		c1.hlvl_call AS cn_a,
		c2.hlvl_call AS cn_b
	FROM selected_genes_pairs sgs
	LEFT JOIN analysis.cnv_by_gene_gatk c1 ON c1.aliquot_barcode = sgs.tumor_barcode_a AND c1.gene_symbol = sgs.gene_symbol
	LEFT JOIN analysis.cnv_by_gene_gatk c2 ON c2.aliquot_barcode = sgs.tumor_barcode_b AND c2.gene_symbol = sgs.gene_symbol
	LEFT JOIN seg_stats_optimized ss1 ON ss1.aliquot_barcode = sgs.tumor_barcode_a
	LEFT JOIN seg_stats_optimized ss2 ON ss2.aliquot_barcode = sgs.tumor_barcode_b
	WHERE abs(c1.hlvl_call) = 2 OR abs(c2.hlvl_call) = 2
),
cnv_by_pair AS
(
	SELECT
		cpg.case_barcode,
		cpg.tumor_pair_barcode,
		cpg.tumor_barcode_a,
		cpg.tumor_barcode_b,
		cpg.idh_codel_subtype,
		bool_and(cnv_change = 'S') AS shared,
		bool_or(cnv_change = 'P') AS private_a,
		bool_or(cnv_change = 'R') AS private_b,
		COUNT(DISTINCT cpg.gene_symbol) AS driver_count,
		COUNT(CASE WHEN cnv_change = 'S' THEN cpg.gene_symbol END) AS driver_count_shared,
		COUNT(CASE WHEN cnv_change = 'P' THEN cpg.gene_symbol END) AS driver_count_private_a,
		COUNT(CASE WHEN cnv_change = 'R' THEN cpg.gene_symbol END) AS driver_count_private_b,
		trim(BOTH ', ' FROM string_agg(CASE WHEN cnv_change = 'S' THEN cpg.gene_symbol || ' ' || effect || ', ' ELSE '' END, '')) AS driver_shared,
		trim(BOTH ', ' FROM string_agg(CASE WHEN cnv_change = 'R' THEN '+' || cpg.gene_symbol || ' ' || effect || ', ' ELSE '' END, '')) AS target_a,
		trim(BOTH ', ' FROM string_agg(CASE WHEN cnv_change = 'P' THEN '-' || cpg.gene_symbol || ' ' || effect || ', ' ELSE '' END, '')) AS target_b,
		(CASE
		 WHEN bool_and(CASE WHEN cnv_change = 'S' THEN gg.idh_codel_subtype IS NOT NULL END) THEN 'In-context'
		 WHEN bool_or(CASE WHEN cnv_change = 'S' THEN gg.idh_codel_subtype IS NULL END) THEN 'Out-of-context'
		 ELSE NULL
		 END) AS driver_context_shared,
		(CASE
		 WHEN bool_and(CASE WHEN cnv_change IN ('P','R') THEN gg.idh_codel_subtype IS NOT NULL END) THEN 'In-context'
		 WHEN bool_or(CASE WHEN cnv_change IN ('P','R')THEN gg.idh_codel_subtype IS NULL END) THEN 'Out-of-context'
		 ELSE NULL
		 END) AS driver_context_change--,
		--(CASE WHEN bool_or(priority = 2 AND NOT is_same_variant) THEN 'Driver convergence' END) AS driver_evolution*/
	FROM cnv_by_pair_gene cpg
	LEFT JOIN ref.gistic_genes gg ON gg.gene_symbol = cpg.gene_symbol AND gg.idh_codel_subtype = cpg.idh_codel_subtype
	GROUP BY 1,2,3,4,5
),
cnv_driver_status AS
(
	SELECT
		case_barcode,
		tumor_pair_barcode,
		tumor_barcode_a,
		tumor_barcode_b,
		driver_count,
		driver_count_shared,
		driver_count_private_a,
		driver_count_private_b,
		driver_shared,
		(CASE WHEN shared THEN 'Driver stable'
		 WHEN NOT shared AND private_a AND NOT private_b THEN 'Driver loss'
		 WHEN NOT shared AND private_b AND NOT private_a THEN 'Driver gain'
		 WHEN NOT shared AND private_b AND private_b THEN 'Driver switch' END) driver_status,
		TRIM(BOTH ', ' FROM target_a || ', ' || target_b) AS target,
		driver_context_shared,
		driver_context_change
	FROM cnv_by_pair
),
cnv_driver_stability AS
(
	SELECT
		stp.tumor_pair_barcode,
		stp.case_barcode,
		stp.tumor_barcode_a,
		stp.tumor_barcode_b,
		idh_codel_subtype,
		(CASE WHEN driver_count > 0 THEN driver_count ELSE 0 END) AS cnv_driver_count,
		driver_count_shared AS cnv_driver_count_shared,
		driver_count_private_a AS cnv_driver_count_private_a,
		driver_count_private_b AS cnv_driver_count_private_b,
		(CASE WHEN driver_shared <> '' THEN driver_shared ELSE NULL END) cnv_driver_shared,
		(CASE WHEN driver_status IS NOT NULL THEN driver_status ELSE 'Driver null' END) cnv_driver_status,
		(CASE
		 WHEN driver_status IN ('Driver switch','Driver loss','Driver gain') THEN 'Driver unstable'
		 WHEN driver_status IN ('Driver stable') OR driver_status IS NULL THEN 'Driver stable' END) cnv_driver_stability,
		(CASE WHEN target <> '' THEN target ELSE NULL END) cnv_driver_change,
		driver_context_shared AS cnv_driver_context_shared,
		driver_context_change AS cnv_driver_context_change
	FROM cnv_driver_status ds
	RIGHT JOIN selected_tumor_pairs stp ON stp.tumor_pair_barcode = ds.tumor_pair_barcode
)
SELECT * FROM cnv_driver_stability ORDER BY 1
			 
--SELECT * FROM variants_by_case
/*SELECT
    (CASE WHEN driver_status IN ('Driver gain','Driver loss','Driver switch') THEN 'Driver change' ELSE 'Driver stable' END),
    recurrence_evolution,
    COUNT(*)
FROM driver_stability ds
LEFT JOIN analysis.neutrality_tumor_pairs ntp ON ntp.tumor_pair_barcode = ds.tumor_pair_barcode
INNER JOIN analysis.silver_set ss ON ss.tumor_pair_barcode = ds.tumor_pair_barcode
WHERE driver_status IS NOT NULL AND recurrence_evolution IS NOT NULL
GROUP BY 1,2*/
			 
-- END --