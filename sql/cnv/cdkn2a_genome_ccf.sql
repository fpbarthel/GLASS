WITH
seg_wccf AS
(
	SELECT
		aliquot_barcode,
		COUNT(*) AS num_seg,
		sum(CASE WHEN cellular_prevalence IS NOT NULL THEN (upper(pos) - lower(pos) -1) * cellular_prevalence END) / sum(CASE WHEN cellular_prevalence IS NOT NULL THEN upper(pos) - lower(pos) -1 END)::decimal AS wccf,
		max(clonal_cluster) AS num_clusters,
		sum(CASE WHEN copy_number <> 2 THEN upper(pos) - lower(pos) -1 END) / sum(upper(pos) - lower(pos) -1)::decimal AS titan_aneuploidy
	FROM variants.titan_seg
	GROUP BY 1
),
cdkn2a_call AS
(
	SELECT aliquot_barcode, hlvl_call AS cdkn2a_call, cellular_prevalence AS cdkn2a_ccf
	FROM analysis.gatk_cnv_by_gene
	WHERE gene_symbol = 'CDKN2A'
),
selected_regions AS
(
	SELECT '10q25-26' AS region, * FROM ref.cytobands WHERE chrom = 10 AND substring(cytoband from 1 for 3) IN ('q25','q26')
),
gene_seg_intersect AS
(
    SELECT aliquot_barcode, region, gs.chrom, (upper(t0.pos * gs.pos) - lower(t0.pos * gs.pos) -1) AS w, 2^log2_copy_ratio::decimal As cr
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
),
calls AS
(
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
		 END) c10q25_26_call,
		gc.wcr,
		wcp AS c10q25_26_ccf
	FROM gene_sample_call gc
	LEFT JOIN seg_stats_optimized ss ON ss.aliquot_barcode = gc.aliquot_barcode
	LEFT JOIN gene_cp_agg cp ON cp.aliquot_barcode = gc.aliquot_barcode AND cp.region = gc.region
	ORDER BY 3
)
SELECT seg.aliquot_barcode, idh_codel_subtype, num_seg, num_clusters, wccf, cnv_exclusion, cdkn2a_call, cdkn2a_ccf, titan_aneuploidy, prop_aneuploidy, aneuploidy_score, c10q25_26_call, c10q25_26_ccf
FROM seg_wccf seg
INNER JOIN cdkn2a_call cc ON cc.aliquot_barcode = seg.aliquot_barcode
INNER JOIN analysis.gatk_aneuploidy ga ON ga.aliquot_barcode = seg.aliquot_barcode
INNER JOIN clinical.subtypes su ON su.case_barcode = substring(seg.aliquot_barcode from 1 for 12)
INNER JOIN analysis.blocklist bl ON bl.aliquot_barcode = seg.aliquot_barcode
INNER JOIN calls ca ON ca.aliquot_barcode = seg.aliquot_barcode
ORDER BY 1