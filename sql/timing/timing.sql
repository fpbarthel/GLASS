WITH
selected_tumor_pairs AS
(
	SELECT tumor_pair_barcode, tumor_barcode_a, tumor_barcode_b, ss.case_barcode, idh_codel_subtype
	FROM analysis.gold_set ss
	INNER JOIN clinical.subtypes st ON st.case_barcode = ss.case_barcode
),
selected_aliquots AS
(
	SELECT tumor_barcode_a AS aliquot_barcode, idh_codel_subtype, case_barcode, 'P' AS sample_type FROM selected_tumor_pairs
	UNION
	SELECT tumor_barcode_b AS aliquot_barcode, idh_codel_subtype, case_barcode, 'R' AS sample_type FROM selected_tumor_pairs
),
selected_snvs AS -- Select specific gene mutations
(
	SELECT DISTINCT sn.gene_symbol, idh_codel_subtype, variant_id, chrom, pos, alt, vc.variant_classification_vep AS variant_classification, variant_classification_priority, protein_change
	FROM variants.passanno sn
	INNER JOIN ref.snv_drivers_subtype ds ON ds.gene_symbol = sn.gene_symbol
	LEFT JOIN variants.variant_classifications vc ON sn.variant_classification = vc.variant_classification
	WHERE
		((sn.gene_symbol NOT IN ('TERT','IDH1') AND variant_classification_priority IS NOT NULL) OR
		(sn.gene_symbol = 'TERT' AND sn.variant_classification = 'FIVE_PRIME_FLANK' AND lower(sn.pos) IN (1295228,1295250)) OR
		(sn.gene_symbol = 'IDH1' AND sn.protein_change IN ('p.R132C','p.R132G','p.R132H','p.R132S')))
),
selected_cnv AS
(
	SELECT gene_symbol, idh_codel_subtype, direction
	FROM ref.cnv_drivers_subtype
),
selected_arms AS
(
	SELECT chrom, arm, idh_codel_subtype, direction
	FROM ref.arm_drivers_subtype
),
filter_snv AS 
(
	SELECT
		pg.case_barcode, pl.aliquot_barcode, sq.idh_codel_subtype, sample_type,
		gene_symbol, variant_classification, protein_change,
		titan_ccf, pyclone_ccf,
		rank() OVER (PARTITION BY pl.aliquot_barcode, gene_symbol ORDER BY cellular_prevalence DESC) AS rnk
	FROM variants.pyclone_loci pl
	INNER JOIN selected_aliquots sq ON sq.aliquot_barcode = pl.aliquot_barcode
	INNER JOIN selected_snvs sg ON sg.variant_id = pl.variant_id AND sg.idh_codel_subtype = sq.idh_codel_subtype
	INNER JOIN variants.passgeno pg ON pg.variant_id = pl.variant_id AND pg.aliquot_barcode = pl.aliquot_barcode
	WHERE ssm2_pass_call
),
timing_snv AS 
(
	SELECT
		case_barcode, idh_codel_subtype, sample_type,
		gene_symbol, variant_classification, protein_change,
		titan_ccf, pyclone_ccf,
		rank() OVER (PARTITION BY aliquot_barcode ORDER BY pyclone_ccf DESC) AS mut_order
	FROM filter_snv
	WHERE rnk = 1
),
timing_cnv AS
(
	SELECT
		sq.case_barcode, sq.idh_codel_subtype, sample_type,
		sc.gene_symbol, hlvl_call, cellular_prevalence, sc.direction
	FROM analysis.gatk_cnv_by_gene cnv
	INNER JOIN selected_aliquots sq ON sq.aliquot_barcode = cnv.aliquot_barcode
	INNER JOIN selected_cnv sc ON sc.gene_symbol = cnv.gene_symbol AND sc.idh_codel_subtype = sq.idh_codel_subtype AND sc.direction = cnv.hlvl_call
),
timing_arm AS -- Computes average TITAN CCF for each arm
(
	SELECT
		sq.case_barcode, sa.idh_codel_subtype, sq.sample_type,
		cnv.chrom, cnv.arm,
		arm_call, ca.pos, arm_num_seg, sa.direction,
		(upper(ca.pos)-lower(ca.pos)-1) AS arm_size,
		sum(CASE WHEN cellular_prevalence IS NOT NULL THEN (upper(ts.pos)-lower(ts.pos)-1) END)::decimal / (upper(ca.pos)-lower(ca.pos)-1) AS sum_seg_size,
		sum(CASE WHEN cellular_prevalence IS NOT NULL THEN (upper(ts.pos)-lower(ts.pos)-1) * cellular_prevalence END) / sum(CASE WHEN cellular_prevalence IS NOT NULL THEN (upper(ts.pos)-lower(ts.pos)-1) END)::decimal AS arm_ccf
	FROM analysis.gatk_cnv_by_arm cnv
	INNER JOIN selected_aliquots sq ON sq.aliquot_barcode = cnv.aliquot_barcode
	INNER JOIN selected_arms sa ON sa.chrom = cnv.chrom AND sa.arm = cnv.arm AND sa.idh_codel_subtype = sq.idh_codel_subtype AND sa.direction = cnv.arm_call
	INNER JOIN ref.chr_arms ca ON ca.chrom = cnv.chrom AND ca.arm = cnv.arm
	INNER JOIN variants.titan_seg ts ON ts.aliquot_barcode = cnv.aliquot_barcode AND ts.chrom = cnv.chrom AND ts.pos && ca.pos
	WHERE (upper(ts.pos)-lower(ts.pos)-1) > 0
	GROUP BY 1,2,3,4,5,6,7,8,9
),
timing_all AS
(
	SELECT case_barcode, idh_codel_subtype, sample_type, variant_classification, gene_symbol || ' mut' AS evnt, pyclone_ccf AS ccf FROM timing_snv WHERE pyclone_ccf IS NOT NULL
	UNION
	SELECT case_barcode, idh_codel_subtype, sample_type, (CASE direction WHEN -2 THEN 'Deletion' WHEN 2 THEN 'Amplification' ELSE NULL END) AS variant_classification, gene_symbol || (CASE direction WHEN -2 THEN ' del' WHEN 2 THEN ' amp' ELSE NULL END) AS evnt, cellular_prevalence AS ccf FROM timing_cnv WHERE cellular_prevalence IS NOT NULL
	UNION
	SELECT case_barcode, idh_codel_subtype, sample_type, (CASE direction WHEN -1 THEN 'Loss' WHEN 1 THEN 'Gain' ELSE NULL END) AS variant_classification, arm || (CASE direction WHEN -1 THEN ' del' WHEN 1 THEN ' amp' ELSE NULL END) AS evnt, arm_ccf AS ccf FROM timing_arm WHERE arm_ccf IS NOT NULL
),
timing_all_ranked AS
(
	SELECT
		*,
		dense_rank() OVER (PARTITION BY case_barcode ORDER BY round(ccf::numeric,1) DESC) AS evnt_rank,
		COUNT(*) OVER (PARTITION BY case_barcode, evnt) AS num_samples
	FROM timing_all
	ORDER BY 1
)/*,
timing_snv_agg AS
(
	SELECT gene_symbol || ' mut' AS evnt, idh_codel_subtype, sample_type, COUNT(CASE WHEN pyclone_ccf > 0.5 THEN 1 END) AS n_clonal, COUNT(pyclone_ccf) AS n_total, COUNT(CASE WHEN pyclone_ccf > 0.5 THEN 1 END) / COUNT(pyclone_ccf)::decimal AS prop_clonal
	FROM timing_snv
	GROUP BY 1,2,3
),
timing_cnv_agg AS
(
	SELECT gene_symbol || (CASE direction WHEN -2 THEN ' del' WHEN 2 THEN ' amp' ELSE NULL END) AS evnt, idh_codel_subtype, sample_type, COUNT(CASE WHEN cellular_prevalence > 0.5 THEN 1 END) AS n_clonal, COUNT(cellular_prevalence) AS n_total, COUNT(CASE WHEN cellular_prevalence > 0.5 THEN 1 END) / COUNT(cellular_prevalence)::decimal AS prop_clonal
	FROM timing_cnv
	GROUP BY 1,2,3
),
timing_arm_agg AS
(
	SELECT arm || (CASE direction WHEN -1 THEN ' del' WHEN 1 THEN ' amp' ELSE NULL END)  AS evnt, idh_codel_subtype, sample_type, COUNT(CASE WHEN arm_ccf > 0.5 THEN 1 END) AS n_clonal, COUNT(arm_ccf) AS n_total, COUNT(CASE WHEN arm_ccf > 0.5 THEN 1 END) / COUNT(arm_ccf)::decimal AS prop_clonal
	FROM timing_arm
	GROUP BY 1,2,3
),
timing_agg AS
(
	SELECT * FROM timing_snv_agg
	UNION
	SELECT * FROM timing_cnv_agg
	UNION
	SELECT * FROM timing_arm_agg
)*/
/*SELECT
	evnt,
	subtype,
	sample_type,
	n_clonal::smallint,
	n_total::smallint,
	prop_clonal
FROM timing_agg*/
/*SELECT subtype,sample_type,evnt,COUNT(evnt_rank)::integer AS num_events,sum(evnt_rank)/COUNT(evnt_rank) AS mean_rank, median(evnt_rank) AS med_rank
FROM timing_all_ranked
GROUP BY 1,2,3
ORDER BY 1,2,6*/
SELECT * FROM timing_all_ranked
		 
-- END --