WITH
selected_tumor_pairs AS
(
	SELECT * FROM analysis.silver_set
),
selected_aliquots AS
(
	SELECT tumor_barcode_a AS aliquot_barcode, case_barcode, 'P' AS sample_type FROM selected_tumor_pairs
	UNION
	SELECT tumor_barcode_b AS aliquot_barcode, case_barcode, 'R' AS sample_type FROM selected_tumor_pairs
),
selected_genes AS
(
	SELECT DISTINCT sn.gene_symbol, ensembl_gene_id, variant_id, chrom, pos, alt, sn.variant_classification, variant_classification_priority, protein_change
	FROM variants.passanno sn
	INNER JOIN ref.driver_genes ds ON ds.gene_symbol = sn.gene_symbol
	INNER JOIN ref.ensembl_gene_mapping gm ON gm.gene_symbol = sn.gene_symbol
	LEFT JOIN variants.variant_classifications vc ON sn.variant_classification = vc.variant_classification
	WHERE
		has_mut IS TRUE AND
		((sn.gene_symbol NOT IN ('TERT','IDH1') AND variant_classification_priority IS NOT NULL) OR
		(sn.gene_symbol = 'TERT' AND sn.variant_classification = 'FIVE_PRIME_FLANK' AND lower(sn.pos) IN (1295228,1295250)) OR
		(sn.gene_symbol = 'IDH1' AND sn.protein_change IN ('p.R132C','p.R132G','p.R132H','p.R132S')))
),
timing_snv AS 
(
	SELECT pl.aliquot_barcode,idh_codel_subtype,sample_type,gene_symbol,variant_classification,protein_change,cellular_prevalence,titan_ccf,pyclone_ccf, rank() OVER (PARTITION BY pl.aliquot_barcode ORDER BY cellular_prevalence DESC) AS mut_order
	FROM variants.pyclone_loci pl
	INNER JOIN selected_genes sg ON sg.variant_id = pl.variant_id
	INNER JOIN selected_aliquots sq ON sq.aliquot_barcode = pl.aliquot_barcode
	INNER JOIN variants.passgeno pg ON pg.variant_id = pl.variant_id AND pg.aliquot_barcode = pl.aliquot_barcode
	INNER JOIN clinical.subtypes st ON st.case_barcode = pg.case_barcode
	WHERE ssm2_pass_call
)
SELECT gene_symbol,idh_codel_subtype,sample_type,SUM(mut_order)/COUNT(mut_order),COUNT(mut_order)
FROM timing_snv GROUP BY 1,2,3 ORDER BY 2,3,4
--SELECT gene_symbol, idh_codel_subtype, sample_type, COUNT(cellular_prevalence) AS num_mut, SUM(cellular_prevalence)/COUNT(cellular_prevalence) AS mean_cp, SUM(cellular_prevalence_sd)/COUNT(cellular_prevalence_sd) AS mean_cp_sd FROM timing_snv
--GROUP BY 1,2,3
--ORDER BY 1,2,3