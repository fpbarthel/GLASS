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
		has_mut IS TRUE AND sn.gene_symbol IN ('IDH1','ATRX','TP53') AND
		((sn.gene_symbol NOT IN ('TERT','IDH1') AND variant_classification_priority IS NOT NULL) OR
		(sn.gene_symbol = 'TERT' AND sn.variant_classification = 'FIVE_PRIME_FLANK' AND lower(sn.pos) IN (1295228,1295250)) OR
		(sn.gene_symbol = 'IDH1' AND sn.protein_change IN ('p.R132C','p.R132G','p.R132H','p.R132S')))
),
t1 AS 
(
	SELECT
		pl.aliquot_barcode,
		sq.case_barcode,
		sg.variant_id,
		idh_codel_subtype,sample_type,gene_symbol,
		--protein_change,
		titan_ccf,pyclone_ccf,
		row_number() OVER (PARTITION BY pl.aliquot_barcode,gene_symbol ORDER BY pyclone_ccf DESC) AS optimal_variant
		--COUNT(DISTINCT gene_symbol) OVER (PARTITION BY aliquot_barcode ORDER BY gene_symbol) = 3 AS has_three_mut,
		--array_agg(DISTINCT(gene_symbol) OVER (PARTITION BY aliquot_barcode ORDER BY gene_symbol)) = '{"ATRX","IDH1","TP53"}' AS all_three,
		--rank() OVER (PARTITION BY pl.aliquot_barcode ORDER BY cellular_prevalence DESC) AS mut_order,
		--COUNT(*) OVER (PARTITION BY gene_symbol,idh_codel_subtype,sample_type) AS n_mut_group
	FROM variants.pyclone_loci pl
	INNER JOIN selected_genes sg ON sg.variant_id = pl.variant_id
	INNER JOIN selected_aliquots sq ON sq.aliquot_barcode = pl.aliquot_barcode
	INNER JOIN variants.passgeno pg ON pg.variant_id = pl.variant_id AND pg.aliquot_barcode = pl.aliquot_barcode
	INNER JOIN clinical.subtypes st ON st.case_barcode = pg.case_barcode
	--WHERE ssm2_pass_call
),
t2 AS
(
	SELECT *
	FROM t1
	WHERE optimal_variant = 1
),
t3 AS
(
	SELECT *, COUNT(*) OVER (PARTITION BY aliquot_barcode) AS num_genes
	FROM t2
),
t4 AS
(
	-- You can edit the number of decimals to round by to change how ties are ranked
	-- By rounding to two decimals, a CCF of 0.977 and 0.983 are considered identical
	-- And two mutations with these CCFs would be given the same rank
	-- If one instead rounds to three decimals, these would not be considered identical
	SELECT
		*,
		(RANK() OVER (PARTITION BY aliquot_barcode ORDER BY round(pyclone_ccf::decimal,2) DESC))::integer AS vrank,
		COUNT(*) OVER (PARTITION BY gene_symbol,sample_type) AS num_samples
	FROM t3
	WHERE num_genes = 3
),
t5 AS
(
	SELECT sample_type,gene_symbol,vrank,num_samples,COUNT(*),COUNT(*)::decimal/num_samples
	FROM t4
	GROUP BY 1,2,3,4
	ORDER BY 1,2,3,4
)
--SELECT gene_symbol,idh_codel_subtype,sample_type,SUM(mut_order)/COUNT(mut_order),COUNT(mut_order)
--FROM timing_snv GROUP BY 1,2,3 ORDER BY 2,3,4
--SELECT gene_symbol,idh_codel_subtype,sample_type,mut_order,n_mut_group,COUNT(mut_order)::decimal/n_mut_group
--FROM timing_snv
--GROUP BY 1,2,3,4,5
--ORDER BY 2,3,4
SELECT * FROM t4