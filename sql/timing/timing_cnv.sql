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
cnv_timing AS
(
	SELECT gc.gene_symbol,idh_codel_subtype,sample_type,hlvl_call, COUNT(cellular_prevalence) AS num_cp, SUM(cellular_prevalence)/COUNT(cellular_prevalence) AS mean_cp
	FROM analysis.gatk_cnv_by_gene gc
	INNER JOIN ref.driver_genes dg ON dg.gene_symbol = gc.gene_symbol
	INNER JOIN selected_aliquots sa ON sa.aliquot_barcode = gc.aliquot_barcode
	INNER JOIN clinical.subtypes st ON st.case_barcode = sa.case_barcode
	GROUP BY 1,2,3,4
)
SELECT * FROM cnv_timing