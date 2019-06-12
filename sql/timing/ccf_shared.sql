SELECT 
		pg.tumor_pair_barcode,
		pg.case_barcode,
		st.idh_codel_subtype,
		pg.tumor_barcode_a,
		pg.tumor_barcode_b,
		hypermutator_status,
		pg.chrom, 
		pg.pos,
		pg.variant_id,
		pg.gene_symbol,
		pg.variant_classification,
		vc.variant_effect,
		vc.variant_classification_vep,
		pl1.cellular_prevalence AS cellular_prevalence_a, 
		pl1.variant_allele_frequency AS variant_allele_frequency_a, 
		(CASE WHEN pl1.cellular_prevalence >= 0.5 THEN 'C' WHEN pl1.cellular_prevalence >= 0.1 AND pl1.cellular_prevalence < 0.5 THEN 'S' ELSE 'ND' END) AS clonality_a,
		pl2.cellular_prevalence AS cellular_prevalence_b, 
		pl2.variant_allele_frequency AS variant_allele_frequency_b, 
		(CASE WHEN pl2.cellular_prevalence >= 0.5 THEN 'C' WHEN pl2.cellular_prevalence >= 0.1 AND pl2.cellular_prevalence < 0.5 THEN 'S' ELSE 'ND' END) AS clonality_b,
		rank() OVER (PARTITION BY pg.tumor_pair_barcode, pg.gene_symbol ORDER BY variant_classification_priority, pl1.cellular_prevalence + pl2.cellular_prevalence DESC)
	FROM variants.pgeno pg
	LEFT JOIN variants.pyclone_loci pl1 ON pl1.variant_id = pg.variant_id AND pl1.aliquot_barcode = pg.tumor_barcode_a 
	LEFT JOIN variants.pyclone_loci pl2 ON pl2.variant_id = pg.variant_id AND pl2.aliquot_barcode= pg.tumor_barcode_b
	LEFT JOIN variants.variant_classifications vc ON vc.variant_classification = pg.variant_classification
	INNER JOIN analysis.gold_set ss ON pg.tumor_pair_barcode = ss.tumor_pair_barcode
	INNER JOIN analysis.tumor_mut_comparison_anno tmc ON tmc.tumor_pair_barcode = ss.tumor_pair_barcode
	INNER JOIN clinical.subtypes st ON st.case_barcode = ss.case_barcode
	WHERE pl1.cellular_prevalence IS NOT NULL AND mutect2_call_a AND mutect2_call_b AND variant_classification_priority IS NOT NULL