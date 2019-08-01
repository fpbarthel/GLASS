/*
Query that creates table of all shared variants with their pyclone ccfs and clonality by pair, as well as their neoantigen status
- This table is directly used to make Extended Data Figure 12A
- Only examines variants that are not labelled as MODIFIER, (coding variants), as neoantigens can only come from coding variants
- Like other queries, this query needs to be redesigned to use VEP classifications once VEP calls for all variants exist
- Each variant is marked as immunogenic if it gives rise to at least one neoantigen
- No coverage filters are applied in this query, only mutect2 calls
*/

WITH start_tab AS
(
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
			pv.variant_classification,
			vc.variant_effect,
			neo.hla_allele,
			neo.netmhcpan_mt_score,
			pl1.cellular_prevalence AS cellular_prevalence_a, 
			pl1.variant_allele_frequency AS variant_allele_frequency_a, 
			(CASE WHEN pl1.cellular_prevalence >= 0.5 THEN 'C' WHEN pl1.cellular_prevalence >= 0.1 AND pl1.cellular_prevalence < 0.5 THEN 'S' ELSE 'ND' END) AS clonality_a,
			pl2.cellular_prevalence AS cellular_prevalence_b, 
			pl2.variant_allele_frequency AS variant_allele_frequency_b, 
			(CASE WHEN pl2.cellular_prevalence >= 0.5 THEN 'C' WHEN pl2.cellular_prevalence >= 0.1 AND pl2.cellular_prevalence < 0.5 THEN 'S' ELSE 'ND' END) AS clonality_b,
			rank() OVER (PARTITION BY pg.tumor_pair_barcode, pg.gene_symbol ORDER BY variant_classification_priority, pl1.cellular_prevalence + pl2.cellular_prevalence DESC),
			row_number() OVER w AS var_number
		FROM variants.pgeno pg
		LEFT JOIN variants.passvep pv ON pv.variant_id = pg.variant_id
		LEFT JOIN variants.pyclone_loci pl1 ON pl1.variant_id = pg.variant_id AND pl1.aliquot_barcode = pg.tumor_barcode_a 
		LEFT JOIN variants.pyclone_loci pl2 ON pl2.variant_id = pg.variant_id AND pl2.aliquot_barcode= pg.tumor_barcode_b
		LEFT JOIN variants.variant_classifications vc ON vc.variant_classification_vep = pv.variant_classification
		LEFT JOIN analysis.neoantigens_by_pair neo ON neo.tumor_pair_barcode = pg.tumor_pair_barcode AND neo.variant_id = pg.variant_id
		INNER JOIN analysis.gold_set ss ON pg.tumor_pair_barcode = ss.tumor_pair_barcode
		INNER JOIN analysis.tumor_mut_comparison_anno tmc ON tmc.tumor_pair_barcode = ss.tumor_pair_barcode
		INNER JOIN clinical.subtypes st ON st.case_barcode = ss.case_barcode
		WHERE pl1.cellular_prevalence IS NOT NULL AND mutect2_call_a AND mutect2_call_b AND variant_classification_impact != 'MODIFIER'
		WINDOW w AS (PARTITION BY pg.tumor_pair_barcode, pg.variant_id ORDER BY neo.netmhcpan_mt_score)
)
SELECT *,
CASE WHEN netmhcpan_mt_score IS NOT NULL THEN 1 WHEN netmhcpan_mt_score IS NULL then 0 END AS is_neoag
FROM start_tab
WHERE var_number = 1
