WITH t1 AS
(
	SELECT gt.tumor_pair_barcode,gt.chrom,gt.pos,gt.alt,gt.gene_symbol,neo.fraction,gt.tumor_barcode_a,gt.tumor_barcode_b,gt.surgical_interval_mo,netmhcpan_mt_score,received_tmz,received_rt,hypermutator_status,surg.treatment_chemotherapy_other,surg.idh_codel_subtype,snv_driver_status,snv_driver_change_b, row_number() OVER w AS neoag_number,
	(CASE WHEN netmhcpan_mt_score IS NULL THEN 0 WHEN netmhcpan_mt_score IS NOT NULL THEN 1 END) AS neoantigen,
	CONCAT('-',gene_name,' p.',substr(mutation,1,1),pvacseq_protein_position,substr(mutation,3,1)) AS aa_change	
	FROM analysis.master_genotype_comparison gt
	INNER JOIN analysis.silver_set ss ON ss.tumor_pair_barcode = gt.tumor_pair_barcode
	INNER JOIN analysis.pvacseq_fraction neo ON 
		gt.tumor_pair_barcode = neo.tumor_pair_barcode AND 
		gt.chrom = neo.chrom AND 
		gt.pos = neo.pos AND
		gt.alt = neo.alt
	LEFT JOIN clinical.clinical_by_tumor_pair clin ON
		gt.tumor_pair_barcode = clin.tumor_pair_barcode
	LEFT JOIN clinical.surgeries surg ON
		substr(gt.tumor_barcode_a,1,15) = surg.sample_barcode
	LEFT JOIN analysis.driver_status_snv driv ON
		gt.tumor_pair_barcode = driv.tumor_pair_barcode
	WHERE (mutect2_call_a OR mutect2_call_b) AND ((ref_count_a + alt_count_a) >= 15 OR (ref_count_b + alt_count_b) >= 15) AND (driv.snv_driver_status = 'Driver loss' OR driv.snv_driver_status = 'Driver switch')
	WINDOW w AS (PARTITION BY gt.tumor_pair_barcode,gt.chrom,gt.pos,gt.alt ORDER BY netmhcpan_mt_score)
)
SELECT * FROM t1 WHERE neoag_number = 1 AND aa_change = snv_driver_change_b ORDER BY tumor_pair_barcode
