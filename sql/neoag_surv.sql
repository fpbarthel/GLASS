WITH t1 AS
(
	SELECT gt.tumor_pair_barcode,gt.chrom,gt.pos,gt.alt,gt.gene_symbol,neo.fraction,gt.tumor_barcode_a,gt.tumor_barcode_b,clin.surgical_interval,netmhcpan_mt_score,received_tmz,received_rt,hypermutator_status,surg.treatment_chemotherapy_other,surg.idh_codel_subtype,case_age_diagnosis_years,case_vital_status,case_overall_survival_mo,row_number() OVER w AS neoag_number,
	(CASE WHEN netmhcpan_mt_score IS NULL THEN 0 WHEN netmhcpan_mt_score IS NOT NULL THEN 1 END) AS neoantigen
	FROM variants.pgeno gt
	INNER JOIN analysis.silver_set ss ON ss.tumor_pair_barcode = gt.tumor_pair_barcode
	INNER JOIN analysis.neoantigens_by_pair neo ON 
		gt.tumor_pair_barcode = neo.tumor_pair_barcode AND 
		gt.chrom = CAST(neo.chrom AS int) AND 
		gt.pos = neo.pos AND
		gt.alt = neo.alt
	LEFT JOIN analysis.tumor_clinical_comparison clin ON
		gt.tumor_pair_barcode = clin.tumor_pair_barcode
	LEFT JOIN clinical.surgeries surg ON
		substr(gt.tumor_barcode_a,1,15) = surg.sample_barcode
	LEFT JOIN clinical.cases surv ON
		gt.case_barcode = surv.case_barcode
	WHERE (mutect2_call_a OR mutect2_call_b) AND ((ref_count_a + alt_count_a) >= 15 OR (ref_count_b + alt_count_b) >= 15) 
	WINDOW w AS (PARTITION BY gt.tumor_pair_barcode,gt.chrom,gt.pos,gt.alt ORDER BY netmhcpan_mt_score)
)
SELECT * FROM t1 WHERE neoag_number = 1 ORDER BY tumor_pair_barcode

