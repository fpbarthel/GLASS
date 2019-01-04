WITH selected_tumor_pairs AS
(
	SELECT
		tumor_pair_barcode,
		case_barcode,
		tumor_barcode_a,
		tumor_barcode_b,
		row_number() OVER (PARTITION BY case_barcode ORDER BY surgical_interval_mo DESC, portion_a ASC, portion_b ASC, substring(tumor_pair_barcode from 27 for 3) ASC) AS priority
	FROM analysis.tumor_pairs ps
	LEFT JOIN analysis.blocklist b1 ON b1.aliquot_barcode = ps.tumor_barcode_a
	LEFT JOIN analysis.blocklist b2 ON b2.aliquot_barcode = ps.tumor_barcode_b
	WHERE
		comparison_type = 'longitudinal' AND
		sample_type_b <> 'M1' AND 													-- exclude metastatic samples here because this is outside the scope of our study
		b1.fingerprint_exclusion = 'allow' AND b2.fingerprint_exclusion = 'allow' AND
		b1.coverage_exclusion = 'allow' AND b2.coverage_exclusion = 'allow'
)
SELECT
	mgt.tumor_pair_barcode,
	(CASE WHEN mutect2_call_a AND mutect2_call_b THEN 'S' WHEN mutect2_call_a AND NOT mutect2_call_b THEN 'P' WHEN mutect2_call_b AND NOT mutect2_call_a THEN 'R' END) AS variant_status,
	mgt.case_barcode,
	mgt.tumor_barcode_a,
	mgt.tumor_barcode_b,
	mgt.gene_symbol,
	mgt.variant_type,
	mgt.variant_classification,
	mgt.chrom::varchar(2),
	lower(mgt.pos) AS start_pos,
	upper(mgt.pos) -1 AS end_pos,
	ref,
	mgt.alt,
	ref_count_a,
	ref_count_b,
	alt_count_a,
	alt_count_b,
	mutect2_call_a,
	mutect2_call_b,
	vaf_corrected_call_a,
	vaf_corrected_call_b,
	logr_copy_number_a,
	logr_copy_number_b,
	corrected_copy_number_a,
	corrected_copy_number_b,
	corrected_call_a::varchar(5),
	corrected_call_b::varchar(5)
FROM analysis.master_genotype_comparison mgt
LEFT JOIN analysis.snvs snvs ON snvs.chrom = mgt.chrom AND snvs.pos = mgt.pos AND snvs.alt = mgt.alt
WHERE (mutect2_call_a OR mutect2_call_b) AND (ref_count_a + alt_count_a) >= 10 AND (ref_count_b + alt_count_b) >= 10
--INNER JOIN selected_tumor_pairs stp ON stp.tumor_pair_barcode = mgt.tumor_pair_barcode