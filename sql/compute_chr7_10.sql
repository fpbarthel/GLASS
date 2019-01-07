/*
Calculate chromosome 7/10 status
*/
WITH
selected_tumor_pairs AS
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
		b1.coverage_exclusion <> 'block' AND b2.cnv_exclusion <> 'block' 
),
t2 AS
(
	SELECT
		*,
		(SELECT log2_copy_ratio FROM analysis.cnv_by_chr_gatk WHERE chrom = '7' AND aliquot_barcode = stp.tumor_barcode_a) AS chr7_logr_a,
		(SELECT log2_copy_ratio FROM analysis.cnv_by_chr_gatk WHERE chrom = '7' AND aliquot_barcode = stp.tumor_barcode_b) AS chr7_logr_b,
		(SELECT log2_copy_ratio FROM analysis.cnv_by_chr_gatk WHERE chrom = '10' AND aliquot_barcode = stp.tumor_barcode_a) AS chr10_logr_a,
		(SELECT log2_copy_ratio FROM analysis.cnv_by_chr_gatk WHERE chrom = '10' AND aliquot_barcode = stp.tumor_barcode_b) AS chr10_logr_b
	FROM selected_tumor_pairs stp
	WHERE priority = 1
),
t3 AS
(
	SELECT
		*,
		(CASE WHEN chr7_logr_a > 0.1 AND chr10_logr_a < -0.1 THEN 1 ELSE 0 END)::boolean AS chr7_10_a,
		(CASE WHEN chr7_logr_b > 0.1 AND chr10_logr_b < -0.1 THEN 1 ELSE 0 END)::boolean AS chr7_10_b
	FROM t2
)
SELECT
	tumor_pair_barcode,
	case_barcode,
	tumor_barcode_a,
	tumor_barcode_b,
	(CASE WHEN chr7_10_a AND chr7_10_b THEN 'shared'
		  WHEN chr7_10_a AND NOT chr7_10_b THEN 'shed'
		  WHEN chr7_10_b AND NOT chr7_10_a THEN 'acquired'
		  WHEN NOT chr7_10_a AND NOT chr7_10_b THEN 'no' END) AS chr7_10_status
FROM t3