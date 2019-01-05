/*
Determine driver changes between selected primary and recurrences
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
		b1.coverage_exclusion = 'allow' AND b2.coverage_exclusion = 'allow' 
),
selected_aliquots AS
(
	SELECT tumor_barcode_a AS aliquot_barcode, case_barcode FROM selected_tumor_pairs WHERE priority = 1
	UNION
	SELECT tumor_barcode_b AS aliquot_barcode, case_barcode FROM selected_tumor_pairs WHERE priority = 1
),
selected_genes AS
(
	SELECT sn.gene_symbol, chrom, pos, alt, sn.variant_classification, variant_classification_priority, hgvs_p
	FROM analysis.snvs sn
	INNER JOIN analysis.dndscv_gene ds ON ds.gene_symbol = sn.gene_symbol AND (ds.qglobal_cv < 0.10 OR ds.gene_symbol IN ('TERT')) -- ,'IDH2','NOTCH1','PDGFRA','PIK3CG','BRAF','H3F3A'))
	LEFT JOIN analysis.variant_classifications vc ON sn.variant_classification = vc.variant_classification
	WHERE
		(sn.gene_symbol NOT IN ('TERT') AND variant_classification_priority IS NOT NULL) OR -- ,'IDH1','IDH2','BRAF','H3F3A') AND variant_classification_priority IS NOT NULL) OR 
		(sn.gene_symbol = 'TERT' AND sn.variant_classification = '5''Flank' AND lower(sn.pos) IN (1295228,1295250)) OR
		(sn.gene_symbol = 'IDH1' AND sn.hgvs_p IN ('p.R132C','p.R132G','p.R132H','p.R132S'))
),
selected_genes_samples AS
(
	SELECT aliquot_barcode, case_barcode, gene_symbol, chrom, lower(pos) AS start_pos, upper(pos)-1 AS end_pos, alt
	FROM selected_aliquots, selected_genes
),
gene_sample_coverage AS
(
	SELECT
		gene_symbol,
		sg.aliquot_barcode,
		sg.case_barcode,
		round(median(alt_count + ref_count)) AS median_cov
	FROM analysis.full_genotypes fgt
	INNER JOIN selected_genes_samples sg ON fgt.aliquot_barcode = sg.aliquot_barcode AND fgt.chrom = sg.chrom AND fgt.start = sg.start_pos AND fgt.end = sg.end_pos AND fgt.alt = sg.alt
	GROUP BY 1,2,3
),
gene_pair_coverage AS
(
	SELECT
		stp.tumor_pair_barcode,
		stp.case_barcode,
		c1.gene_symbol,
		c1.median_cov AS cov_a,
		c2.median_cov AS cov_b
	FROM selected_tumor_pairs stp
	LEFT JOIN gene_sample_coverage c1 ON c1.aliquot_barcode = stp.tumor_barcode_a
	LEFT JOIN gene_sample_coverage c2 ON c2.aliquot_barcode = stp.tumor_barcode_b AND c1.gene_symbol = c2.gene_symbol
	WHERE priority = 1
),
variants_by_case_and_gene AS
(
	SELECT
		gtc.gene_symbol,
		gtc.case_barcode,
		gtc.tumor_pair_barcode,
		gtc.tumor_barcode_a,
		gtc.tumor_barcode_b,
		gtc.chrom,
		gtc.pos,
		gtc.alt,
		gtc.variant_classification,
		sg.hgvs_p,
		(alt_count_a > 0) AS selected_call_a,
		(alt_count_b > 0) AS selected_call_b,
		row_number() OVER (PARTITION BY gtc.gene_symbol, gtc.case_barcode ORDER BY mutect2_call_a::integer + mutect2_call_b::integer DESC, variant_classification_priority, mutect2_call_a::integer + mutect2_call_b::integer DESC, read_depth_a + read_depth_b DESC) AS priority
	FROM analysis.master_genotype_comparison gtc
	INNER JOIN selected_tumor_pairs stp ON stp.tumor_pair_barcode = gtc.tumor_pair_barcode AND stp.priority = 1
	INNER JOIN selected_genes sg ON sg.chrom = gtc.chrom AND sg.pos = gtc.pos AND sg.alt = gtc.alt
	WHERE 
		(mutect2_call_a OR mutect2_call_b) AND read_depth_a >= 15 AND read_depth_b >= 15
),
variants_by_case AS
(
	SELECT
		case_barcode,
		tumor_pair_barcode,
		tumor_barcode_a,
		tumor_barcode_b,
		bool_and(selected_call_a AND selected_call_b) AS shared,
		bool_or(selected_call_a AND NOT selected_call_b) AS private_a,
		bool_or(NOT selected_call_a AND selected_call_b) AS private_b,
		COUNT(DISTINCT gene_symbol) AS driver_count,
		trim(BOTH ', ' FROM string_agg(CASE WHEN selected_call_b AND NOT selected_call_a THEN '+' || gene_symbol || ', ' ELSE '' END, '')) AS target_a,
		trim(BOTH ', ' FROM string_agg(CASE WHEN selected_call_a AND NOT selected_call_b THEN '-' || gene_symbol || ', ' ELSE '' END, '')) AS target_b
	FROM variants_by_case_and_gene
	WHERE priority = 1
	GROUP BY case_barcode, tumor_pair_barcode, tumor_barcode_a, tumor_barcode_b
),
driver_status AS
(
	SELECT
		case_barcode,
		tumor_pair_barcode,
		tumor_barcode_a,
		tumor_barcode_b,
		driver_count,
		(CASE WHEN shared THEN 'Driver stable'
		 WHEN NOT shared AND private_a AND NOT private_b THEN 'Driver loss'
		 WHEN NOT shared AND private_b AND NOT private_a THEN 'Driver gain'
		 WHEN NOT shared AND private_b AND private_b THEN 'Driver switch' END) driver_status,
		TRIM(BOTH ', ' FROM target_a || ', ' || target_b) AS target
	FROM variants_by_case
)
SELECT * FROM driver_status --driver_status,COUNT(*) FROM driver_status GROUP BY driver_status