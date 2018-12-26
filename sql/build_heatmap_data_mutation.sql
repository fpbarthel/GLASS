/*
Build a squared gene x subject matrix for heatmap plotting
1. Get selected tumor pairs as described elsewhere (`selected_tumor_pairs`)
	- reduce this list to a list of aliquots, combining (a) and (b) into a single column (`selected_aliquots`). This list is only used for coverage
2. Select a list of genes and specific variants
	- Take genes deemed significant using the dNdS-CV method
	- Manually adding a genes that are not significant by dNdS but are known glioma genes (handpicked from Anzhela's list)
	- Filter the list of variants by subsetting hotspot variants only for those genes where they are known
3. Take the cartesian product between the genes and aliquots table, generating a table with (genes x aliquots) numbers of rows (`selected_genes_samples`)
	- Get coverage statistics from all these loci across all aliquots and compute median for each gene and aliquot (`gene_sample_coverage`)
	- Use the tumor pairs list from step 1 to align tumor pairs to genes and coverage statistics for sample (a) and (b) from each pair (`gene_pair_coverage`)
4. Take variants from each tumor_barcode, subsetting selected tumor pairs from step 1 and by variants from step 2 (`variants_by_case_and_gene`)
	- For each pair, there exists an optimal first tumor (a) and subsequent tumor sample (b)
	- Aggregate over variants by case_barcode and gene_symbol, selecting only the top variant for each subject/gene combination
	- Restrict to events with coverage >= 5 in both (a) and (b)
	- Variants for each tumor pair/gene combination are ordered according to variant_classification_priority (see new table analysis.variant_classifications) and whether or not the mutation was called in a/b and finally based on read_depth
	- Manually calling all variants based on a 0.05 VAF threshold, since we already applied stringent filtering prior to this step (eg. only variants with a mutect call in either (a) or (b) are included in the data at this point)
5. Right join the limited list of variants from step 4 with the extensive coverage statistics from step 3, generating a squared table with all genes and all subjects (`squared_variants`)
6. Perform a subject level call of whether the variant was called in sample A only (initial only), B only (recurrence only), or both (shared), or in neither mark as wild-type.
	- If there is not enough coverage in either A or B, mark the gene/subject combination as NULL to indicate insufficient coverage.
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
	INNER JOIN analysis.dndscv_gene ds ON ds.gene_symbol = sn.gene_symbol AND (ds.qglobal_cv < 0.10 OR ds.gene_symbol IN ('TERT','IDH2','NOTCH1','PDGFRA','PIK3CG','BRAF','H3F3A'))
	LEFT JOIN analysis.variant_classifications vc ON sn.variant_classification = vc.variant_classification
	WHERE
		(sn.gene_symbol NOT IN ('TERT','IDH1','IDH2','BRAF','H3F3A') AND variant_classification_priority IS NOT NULL) OR 
		(sn.gene_symbol = 'TERT' AND sn.variant_classification = '5''Flank' AND lower(sn.pos) IN (1295228,1295250)) OR
		(sn.gene_symbol = 'IDH1' AND sn.hgvs_p IN ('p.R132C','p.R132G','p.R132H','p.R132S')) OR
		(sn.gene_symbol = 'IDH2' AND sn.hgvs_p = 'p.R172K') OR
		(sn.gene_symbol = 'BRAF' AND sn.hgvs_p = 'p.V600E') OR
		(sn.gene_symbol = 'H3F3A' AND sn.hgvs_p = 'p.G35R')
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
		gtc.chrom,
		gtc.pos,
		gtc.alt,
		gtc.variant_classification,
		sg.hgvs_p,
		(alt_count_a::decimal / (alt_count_a+ref_count_a) > 0.05) AS selected_call_a,
		(alt_count_b::decimal / (alt_count_b+ref_count_b) > 0.05) AS selected_call_b,
		row_number() OVER (PARTITION BY gtc.gene_symbol, gtc.case_barcode ORDER BY mutect2_call_a::integer + mutect2_call_b::integer = 2, variant_classification_priority, mutect2_call_a::integer + mutect2_call_b::integer DESC, read_depth_a + read_depth_b DESC) AS priority
	FROM analysis.master_genotype_comparison gtc
	INNER JOIN selected_tumor_pairs stp ON stp.tumor_pair_barcode = gtc.tumor_pair_barcode AND stp.priority = 1
	INNER JOIN selected_genes sg ON sg.chrom = gtc.chrom AND sg.pos = gtc.pos AND sg.alt = gtc.alt
	WHERE
		(mutect2_call_a OR mutect2_call_b) AND read_depth_a >= 5 AND read_depth_b >= 5
),
squared_variants AS
(
	SELECT
		gc.gene_symbol,
		gc.case_barcode,
		vcg.chrom,
		vcg.pos,
		vcg.alt,
		vcg.variant_classification,
		vcg.hgvs_p,
		vcg.selected_call_a,
		vcg.selected_call_b,
		gc.cov_a,
		gc.cov_b  
	FROM (SELECT * FROM variants_by_case_and_gene WHERE priority = 1) vcg
	RIGHT JOIN gene_pair_coverage gc ON gc.tumor_pair_barcode = vcg.tumor_pair_barcode AND gc.gene_symbol = vcg.gene_symbol
)
SELECT
	gene_symbol,
	var.case_barcode,
	variant_classification,
	hgvs_p,
	(CASE
	 WHEN cov_a >= 5 AND cov_b >= 5 AND selected_call_a AND selected_call_b 				THEN 'Shared'
	 WHEN cov_a >= 5 AND cov_b >= 5 AND selected_call_a AND NOT selected_call_b 			THEN 'Shed'
	 WHEN cov_a >= 5 AND cov_b >= 5 AND selected_call_b AND NOT selected_call_a 			THEN 'Acquired'
	 ELSE NULL END) AS variant_call,
	(CASE WHEN cov_a >= 5 AND cov_b >= 5 THEN 'Coverage >= 5' ELSE 'Coverage < 5' END) AS covered--,
	--su.idh_codel_subtype 
FROM squared_variants var
--LEFT JOIN selected_tumor_pairs stp ON stp.case_barcode = var.case_barcode AND priority = 1
--LEFT JOIN clinical.surgeries su ON su.sample_barcode = substring(tumor_barcode_a from 1 for 15)