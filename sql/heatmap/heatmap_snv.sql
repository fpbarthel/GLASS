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
	SELECT * FROM analysis.gold_set
),
selected_aliquots AS
(
	SELECT tumor_barcode_a AS aliquot_barcode, case_barcode FROM selected_tumor_pairs
	UNION
	SELECT tumor_barcode_b AS aliquot_barcode, case_barcode FROM selected_tumor_pairs
),
selected_genes AS
(
	SELECT DISTINCT sn.gene_symbol, ensembl_gene_id, variant_id, chrom, pos, alt, sn.variant_classification, variant_classification_priority, protein_change
	FROM variants.anno sn
	INNER JOIN ref.driver_genes ds ON ds.gene_symbol = sn.gene_symbol
	INNER JOIN ref.ensembl_gene_mapping gm ON gm.gene_symbol = sn.gene_symbol
	LEFT JOIN variants.variant_classifications vc ON sn.variant_classification = vc.variant_classification
	WHERE
		has_mut IS TRUE AND
		((sn.gene_symbol NOT IN ('TERT','IDH1') AND variant_classification_priority IS NOT NULL) OR
		(sn.gene_symbol = 'TERT' AND sn.variant_classification = 'FIVE_PRIME_FLANK' AND lower(sn.pos) IN (1295228,1295250)) OR
		(sn.gene_symbol = 'IDH1' AND sn.protein_change IN ('p.R132C','p.R132G','p.R132H','p.R132S')))
),
selected_genes_uq AS
(
	SELECT DISTINCT sg.gene_symbol, ensembl_gene_id
	FROM selected_genes sg
),
selected_genes_samples AS
(
	SELECT aliquot_barcode, case_barcode, gene_symbol, ensembl_gene_id, chrom, lower(pos) AS start_pos, upper(pos)-1 AS end_pos, alt
	FROM selected_aliquots, selected_genes
),
selected_genes_samples_uq AS
(
	SELECT DISTINCT aliquot_barcode, case_barcode, gene_symbol, ensembl_gene_id
	FROM selected_genes_samples
),
selected_variants_samples AS
(
	SELECT variant_id, tumor_pair_barcode, variant_classification_priority, protein_change
	FROM selected_genes
	CROSS JOIN selected_tumor_pairs
),
hotspot_coverage AS -- For IDH1 and TERT we don't want genic coverage but coverage at the hotspot sites
(
	SELECT aliquot_barcode, case_barcode, gene_symbol, sum(ad_alt + ad_ref)::decimal / COUNT(*) AS gene_cov
	FROM variants.geno pg
	INNER JOIN variants.passanno pa ON pa.variant_id = pg.variant_id
	LEFT JOIN variants.variant_classifications vc ON vc.variant_classification = pa.variant_classification
	WHERE
		(pa.gene_symbol = 'TERT' AND pa.variant_classification = 'FIVE_PRIME_FLANK' AND lower(pa.pos) IN (1295228,1295250)) OR
		(pa.gene_symbol = 'IDH1' AND pa.protein_change IN ('p.R132C','p.R132G','p.R132H','p.R132S'))
	GROUP BY 1,2,3
	ORDER BY 1,2,3
),
gene_sample_coverage AS -- For all other genes get genic coverage
(
	SELECT sgu.aliquot_barcode, case_barcode, sgu.gene_symbol, gene_coverage::double precision / gene_size AS gene_cov
	FROM selected_genes_samples_uq sgu
	INNER JOIN analysis.gencode_coverage gc ON gc.ensembl_gene_id = sgu.ensembl_gene_id AND gc.aliquot_barcode = sgu.aliquot_barcode
	INNER JOIN ref.ensembl_genes eg ON eg.ensembl_gene_id = sgu.ensembl_gene_id
	WHERE sgu.gene_symbol NOT IN ('IDH1','TERT')

	UNION

	SELECT * FROM hotspot_coverage
),
gene_pair_coverage AS
(
	SELECT
		stp.tumor_pair_barcode,
		stp.case_barcode,
		sg.gene_symbol,
		c1.gene_cov AS cov_a,
		c2.gene_cov AS cov_b
	FROM selected_tumor_pairs stp
	CROSS JOIN selected_genes_uq sg
	LEFT JOIN gene_sample_coverage c1 ON c1.aliquot_barcode = stp.tumor_barcode_a AND c1.gene_symbol = sg.gene_symbol
	LEFT JOIN gene_sample_coverage c2 ON c2.aliquot_barcode = stp.tumor_barcode_b AND c2.gene_symbol = sg.gene_symbol
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
		vc.variant_classification_vep AS variant_classification,
		protein_change,
		alt_count_a > 0 AS selected_call_a, --(alt_count_a::decimal / (alt_count_a+ref_count_a) > 0.05)
		alt_count_b > 0 AS selected_call_b, --(alt_count_b::decimal / (alt_count_b+ref_count_b) > 0.05)
		row_number() OVER (PARTITION BY gtc.gene_symbol, gtc.case_barcode ORDER BY mutect2_call_a::integer + mutect2_call_b::integer = 2, vc.variant_classification_priority, mutect2_call_a::integer + mutect2_call_b::integer DESC, ref_count_a + ref_count_b + alt_count_a + alt_count_b DESC) AS priority
	FROM variants.pgeno gtc
	INNER JOIN selected_variants_samples svs ON svs.tumor_pair_barcode = gtc.tumor_pair_barcode AND svs.variant_id = gtc.variant_id
	INNER JOIN variants.variant_classifications vc ON vc.variant_classification = gtc.variant_classification
	WHERE
		(mutect2_call_a OR mutect2_call_b) AND (alt_count_a+ref_count_a) >= 5 AND (alt_count_b+ref_count_b) >= 5
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
		vcg.protein_change,
		vcg.selected_call_a,
		vcg.selected_call_b,
		gc.cov_a,
		gc.cov_b
	FROM (SELECT * FROM variants_by_case_and_gene WHERE priority = 1) vcg
	RIGHT JOIN gene_pair_coverage gc ON gc.tumor_pair_barcode = vcg.tumor_pair_barcode AND gc.gene_symbol = vcg.gene_symbol
)
--SELECT * --sa.aliquot_barcode, case_barcode, gm.gene_symbol, gene_coverage::double precision / gene_size AS gene_cov
--FROM selected_genes_samples_uq sgu
--INNER JOIN analysis.gencode_coverage gc ON gc.ensembl_gene_id = sgu.ensembl_gene_id AND gc.aliquot_barcode = sgu.aliquot_barcode-- 
--SELECT * from gene_sample_coverage --gene_pair_coverage gc
--SELECT * FROM variants_by_case_and_gene vcg
--RIGHT JOIN gene_pair_coverage gc ON gc.tumor_pair_barcode = vcg.tumor_pair_barcode AND gc.gene_symbol = vcg.gene_symbol
SELECT
	gene_symbol,
	var.case_barcode,
	idh_codel_subtype,
	variant_classification,
	protein_change,
	(CASE
	 WHEN cov_a >= 5 AND cov_b >= 5 AND selected_call_a AND selected_call_b 				THEN 'S'
	 WHEN cov_a >= 5 AND cov_b >= 5 AND selected_call_a AND NOT selected_call_b 			THEN 'P'
	 WHEN cov_a >= 5 AND cov_b >= 5 AND selected_call_b AND NOT selected_call_a 			THEN 'R'
	 ELSE NULL END) AS variant_call,
	(CASE
		WHEN cov_a >= 30 AND cov_b >= 30 THEN 'Coverage >= 30x' 
		--WHEN cov_a >= 20 AND cov_b >= 20 THEN 'Coverage >= 20x' 
		--WHEN cov_a >= 10 AND cov_b >= 5 THEN 'Coverage >= 10x' 
		WHEN cov_a >= 15 AND cov_b >= 15 THEN 'Coverage >= 15x' 
		WHEN cov_a >= 5 AND cov_b >= 5 THEN 'Coverage >= 5x' 
		ELSE 'Coverage < 5x' END) AS covered
FROM squared_variants var
LEFT JOIN clinical.subtypes st ON st.case_barcode = var.case_barcode
--LEFT JOIN selected_tumor_pairs stp ON stp.case_barcode = var.case_barcode AND priority = 1
--LEFT JOIN clinical.surgeries su ON su.sample_barcode = substring(tumor_barcode_a from 1 for 15)