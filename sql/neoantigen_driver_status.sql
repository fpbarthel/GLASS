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
neoag_mutation AS
(
	SELECT DISTINCT chrom, pos, alt, mutation AS aa_change
	FROM analysis.neoantigens_by_pair
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
		vc.variant_classification_vep AS variant_classification,
		protein_change,
		pv.mutation AS aa_change,
		pv.netmhcpan_mt_score,
		(alt_count_a) > 0 AS selected_call_a, --(alt_count_a::decimal / (alt_count_a+ref_count_a) > 0.05)
		(alt_count_b) > 0 AS selected_call_b, --(alt_count_b::decimal / (alt_count_b+ref_count_b) > 0.05)
		row_number() OVER (PARTITION BY gtc.gene_symbol, gtc.case_barcode ORDER BY mutect2_call_a::integer + mutect2_call_b::integer = 2, vc.variant_classification_priority, mutect2_call_a::integer + mutect2_call_b::integer DESC, ref_count_a + ref_count_b + alt_count_a + alt_count_b DESC) AS priority
	FROM variants.pgeno gtc
	INNER JOIN selected_variants_samples svs ON svs.tumor_pair_barcode = gtc.tumor_pair_barcode AND svs.variant_id = gtc.variant_id
	INNER JOIN variants.variant_classifications vc ON vc.variant_classification = gtc.variant_classification
	LEFT JOIN analysis.neoantigens_by_pair pv ON pv.tumor_pair_barcode = gtc.tumor_pair_barcode AND pv.chrom = gtc.chrom AND pv.pos = gtc.pos AND pv.alt = gtc.alt
	WHERE
		(mutect2_call_a OR mutect2_call_b) AND (alt_count_a+ref_count_a) >= 15 AND (alt_count_b+ref_count_b) >= 15
),
variants_agg AS
(
	SELECT
		gene_symbol, case_barcode, tumor_pair_barcode, tumor_barcode_a, tumor_barcode_b, chrom, pos, alt, variant_classification, protein_change, aa_change, netmhcpan_mt_score,
		aa_change IS NOT NULL AS is_immunogenic,
		(CASE
		 WHEN selected_call_b AND selected_call_a THEN 'S'
		 WHEN selected_call_a AND NOT selected_call_b THEN 'P'
		 WHEN selected_call_b AND NOT selected_call_a THEN 'R' END) fraction
	FROM variants_by_case_and_gene
	WHERE priority = 1
)
SELECT va.* FROM variants_agg va
