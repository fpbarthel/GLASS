WITH
/*
Define a list of 'selected tumor pairs': select a single primary/recurrent pair per subject/case after removing any aliquots from the blocklist and sorting samples by surgical interval,
so that pairs with the highest surgical interval are prioritized, and pairs with WGS are prioritized over WXS when both are available
*/
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
		b1.cnv_exclusion <> 'block' AND b2.cnv_exclusion <> 'block' 
),
selected_genes AS
(
	SELECT gene_symbol::varchar(255)
	FROM (VALUES ('CDKN2A'),('PTEN'),('PDGFRA'),('EGFR'),('MDM2'),('MDM4'),('RB1'),('TP53'),('CDK4'),('CDK6'),('MYC'),('MET')) AS v (gene_symbol)
),
selected_genes_samples AS
(
	SELECT * FROM selected_tumor_pairs, selected_genes WHERE priority = 1
),
cnv_by_pair_gene AS
(
	SELECT
		sgs.case_barcode,
		sgs.gene_symbol,
		c1.cn_call AS cn_a,
		c2.cn_call AS cn_b
	FROM selected_genes_samples sgs
	LEFT JOIN analysis.cnv_by_gene_gatk c1 ON c1.aliquot_barcode = sgs.tumor_barcode_a AND c1.gene_symbol = sgs.gene_symbol
	LEFT JOIN analysis.cnv_by_gene_gatk c2 ON c2.aliquot_barcode = sgs.tumor_barcode_b AND c2.gene_symbol = sgs.gene_symbol
)
SELECT case_barcode, gene_symbol, cnv_retention, cnv_class
FROM cnv_by_pair_gene cn
LEFT JOIN analysis.gatk_seg_retention_states gs ON gs.cn_a = cn.cn_a AND gs.cn_b = cn.cn_b