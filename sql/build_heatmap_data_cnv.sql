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
		b1.cnv_exclusion = 'allow' AND b2.cnv_exclusion = 'allow' 
),
selected_genes AS
(
	SELECT gene_symbol::varchar(255)
	FROM (VALUES ('CDKN2A'),('PTEN'),('PDGFRA'),('EGFR'),('MDM2'),('MDM4'),('RB1'),('TP53'),('CCNE1'),('CCND1'),('CDK4'),('CDK6'),('MYC'),('MET'),('APC'),('FGFR2')) AS v (gene_symbol)
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
		round(c1.corrected_cn) AS cn_a,
		round(c2.corrected_cn) AS cn_b
	FROM selected_genes_samples sgs
	LEFT JOIN analysis.pairs p1 ON p1.tumor_barcode = sgs.tumor_barcode_a
	LEFT JOIN analysis.pairs p2 ON p2.tumor_barcode = sgs.tumor_barcode_b
	LEFT JOIN analysis.cnv_by_gene c1 ON c1.pair_barcode = p1.pair_barcode AND c1.gene_symbol = sgs.gene_symbol
	LEFT JOIN analysis.cnv_by_gene c2 ON c2.pair_barcode = p2.pair_barcode AND c2.gene_symbol = sgs.gene_symbol
),
cnv_anno AS
(
	SELECT
		*,
		(CASE WHEN cn_a < 1 OR cn_a IS NULL THEN -2 WHEN cn_a < 2 THEN -1 WHEN cn_a = 2 THEN 0 WHEN cn_a > 3 THEN 2 WHEN cn_a > 2 THEN 1 END) cn_class_a,
		(CASE WHEN cn_b < 1 OR cn_b IS NULL THEN -2 WHEN cn_b < 2 THEN -1 WHEN cn_b = 2 THEN 0 WHEN cn_b > 3 THEN 2 WHEN cn_b > 2 THEN 1 END) cn_class_b
	FROM cnv_by_pair_gene
),
cnv_anno_class_change AS
(
	SELECT case_barcode, gene_symbol, cnv_call, cnv_class, cnv_color
	FROM cnv_anno ca
	LEFT JOIN analysis.cn_class_changes cc ON cc.cn_class_a = ca.cn_class_a AND cc.cn_class_b = ca.cn_class_b
)
SELECT * FROM cnv_anno_class_change
/*SELECT cn_class_change,COUNT(*)
FROM cnv_anno_class_change
GROUP BY cn_class_change
ORDER BY 2 DESC*/