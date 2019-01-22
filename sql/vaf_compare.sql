WITH
variants_by_case_and_gene AS
(
	SELECT
		gtc.gene_symbol,
		gtc.case_barcode,
		gtc.variant_classification,
		sn.hgvs_p,
		ROUND(alt_count_a::decimal / (alt_count_a + ref_count_a),4) AS vaf_a,
		ROUND(alt_count_b::decimal / (alt_count_b + ref_count_b),4) AS vaf_b,
		row_number() OVER (PARTITION BY gtc.gene_symbol, gtc.case_barcode ORDER BY vc.variant_classification_priority, mutect2_call_a::integer + mutect2_call_b::integer DESC, (alt_count_a + ref_count_a) + (alt_count_b + ref_count_b) DESC) AS priority
	FROM analysis.master_genotype_comparison gtc
	INNER JOIN analysis.silver_set stp ON stp.tumor_pair_barcode = gtc.tumor_pair_barcode
	INNER JOIN analysis.dnds_fraction_sel_cv ds ON ds.gene_symbol = gtc.gene_symbol --AND (ds.qglobal_cv < 0.05 OR gtc.gene_symbol IN ('TERT','IDH2','NOTCH1','PDGFRA','PIK3CG','BRAF','H3F3A'))
	LEFT JOIN analysis.variant_classifications vc ON gtc.variant_classification = vc.variant_classification
	INNER JOIN analysis.snvs sn ON sn.chrom = gtc.chrom AND sn.pos = gtc.pos AND sn.alt = gtc.alt
	WHERE
		(mutect2_call_a OR mutect2_call_b) AND
	(ds.qglobal_cv < 0.05 OR ds.gene_symbol IN ('TERT','IDH2','NOTCH1','PDGFRA','PIK3CG','BRAF','H3F3A')) AND
		(alt_count_a + ref_count_a) >= 5 AND (alt_count_b + ref_count_b) >= 5 AND
		(gtc.gene_symbol NOT IN ('TERT','IDH1','IDH2','BRAF','H3F3A') AND variant_classification_priority IS NOT NULL) OR 
		(gtc.gene_symbol = 'TERT' AND gtc.variant_classification = '5''Flank' AND lower(sn.pos) IN (1295228,1295250)) OR
		(gtc.gene_symbol = 'IDH1' AND sn.hgvs_p IN ('p.R132C','p.R132G','p.R132H','p.R132S')) OR
		(gtc.gene_symbol = 'IDH2' AND sn.hgvs_p = 'p.R172K') OR
		(gtc.gene_symbol = 'BRAF' AND sn.hgvs_p = 'p.V600E') OR
		(gtc.gene_symbol = 'H3F3A' AND sn.hgvs_p = 'p.G35R')
)
SELECT gene_symbol, case_barcode, variant_classification, hgvs_p, vaf_a, vaf_b
FROM variants_by_case_and_gene vg
WHERE priority = 1