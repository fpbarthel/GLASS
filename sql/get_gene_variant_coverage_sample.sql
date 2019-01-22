/*
Get coverage statistics for each gene and sample
1. Get selected tumor pairs as described elsewhere (`selected_tumor_pairs`)
2. Reduce this list to a list of aliquots, combining (a) and (b) into a single column (`selected_aliquots`)
3. Select a list of genes and specific variants
	- Take genes deemed significant using the dNdS-CV method
	- Manually adding a genes that are not significant by dNdS but are known glioma genes (handpicked from Anzhela's list)
	- Filter the list of variants by subsetting hotspot variants only for those genes where they are known
4. Take the cartesian product between the genes and aliquots table, generating a table with (genes x aliquots) numbers of rows
5. Get coverage statistics from all these loci across all aliquots
6. Compute min, max and median for each gene and aliquot
*/
WITH
selected_tumor_pairs AS
(
	SELECT
		tumor_pair_barcode,
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
	SELECT tumor_barcode_a AS aliquot_barcode FROM selected_tumor_pairs WHERE priority = 1
	UNION
	SELECT tumor_barcode_b AS aliquot_barcode FROM selected_tumor_pairs WHERE priority = 1
),
selected_genes AS
(
	SELECT sn.gene_symbol, chrom, pos, alt
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
	SELECT aliquot_barcode, gene_symbol, chrom, lower(pos) AS start_pos, upper(pos)-1 AS end_pos, alt
	FROM selected_aliquots, selected_genes
)
SELECT
	gene_symbol,
	sg.aliquot_barcode,
	min(read_depth) AS min_read_depth,
	min(alt_count + ref_count) AS min_counts,
	round(median(read_depth)) AS median_read_depth,
	round(median(alt_count + ref_count)) AS median_counts,
	max(read_depth) AS max_read_depth,
	max(alt_count + ref_count) AS max_counts
FROM analysis.full_genotypes fgt
INNER JOIN selected_genes_samples sg ON fgt.aliquot_barcode = sg.aliquot_barcode AND fgt.chrom = sg.chrom AND fgt.start = sg.start_pos AND fgt.end = sg.end_pos AND fgt.alt = sg.alt
GROUP BY gene_symbol, sg.aliquot_barcode	 