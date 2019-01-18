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
	SELECT DISTINCT sn.gene_symbol, chrom, pos, alt, sn.variant_classification, variant_classification_priority, hgvs_p
	FROM analysis.snvs sn
	INNER JOIN analysis.dnds_fraction_sel_cv ds ON ds.gene_symbol = sn.gene_symbol AND (ds.qglobal_cv < 0.10 OR ds.gene_symbol IN ('TERT','IDH2','NOTCH1','PDGFRA','PIK3CG','BRAF','H3F3A'))
	LEFT JOIN analysis.variant_classifications vc ON sn.variant_classification = vc.variant_classification
	WHERE
		(sn.gene_symbol NOT IN ('TERT','IDH1','IDH2','BRAF','H3F3A') AND variant_classification_priority IS NOT NULL) OR 
		(sn.gene_symbol = 'TERT' AND sn.variant_classification = '5''Flank' AND lower(sn.pos) IN (1295228,1295250)) OR
		(sn.gene_symbol = 'IDH1' AND sn.hgvs_p IN ('p.R132C','p.R132G','p.R132H','p.R132S')) OR
		(sn.gene_symbol = 'IDH2' AND sn.hgvs_p = 'p.R172K') OR
		(sn.gene_symbol = 'BRAF' AND sn.hgvs_p = 'p.V600E') OR
		(sn.gene_symbol = 'H3F3A' AND sn.hgvs_p = 'p.G35R')
),
/*
Aggregate counts over tumor pairs and genes
Performs an INNER JOIN with selected_tumor_pairs from above so that only those pairs are included in the analysis.
Restrict to events with coverage >= 15 in both A and B
Variants for each tumor pair/gene combination are ordered according to variant_classification_priority (see new table analysis.variant_classifications) and whether or not the mutation was called in a/b and finally based on read_depth
Note that I am adding `mutect2_call_a` to `mutect2_call_b` (true = 1, false = 0) as to avoid prioritizing mutations in either A or B over the other
The row_number() function asigns a row number to each row within each group of gene_symbol and case_barcode, after ordering by the given parameters
*/
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
/*
The previous step aggregates results by gene and subject using WINDOW functions (the PostgreSQL equivalent of mutate() in dplyr)
Here we are aggregating results by gene using GROUP BY (the equivalent of group_by() in dplyr), by counting the number of subjects
*/
variants_by_gene AS
(
	SELECT
		gene_symbol,
		SUM(CASE WHEN selected_call_a THEN 1 ELSE 0 END)::integer		AS count_a,											-- count of mutations in A
		SUM(CASE WHEN selected_call_b THEN 1 ELSE 0 END)::integer		AS count_b,											-- count of mutations in B
		SUM(CASE WHEN selected_call_a AND selected_call_b THEN 1 ELSE 0 END)::integer 		AS shared, 						-- if both A and B are true then a variant is shared
		SUM(CASE WHEN selected_call_a AND NOT selected_call_b THEN 1 ELSE 0 END)::integer 	AS private_a,					-- if A is true and B is not then it is unique to A
		SUM(CASE WHEN selected_call_b AND NOT selected_call_a THEN 1 ELSE 0 END)::integer 	AS private_b,					-- if B is true and A is not then it is unique to B,
		SUM(CASE WHEN selected_call_b OR selected_call_a THEN 1 ELSE 0 END)::integer 		AS total,						-- if a is true OR B is true, sums to total variants between A and B
		ROUND(SUM(CASE WHEN selected_call_a AND selected_call_b THEN 1 ELSE 0 END) / 
			  SUM(CASE WHEN selected_call_b OR selected_call_a THEN 1 ELSE 0 END)::decimal,2) AS prop_shared,		-- divides `shared` by total to get the proportion shared
		ROUND(SUM(CASE WHEN selected_call_a AND NOT selected_call_b THEN 1 ELSE 0 END) /
			  SUM(CASE WHEN selected_call_b OR selected_call_a THEN 1 ELSE 0 END)::decimal,2) AS prop_private_a,	-- divides `private_a` by total to get the proportion private to A
		ROUND(SUM(CASE WHEN selected_call_b AND NOT selected_call_a THEN 1 ELSE 0 END) /
			  SUM(CASE WHEN selected_call_b OR selected_call_a THEN 1 ELSE 0 END)::decimal,2) AS prop_private_b		-- divides `private_b` by total to get the proportion private to B
	FROM variants_by_case_and_gene
	WHERE priority = 1 		-- select only a single variant per gene/subject combination
	GROUP BY gene_symbol	-- aggregate across genes
	-- HAVING SUM(CASE WHEN selected_call_b OR selected_call_a THEN 1 ELSE 0 END) > 4 -- can remove this line, but basically just removes any rows where `total` < 5
	ORDER BY 7 DESC -- order by column 7 (= total)
	LIMIT 250 -- only show first 250 rows
)
SELECT vg.* --,cds_size,expr,reptime,hic,ROUND((total::decimal/cds_size)*1e3,4) AS normalized_total
FROM variants_by_gene vg
--LEFT JOIN ref.genes rg ON vg.gene_symbol = rg.gene_symbol
--WHERE cds_size > 0 AND (total::decimal/cds_size)*1e3 > 1.5
ORDER BY 7 DESC, 1 DESC