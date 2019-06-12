/*
Fractionated mutation proportions per gene (and subtype)
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
		sg.protein_change,
		alt_count_a > 0 AS selected_call_a, --(alt_count_a::decimal / (alt_count_a+ref_count_a) > 0.05)
		alt_count_b > 0 AS selected_call_b, --(alt_count_b::decimal / (alt_count_b+ref_count_b) > 0.05)
		row_number() OVER (PARTITION BY gtc.gene_symbol, gtc.case_barcode ORDER BY mutect2_call_a::integer + mutect2_call_b::integer = 2, variant_classification_priority, mutect2_call_a::integer + mutect2_call_b::integer DESC, ref_count_a + ref_count_b + alt_count_a + alt_count_b DESC) AS priority
	FROM variants.pgeno gtc
	INNER JOIN selected_tumor_pairs stp ON stp.tumor_pair_barcode = gtc.tumor_pair_barcode
	INNER JOIN selected_genes sg ON sg.chrom = gtc.chrom AND sg.pos = gtc.pos AND sg.alt = gtc.alt
	WHERE
		(mutect2_call_a OR mutect2_call_b) AND (alt_count_a+ref_count_a) >= 5 AND (alt_count_b+ref_count_b) >= 5
),
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
	GROUP BY 1	-- aggregate across genes
	ORDER BY 7 DESC -- order by column 7 (= total)
),
variants_by_gene_subtype AS
(
	SELECT
		idh_codel_subtype,
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
	FROM variants_by_case_and_gene vg
	LEFT JOIN clinical.subtypes st ON st.case_barcode = vg.case_barcode
	WHERE priority = 1 		-- select only a single variant per gene/subject combination
	GROUP BY 1,2	-- aggregate across genes
	ORDER BY 7 DESC -- order by column 7 (= total)
)
(SELECT 'all' AS idh_codel_subtype, vg.*
FROM variants_by_gene vg
--ORDER BY 8 DESC, 2 DESC
		 
UNION
		 
SELECT vs.*
FROM variants_by_gene_subtype vs)
ORDER BY 1 ASC, 8 DESC