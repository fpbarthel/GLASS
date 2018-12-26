library(DBI)
library(tidyverse)
library(ggplot2)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

q <- "WITH
/*
Define a list of 'selected tumor pairs': select a single primary/recurrent pair per subject/case after removing any aliquots from the blocklist and sorting samples by surgical interval,
so that pairs with the highest surgical interval are prioritized, and pairs with WGS are prioritized over WXS when both are available
*/
selected_tumor_pairs AS
(
  SELECT
  tumor_pair_barcode,
  row_number() OVER (PARTITION BY case_barcode ORDER BY surgical_interval_mo DESC, portion_a ASC, portion_b ASC, substring(tumor_pair_barcode from 27 for 3) ASC) AS priority
  FROM analysis.tumor_pairs ps
  LEFT JOIN analysis.blocklist b1 ON b1.aliquot_barcode = ps.tumor_barcode_a
  LEFT JOIN analysis.blocklist b2 ON b2.aliquot_barcode = ps.tumor_barcode_b
  LEFT JOIN analysis.mutation_freq mf1 ON mf1.aliquot_barcode = ps.tumor_barcode_a  -- join with mutation freq to remove hypermutators
  LEFT JOIN analysis.mutation_freq mf2 ON mf2.aliquot_barcode = ps.tumor_barcode_b
  WHERE
  comparison_type = 'longitudinal' AND
  sample_type_b <> 'M1' AND 													-- exclude metastatic samples here because this is outside the scope of our study
  --b1.cnv_exclusion = 'allow' AND b2.cnv_exclusion = 'allow' AND
  b1.coverage_exclusion = 'allow' AND b2.coverage_exclusion = 'allow' AND
  b1.clinical_exclusion = 'allow' AND b2.clinical_exclusion = 'allow' --AND
  --mf1.coverage_adj_mut_freq < 10 AND mf2.coverage_adj_mut_freq < 10			-- filter hypermutators
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
  gtc.chrom,
  gtc.pos,
  gtc.alt,
  gtc.variant_classification,
  snvs.hgvs_p,
  ds.qglobal_cv,
  (alt_count_a::decimal / (alt_count_a+ref_count_a) > 0.05) AS selected_call_a, -- manually calling all variants based on VAF for known genes, rather than using 
  (alt_count_b::decimal / (alt_count_b+ref_count_b) > 0.05) AS selected_call_b,
  row_number() OVER (PARTITION BY gtc.gene_symbol, gtc.case_barcode ORDER BY mutect2_call_a::integer + mutect2_call_b::integer = 1, vc.variant_classification_priority, mutect2_call_a::integer + mutect2_call_b::integer DESC, read_depth_a + read_depth_b DESC) AS priority
  FROM analysis.master_genotype_comparison gtc
  INNER JOIN selected_tumor_pairs stp ON stp.tumor_pair_barcode = gtc.tumor_pair_barcode AND stp.priority = 1
  LEFT JOIN analysis.variant_classifications vc ON gtc.variant_classification = vc.variant_classification
  LEFT JOIN analysis.snvs ON snvs.chrom = gtc.chrom AND snvs.pos = gtc.pos AND snvs.alt = gtc.alt
  INNER JOIN analysis.dndscv_gene ds ON ds.gene_symbol = gtc.gene_symbol AND (ds.qglobal_cv < 0.10 OR ds.gene_symbol IN ('TERT','IDH2','NOTCH1','PDGFRA','PIK3CG','BRAF','H3F3A'))
  WHERE
  (mutect2_call_a OR mutect2_call_b) AND read_depth_a >= 5 AND read_depth_b >= 5
  AND ((gtc.gene_symbol NOT IN ('TERT','IDH1','IDH2','BRAF') AND variant_classification_priority IS NOT NULL) OR 
  (gtc.gene_symbol = 'TERT' AND gtc.variant_classification = '5''Flank' AND lower(gtc.pos) IN (1295228,1295250)) OR
  (gtc.gene_symbol = 'IDH1' AND snvs.hgvs_p IN ('p.R132C','p.R132G','p.R132H','p.R132S')) OR
  (gtc.gene_symbol = 'IDH2' AND snvs.hgvs_p = 'p.R172K') OR
  (gtc.gene_symbol = 'BRAF' AND snvs.hgvs_p = 'p.V600E') OR
  (gtc.gene_symbol = 'H3F3A' AND snvs.hgvs_p = 'p.G35R')) -- this removes any variant types we don't care about, eg. Silent and Intronic mutations, see the analysis.variant_classifications table for more details
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
  --HAVING SUM(CASE WHEN selected_call_b OR selected_call_a THEN 1 ELSE 0 END) > 4 -- can remove this line, but basically just removes any rows where `total` < 5
  ORDER BY 7 DESC -- order by column 7 (= total)
  LIMIT 250 -- only show first 250 rows
  )
  SELECT vg.*,cds_size,expr,reptime,hic,ROUND((total::decimal/cds_size)*1e3,4) AS normalized_total FROM variants_by_gene vg
  LEFT JOIN ref.genes rg ON vg.gene_symbol = rg.gene_symbol
  --WHERE cds_size > 0 AND (total::decimal/cds_size)*1e3 > 1.0
  ORDER BY 4 DESC
  --SELECT vcg.*--, su.idh_codel_subtype
  --FROM variants_by_case_and_gene vcg
  --LEFT JOIN clinical.surgeries su ON vcg.case_barcode = su.case_barcode 
  --WHERE gene_symbol = 'BRAF'
  --ORDER BY pos"
  
qres <- dbGetQuery(con, q)

df <- qres %>% arrange(desc(shared)) %>% mutate(gene_symbol = factor(gene_symbol, levels = gene_symbol))

ggplot(df, aes(x=gene_symbol)) +
  geom_bar(aes(y=count_a), stat="identity") + 
  geom_bar(aes(y=shared), stat="identity", fill = "yellow") + 
  geom_bar(aes(y=-count_b), stat="identity") +
  geom_bar(aes(y=-shared), stat="identity", fill = "yellow") +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw() + 
  labs(y = "Number of Variants\n<-- recurrence -- primary -->")
