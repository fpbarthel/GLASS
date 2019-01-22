/*
Aggregate counts over tumor pairs.
Performs an INNER JOIN with selected_tumor_pairs from above so that only those pairs are included in the analysis.
Restrict to events with coverage >= 15 in both A and B
*/
SELECT DISTINCT
  gtc.tumor_pair_barcode,
  gtc.chrom::character varying(2),
  lower(gtc.pos) AS start_pos,
  upper(gtc.pos)-1 AS end_pos,
  snvs.ref,
  gtc.alt,
  (CASE WHEN mutect2_call_a AND mutect2_call_b THEN 'S' WHEN mutect2_call_a AND NOT mutect2_call_b THEN 'P' WHEN mutect2_call_b AND NOT mutect2_call_a THEN 'R' END) AS fraction
FROM analysis.master_genotype_comparison gtc
LEFT JOIN analysis.snvs snvs ON snvs.chrom = gtc.chrom AND snvs.pos = gtc.pos AND snvs.alt = gtc.alt
WHERE
  (mutect2_call_a OR mutect2_call_b) AND 
  (gtc.alt_count_a + gtc.ref_count_a) >= 15 AND 
  (gtc.alt_count_b + gtc.ref_count_b) >= 15 AND
  snvs.variant_type = 'SNP'

-- END --