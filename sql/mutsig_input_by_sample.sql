/*
Aggregate counts over tumor pairs.
Performs an INNER JOIN with selected_tumor_pairs from above so that only those pairs are included in the analysis.
Restrict to events with coverage >= 15 in both A and B
*/
SELECT DISTINCT
  gt.aliquot_barcode,
  gt.chrom::character varying(2),
  lower(gt.pos) AS start_pos,
  upper(gt.pos)-1 AS end_pos,
  snvs.ref,
  gt.alt
FROM analysis.genotypes gt
LEFT JOIN analysis.snvs snvs ON snvs.chrom = gt.chrom AND snvs.pos = gt.pos AND snvs.alt = gt.alt
WHERE
  mutect2_call AND 
  (gt.alt_count + gt.ref_count) >= 15 AND
  snvs.variant_type = 'SNP'

-- END --