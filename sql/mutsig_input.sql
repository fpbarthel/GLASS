WITH
/*
Define a list of 'selected tumor pairs': select a single primary/recurrent pair per subject/case after removing any aliquots from the blocklist and sorting samples by surgical interval,
so that pairs with the highest surgical interval are prioritized, and pairs with WGS are prioritized over WXS when both are available
*/
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
    sample_type_b <> 'M1' AND
    b1.coverage_exclusion = 'allow' AND b2.coverage_exclusion = 'allow'
    --b1.cnv_exclusion = 'allow' AND b2.cnv_exclusion = 'allow' AND
    --b1.clinical_exclusion = 'allow' AND b2.clinical_exclusion = 'allow'
)
/*
Aggregate counts over tumor pairs.
Performs an INNER JOIN with selected_tumor_pairs from above so that only those pairs are included in the analysis.
Restrict to events with coverage >= 15 in both A and B
*/
SELECT DISTINCT -- remove duplicate entries
  gt.aliquot_barcode,
  gt.chrom::character varying(2),
  lower(gt.pos) AS start_pos,
  upper(gt.pos) AS end_pos,
  snvs.ref,
  gt.alt,
  (CASE WHEN mutect2_call_a AND mutect2_call_b THEN 'S' WHEN mutect2_call_a AND NOT mutect2_call_b THEN 'P' WHEN mutect2_call_b AND NOT mutect2_call_a THEN 'R' END) AS status
FROM analysis.genotypes gt
INNER JOIN selected_tumor_pairs stp ON stp.priority = 1 AND (stp.tumor_barcode_a = gt.aliquot_barcode OR stp.tumor_barcode_b = gt.aliquot_barcode)
LEFT JOIN analysis.snvs ON snvs.chrom = gt.chrom AND snvs.pos = gt.pos AND snvs.alt = gt.alt
WHERE mutect2_call AND snvs.variant_type = 'SNP'