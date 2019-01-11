/*
Prepare input for use with dNdS CV
Remove contiguous sites using EXISTS
*/
WITH selected_tumor_pairs AS
(
    SELECT * FROM analysis.silver_set
)

SELECT DISTINCT -- remove duplicate entries
    gt.case_barcode,
    gt.chrom::character varying(2),
    lower(gt.pos) AS pos,
    snvs.ref,
    gt.alt AS mut,
    na.recurrence_evolution AS evolution,
    su.idh_codel_subtype AS subtype,
    (CASE WHEN mutect2_call_a AND mutect2_call_b     THEN 'S'
          WHEN mutect2_call_a AND NOT mutect2_call_b THEN 'P'
          WHEN mutect2_call_b AND NOT mutect2_call_a THEN 'R' END) AS fraction
    /*(CASE WHEN vaf_corrected_call_a AND vaf_corrected_call_b     THEN 'S'
          WHEN vaf_corrected_call_a AND NOT vaf_corrected_call_b THEN 'P'
          WHEN vaf_corrected_call_b AND NOT vaf_corrected_call_a THEN 'R' END) AS fraction*/
    /*(CASE WHEN gt.alt_count_a > 0 AND  gt.alt_count_b > 0    THEN 'S'
          WHEN gt.alt_count_a > 0 AND NOT gt.alt_count_b > 0 THEN 'P'
          WHEN gt.alt_count_b > 0 AND NOT gt.alt_count_a > 0 THEN 'R' END) AS fraction*/
FROM analysis.master_genotype_comparison gt
INNER JOIN selected_tumor_pairs stp ON stp.tumor_pair_barcode = gt.tumor_pair_barcode 
LEFT JOIN analysis.snvs ON snvs.chrom = gt.chrom AND snvs.pos = gt.pos AND snvs.alt = gt.alt
LEFT JOIN analysis.neutrality_tumor_pairs na ON na.tumor_pair_barcode = gt.tumor_pair_barcode
LEFT JOIN clinical.surgeries su ON su.sample_barcode = substring(gt.tumor_pair_barcode from 1 for 15)
WHERE
    stp.priority = 1 AND
    (mutect2_call_a OR mutect2_call_b) AND
    --(gt.alt_count_a + gt.ref_count_a) >= 5 AND 
    --(gt.alt_count_b + gt.ref_count_b) >= 5 AND
    (gt.alt_count_a > 0 OR gt.alt_count_b > 0) AND
    NOT EXISTS (SELECT 1 FROM analysis.snvs snvs WHERE snvs.chrom = gt.chrom AND lower(snvs.pos) = lower(gt.pos) - 1)
    
-- END --