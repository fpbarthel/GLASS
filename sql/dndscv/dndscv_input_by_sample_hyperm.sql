/*
Prepare input for use with dNdS CV
Remove contiguous sites using EXISTS
Modified for per sample analysis
*/
WITH
selected_aliquots AS
(
    SELECT gs.tumor_barcode_a AS aliquot_barcode, 'P' AS sample_type FROM analysis.gold_set gs
    INNER JOIN analysis.tumor_clinical_comparison tcc ON tcc.tumor_pair_barcode = gs.tumor_pair_barcode
    WHERE hypermutator_status IS TRUE
    UNION
    SELECT gs.tumor_barcode_b AS aliquot_barcode, 'R' AS sample_type FROM analysis.gold_set gs
    INNER JOIN analysis.tumor_clinical_comparison tcc ON tcc.tumor_pair_barcode = gs.tumor_pair_barcode
    WHERE hypermutator_status IS TRUE
)

SELECT DISTINCT -- remove duplicate entries
    gt.aliquot_barcode,
    (CASE WHEN gt.chrom = 23 THEN 'X' ELSE gt.chrom::varchar(2) END) AS chrom,
    lower(gt.pos) AS pos,
    pa.ref,
    gt.alt AS mut,
    sample_type,
    idh_codel_subtype AS subtype
FROM variants.passgeno gt
INNER JOIN selected_aliquots sa ON sa.aliquot_barcode = gt.aliquot_barcode 
LEFT JOIN variants.passanno pa ON pa.variant_id = gt.variant_id
LEFT JOIN clinical.subtypes su ON su.case_barcode = substring(gt.aliquot_barcode from 1 for 12)
WHERE
    ssm2_pass_call
    
-- END --