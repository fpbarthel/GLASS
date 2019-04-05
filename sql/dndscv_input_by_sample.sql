/*
Prepare input for use with dNdS CV
Remove contiguous sites using EXISTS
Modified for per sample analysis
*/
WITH
selected_aliquots AS
(
    SELECT tumor_barcode_a AS aliquot_barcode, 'P' AS sample_type FROM analysis.silver_set
    UNION
    SELECT tumor_barcode_b AS aliquot_barcode, 'R' AS sample_type FROM analysis.silver_set
)

SELECT DISTINCT -- remove duplicate entries
    gt.aliquot_barcode,
    gt.chrom::character varying(2),
    lower(gt.pos) AS pos,
    anno.ref,
    gt.alt AS mut,
    sample_type,
    --evolution,
    idh_codel_subtype AS subtype
FROM variants.pgeno gt
INNER JOIN selected_aliquots sa ON sa.aliquot_barcode = gt.aliquot_barcode 
LEFT JOIN variants.anno ON anno.chrom = gt.chrom AND anno.pos = gt.pos AND anno.alt = gt.alt
--LEFT JOIN analysis.neutrality_aliquots na ON na.aliquot_barcode = gt.aliquot_barcode
LEFT JOIN clinical.subtypes su ON su.case_barcode = substring(gt.aliquot_barcode from 1 for 12)
WHERE
    ssm2_pass_call AND
    --(gt.alt_count + gt.ref_count) >= 30 AND 
    gt.ad_alt > 0 AND
    NOT EXISTS (SELECT 1 FROM variants.anno WHERE anno.chrom = gt.chrom AND lower(anno.pos) = lower(gt.pos) - 1)
    
-- END --