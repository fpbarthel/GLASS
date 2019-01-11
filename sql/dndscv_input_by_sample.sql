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
    snvs.ref,
    gt.alt AS mut,
    sample_type,
	evolution,
	idh_codel_subtype AS subtype
FROM analysis.genotypes gt
INNER JOIN selected_aliquots sa ON sa.aliquot_barcode = gt.aliquot_barcode 
LEFT JOIN analysis.snvs ON snvs.chrom = gt.chrom AND snvs.pos = gt.pos AND snvs.alt = gt.alt
LEFT JOIN analysis.neutrality_aliquots na ON na.aliquot_barcode = gt.aliquot_barcode
LEFT JOIN clinical.surgeries su ON su.sample_barcode = substring(gt.aliquot_barcode from 1 for 15)
WHERE
    mutect2_call AND
    --(gt.alt_count + gt.ref_count) >= 30 AND 
    gt.alt_count > 0 AND
    NOT EXISTS (SELECT 1 FROM analysis.snvs snvs WHERE snvs.chrom = gt.chrom AND lower(snvs.pos) = lower(gt.pos) - 1)
    
-- END --