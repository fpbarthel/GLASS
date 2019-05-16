/*
Prepare input for use with neutralitytestr using SNVs in non-altered regions.
Also, retain information about clonality AND variant_classification and variant_type.
*/

WITH t1 AS (SELECT 
    pg.tumor_pair_barcode,
    pg.case_barcode,
    pg.tumor_barcode_a,
    pg.tumor_barcode_b,
    pg.chrom, 
    pg.pos,
    pg.variant_id,
    pg.variant_type,
    pg.variant_classification,
    pg.mutect2_call_a,
    pg.mutect2_call_b, 
    pl1.cellular_prevalence AS cellular_prevalence_a, 
    pl1.variant_allele_frequency AS variant_allele_frequency_a, 
    (CASE WHEN pl1.cellular_prevalence >= 0.5 THEN 'C' WHEN pl1.cellular_prevalence < 0.5 THEN 'S' END) AS clonality_a,
    pl2.cellular_prevalence AS cellular_prevalence_b, 
    pl2.variant_allele_frequency AS variant_allele_frequency_b, 
    (CASE WHEN pl2.cellular_prevalence >= 0.5 THEN 'C' WHEN pl2.cellular_prevalence < 0.5 THEN 'S' END) AS clonality_b,
    (CASE WHEN mutect2_call_a AND mutect2_call_b THEN 'S' 
    WHEN mutect2_call_a AND NOT mutect2_call_b THEN 'P' 
    WHEN mutect2_call_b AND NOT mutect2_call_a THEN 'R' END) AS fraction
FROM variants.pgeno pg
LEFT JOIN variants.pyclone_loci pl1 ON pl1.variant_id = pg.variant_id AND pl1.aliquot_barcode = pg.tumor_barcode_a 
LEFT JOIN variants.pyclone_loci pl2 ON pl2.variant_id = pg.variant_id AND pl2.aliquot_barcode= pg.tumor_barcode_b
INNER JOIN analysis.silver_set ss ON pg.tumor_pair_barcode = ss.tumor_pair_barcode
WHERE pl1.cellular_prevalence IS NOT NULL) 

SELECT
    t1.tumor_pair_barcode,
    t1.case_barcode,
    t1.tumor_barcode_a,
    t1.tumor_barcode_b,
    t1.chrom, 
    t1.pos,
    t1.variant_id,
    t1.variant_type,
    t1.variant_classification,
    t1.mutect2_call_a,
    t1.mutect2_call_b,
    t1.cellular_prevalence_a,
    t1.cellular_prevalence_b,
    t1.variant_allele_frequency_a,
    t1.variant_allele_frequency_b,
    t1.clonality_a,
    t1.clonality_b,
    t1.fraction,
    gs1.cnv_call AS cnv_call_a, 
    gs2.cnv_call AS cnv_call_b
FROM t1 
LEFT JOIN variants.gatk_seg gs1 ON gs1.aliquot_barcode = t1.tumor_barcode_a AND gs1.chrom = t1.chrom AND gs1.pos && t1.pos 
LEFT JOIN variants.gatk_seg gs2 ON gs2.aliquot_barcode = t1.tumor_barcode_b AND gs2.chrom = t1.chrom AND gs2.pos && t1.pos 
WHERE t1.fraction IS NOT NULL

-- END --