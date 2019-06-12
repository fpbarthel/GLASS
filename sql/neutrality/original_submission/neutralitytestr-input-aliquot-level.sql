/*
Prepare aliquot-level input for use with Neutrality Testr
*/
SELECT
    gtc.aliquot_barcode,
    gtc.case_barcode,
    gtc.chrom::character varying(2),
    gtc.pos,
    gtc.alt,
    gtc.ref_count,
    gtc.alt_count,
    ROUND(gtc.alt_count::decimal / (gtc.alt_count + gtc.ref_count),4) AS vaf
FROM analysis.genotypes gtc
WHERE 
    (mutect2_call) AND 
    (gtc.alt_count + gtc.ref_count) >= 30
-- END -- 
