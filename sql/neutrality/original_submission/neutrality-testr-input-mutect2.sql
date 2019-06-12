/*
Prepare input for use with Neutrality Testr
*/
SELECT
    gtc.case_barcode,
    gtc.tumor_barcode_a,
    gtc.tumor_barcode_b,
    gtc.chrom,
    gtc.pos,
    gtc.alt,
    gtc.ref_count_a,
    gtc.ref_count_b,
    gtc.alt_count_a,
    gtc.alt_count_b,
    gtc.ref_count_a + gtc.ref_count_b AS ref_count_ab,
    gtc.alt_count_a + gtc.alt_count_b AS alt_count_ab,
    ROUND(gtc.alt_count_a::decimal / (gtc.alt_count_a + gtc.ref_count_a),4) AS vaf_a,
    ROUND(gtc.alt_count_b::decimal / (gtc.alt_count_b + gtc.ref_count_b),4) AS vaf_b,
    ROUND((gtc.alt_count_a::decimal + gtc.alt_count_b::decimal) / (gtc.alt_count_a + gtc.alt_count_b + gtc.ref_count_a + gtc.ref_count_b),4) AS vaf_ab,
    (CASE WHEN mutect2_call_a AND mutect2_call_b THEN 'S' WHEN mutect2_call_a AND NOT mutect2_call_b THEN 'P' WHEN mutect2_call_b AND NOT mutect2_call_a THEN 'R' END) AS status
FROM analysis.master_genotype_comparison gtc
LEFT JOIN analysis.snvs snvs ON snvs.chrom = gtc.chrom AND snvs.pos = gtc.pos AND snvs.alt = gtc.alt
WHERE 
    (mutect2_call_a OR mutect2_call_b) AND 
    (gtc.alt_count_a + gtc.ref_count_a) >= 30 AND 
    (gtc.alt_count_b + gtc.ref_count_b) >= 30  AND
    (gtc.alt_count_a > 0 OR gtc.alt_count_b > 0)

-- END --