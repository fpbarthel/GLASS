/*
Prepare input for use with dNdS CV
Use triplets only
*/
WITH
t1 AS (
	SELECT
		s1.case_barcode,
		a1.aliquot_barcode AS tumor_barcode_a,
		a2.aliquot_barcode AS tumor_barcode_b,
		a3.aliquot_barcode AS tumor_barcode_c,
		s1.sample_type AS sample_type_a,
		s2.sample_type AS sample_type_b,
		s3.sample_type AS sample_type_c,
		a1.aliquot_portion AS portion_a,
		a2.aliquot_portion AS portion_b,
		a3.aliquot_portion AS portion_c,
		a1.aliquot_analysis_type AS type_a,
		a2.aliquot_analysis_type AS type_b,
		a3.aliquot_analysis_type AS type_c,
		u1.surgery_number AS surgery_a,
		u2.surgery_number AS surgery_b,
		u3.surgery_number AS surgery_c,
		rank() OVER (PARTITION BY s1.case_barcode ORDER BY s1.case_barcode, u1.surgery_number, u2.surgery_number, u3.surgery_number, a1.aliquot_analysis_type, a2.aliquot_analysis_type, a3.aliquot_analysis_type) AS rnk
	FROM clinical.surgeries u1
	INNER JOIN clinical.surgeries u2 ON u2.case_barcode = u1.case_barcode
	INNER JOIN clinical.surgeries u3 ON u3.case_barcode = u1.case_barcode
	INNER JOIN biospecimen.aliquots a1 ON a1.sample_barcode = u1.sample_barcode
	INNER JOIN biospecimen.aliquots a2 ON a2.sample_barcode = u2.sample_barcode
	INNER JOIN biospecimen.aliquots a3 ON a3.sample_barcode = u3.sample_barcode
	INNER JOIN biospecimen.samples s1 ON a1.sample_barcode = s1.sample_barcode
	INNER JOIN biospecimen.samples s2 ON a2.sample_barcode = s2.sample_barcode
	INNER JOIN biospecimen.samples s3 ON a3.sample_barcode = s3.sample_barcode
	LEFT JOIN analysis.blocklist b1 ON b1.aliquot_barcode = a1.aliquot_barcode
	LEFT JOIN analysis.blocklist b2 ON b2.aliquot_barcode = a2.aliquot_barcode
	LEFT JOIN analysis.blocklist b3 ON b3.aliquot_barcode = a3.aliquot_barcode
	WHERE
		a2.aliquot_analysis_type = a1.aliquot_analysis_type AND
		a3.aliquot_analysis_type = a1.aliquot_analysis_type AND
		u1.surgery_number < u2.surgery_number AND 
		u2.surgery_number < u3.surgery_number AND 
		b1.fingerprint_exclusion = 'allow'::bpchar AND
		b2.fingerprint_exclusion = 'allow'::bpchar AND 
		b1.coverage_exclusion = 'allow'::bpchar AND 
		b2.coverage_exclusion = 'allow'::bpchar AND 
		b3.coverage_exclusion = 'allow'::bpchar AND 
		b3.coverage_exclusion = 'allow'::bpchar
),
--SELECT t1.*, st.idh_codel_subtype
--FROM t1
--INNER JOIN clinical.subtypes st ON st.case_barcode = t1.case_barcode
t2 AS (
	SELECT
	    anno.variant_id,
	    t1.case_barcode,
	    t1.tumor_barcode_a,
	    t1.tumor_barcode_b,
	    t1.tumor_barcode_c,
	    anno.gene_symbol,
	    anno.variant_type,
	    anno.variant_classification,
	    anno.chrom,
	    anno.pos,
	    anno.ref,
	    anno.alt,
	    pg1.ad_ref AS ref_count_a,
	    pg2.ad_ref AS ref_count_b,
	    pg3.ad_ref AS ref_count_c,
	    pg1.ad_alt AS alt_count_a,
	    pg2.ad_alt AS alt_count_b,
	    pg3.ad_alt AS alt_count_c,
	    pg1.af AS af_a,
	    pg2.af AS af_b,
	    pg3.af AS af_c,
	    pg1.ssm2_pass_call AS mutect2_call_a,
	    pg2.ssm2_pass_call AS mutect2_call_b,
	    pg3.ssm2_pass_call AS mutect2_call_c
	FROM t1
	INNER JOIN variants.passgeno pg1 ON t1.tumor_barcode_a = pg1.aliquot_barcode
	INNER JOIN variants.passanno anno ON anno.variant_id = pg1.variant_id
	INNER JOIN variants.passgeno pg2 ON t1.tumor_barcode_b = pg2.aliquot_barcode AND pg2.variant_id = pg1.variant_id
	INNER JOIN variants.passgeno pg3 ON t1.tumor_barcode_c = pg3.aliquot_barcode AND pg3.variant_id = pg1.variant_id
	WHERE t1.rnk = 1
)
SELECT
  t2.case_barcode,
  (CASE WHEN t2.chrom = 23 THEN 'X' ELSE t2.chrom::varchar(2) END) AS chrom,
  lower(t2.pos) AS pos,
  t2.ref,
  t2.alt AS mut,
  st.idh_codel_subtype AS subtype,
    (CASE
    	WHEN mutect2_call_a AND mutect2_call_b AND mutect2_call_c THEN 'ABC'
    	WHEN mutect2_call_a AND mutect2_call_b AND NOT mutect2_call_c THEN 'AB'
    	WHEN mutect2_call_a AND NOT mutect2_call_b AND mutect2_call_c THEN 'AC'
    	WHEN NOT mutect2_call_a AND mutect2_call_b AND mutect2_call_c THEN 'BC'
    	WHEN mutect2_call_a AND NOT mutect2_call_b AND NOT mutect2_call_c THEN 'A'
    	WHEN NOT mutect2_call_a AND mutect2_call_b AND NOT mutect2_call_c THEN 'B'
    	WHEN NOT mutect2_call_a AND NOT mutect2_call_b AND mutect2_call_c THEN 'C'
    END) AS fraction
FROM t2
LEFT JOIN clinical.subtypes st ON st.case_barcode = t2.case_barcode
WHERE
    mutect2_call_a OR mutect2_call_b OR mutect2_call_c
-- END --