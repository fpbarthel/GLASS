WITH selected_aliquots
AS (
  SELECT
    sa.case_barcode,
    al.aliquot_analysis_type,
    al.aliquot_barcode,
    round(tp.purity::numeric, 2) AS purity,
    (su.case_barcode::text || '-'::text) || al.aliquot_analysis_type::text AS short_name,
    count(*) OVER (PARTITION BY su.case_barcode, al.aliquot_analysis_type) AS num_samples
  FROM biospecimen.aliquots al
  LEFT JOIN biospecimen.samples sa ON sa.sample_barcode = al.sample_barcode
  LEFT JOIN analysis.blocklist bl ON bl.aliquot_barcode = al.aliquot_barcode
  LEFT JOIN clinical.surgeries su ON al.sample_barcode = su.sample_barcode
  LEFT JOIN analysis.pairs pa ON al.aliquot_barcode = pa.tumor_barcode
  LEFT JOIN variants.titan_params tp ON tp.pair_barcode = pa.pair_barcode
  WHERE 
    bl.fingerprint_exclusion = 'allow' AND
    bl.coverage_exclusion = 'allow' AND
    bl.cnv_exclusion = 'allow' AND
    sa.sample_type NOT IN ('NB','NM')
  ORDER BY su.case_barcode, su.surgery_number, al.aliquot_portion
)
SELECT
  case_barcode,
  aliquot_analysis_type,
  aliquot_barcode,
  purity,
  short_name,
  num_samples
FROM selected_aliquots
WHERE num_samples > 1;