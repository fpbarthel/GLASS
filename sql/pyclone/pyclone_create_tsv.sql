WITH
t1 AS (
      SELECT
              gt.case_barcode,
              gt.aliquot_barcode,
              gt.variant_id::integer AS mutation_id,
              ad_ref AS ref_counts,
              ad_alt AS var_counts,
              (CASE WHEN case_sex = 'male' AND gt.chrom = 23 THEN 1 ELSE 2 END) AS normal_cn,
              minor_cn,
              major_cn,
              (COUNT(*) OVER (PARTITION BY gt.case_barcode, gt.variant_id)) AS num_aliquots_variants
      FROM variants.passgeno gt
      INNER JOIN biospecimen.aliquots al ON al.aliquot_barcode = gt.aliquot_barcode
      INNER JOIN biospecimen.samples sa ON sa.sample_barcode = al.sample_barcode
    INNER JOIN analysis.blocklist bl ON bl.aliquot_barcode = al.aliquot_barcode
      INNER JOIN analysis.pairs ps ON ps.tumor_barcode = gt.aliquot_barcode
      INNER JOIN clinical.cases cs ON cs.case_barcode = gt.case_barcode
      INNER JOIN variants.titan_seg ts ON ts.pair_barcode = ps.pair_barcode AND ts.chrom = gt.chrom AND ts.pos && gt.pos
      WHERE
              gt.case_barcode = ? AND 
              (case_sex IS NOT NULL OR gt.chrom <> 23) AND
              major_cn > 0 AND
              ad_ref + ad_alt >= 30 AND
              minor_cn IS NOT NULL AND 
              major_cn IS NOT NULL AND
              bl.fingerprint_exclusion = 'allow' AND
              bl.coverage_exclusion = 'allow' AND
              sa.sample_type NOT IN ('NB','NM')
),
t2 AS (
      SELECT case_barcode,(COUNT(DISTINCT aliquot_barcode)) AS num_aliquots FROM t1 GROUP BY 1
),
t3 AS (
      SELECT t1.case_barcode,aliquot_barcode,mutation_id,ref_counts,var_counts,normal_cn,minor_cn,major_cn
      FROM t1
      LEFT JOIN t2 ON t1.case_barcode = t2.case_barcode
      WHERE num_aliquots_variants = num_aliquots AND num_aliquots > 1
),
t4 AS (
      SELECT aliquot_barcode,COUNT(*)
      FROM t3
      GROUP BY 1
)
SELECT * FROM t3-- ORDER BY 1,2 DESC