SELECT
      gt.aliquot_barcode,
      gt.chrom || ':' || lower(gt.pos) || '-' || snvs.ref || '/' || gt.alt AS mutation_id,
      ref_count AS ref_counts,
      alt_count AS var_counts,
      (CASE WHEN case_sex = 'male' AND gt.chrom = 'X' THEN 1 ELSE 2 END) AS normal_cn,
      minor_cn,
      major_cn,
      (COUNT(*) OVER window_variant)::integer AS num_aliquots_total,
      (COUNT(CASE WHEN ref_count + alt_count >= 30 THEN 1 END) OVER window_variant)::integer AS num_aliquots_variant_30x,
      (COUNT(CASE WHEN ref_count + alt_count >= 50 THEN 1 END) OVER window_variant)::integer AS num_aliquots_variant_50x,
      (COUNT(CASE WHEN minor_cn IS NOT NULL AND major_cn IS NOT NULL THEN 1 END) OVER window_variant)::integer AS num_aliquots_non_null_cn
FROM analysis.genotypes gt
LEFT JOIN analysis.pairs ps ON ps.tumor_barcode = gt.aliquot_barcode
LEFT JOIN clinical.cases cs ON cs.case_barcode = gt.case_barcode
LEFT JOIN analysis.titan_seg ts ON ts.pair_barcode = ps.pair_barcode AND ts.chrom = gt.chrom AND ts.pos && gt.pos
LEFT JOIN analysis.snvs snvs ON snvs.chrom = gt.chrom AND snvs.pos = gt.pos AND snvs.alt = gt.alt
WHERE
      gt.case_barcode = ? AND 
      aliquot_analysis_type = ? AND 
      (case_sex IS NOT NULL OR gt.chrom <> 'X')
WINDOW window_variant AS (PARTITION BY gt.chrom, gt.pos, gt.alt)
ORDER BY 9,8