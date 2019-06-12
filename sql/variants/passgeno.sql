/*
    - - - - - - - - - -
    variants.passgeno
    - - - - - - - - - -
    List of genotypes of PASS variants only
*/
WITH
t1 AS (
  SELECT DISTINCT ON (gt.aliquot_barcode,gt.variant_id)
    gt.aliquot_barcode,
    gt.variant_id,
    gt.case_barcode,
    gt.chrom,
    gt.pos,
    gt.alt,
    ad_ref,
    ad_alt,
    af,
    ssm2_pass_call,
    ssm2_saaf_none,
    ts.major_cn,
    ts.minor_cn,
    ts.clonal_cluster AS titan_cluster,
    ts.cellular_prevalence AS titan_ccf,
    pl.cluster_id AS pyclone_cluster,
    pl.cellular_prevalence AS pyclone_ccf,
    pl.cellular_prevalence_sd AS pyclone_ccf_sd,
    pl.variant_allele_frequency AS pyclone_vaf,
    SUM(ad_alt) OVER (PARTITION BY gt.variant_id, gt.case_barcode) AS sum_alt_case
  FROM variants.info
  INNER JOIN variants.geno gt ON gt.case_barcode = info.case_barcode AND gt.variant_id = info.variant_id
  LEFT JOIN variants.pyclone_loci pl ON pl.variant_id = info.variant_id AND pl.aliquot_barcode = gt.aliquot_barcode
  LEFT JOIN variants.titan_seg ts ON ts.aliquot_barcode = gt.aliquot_barcode AND ts.chrom = gt.chrom AND ts.pos && gt.pos
  WHERE
    (info.chrom = 2 AND info.pos IN ('[209113112,209113112]', '[209113113,209113113]')) OR
    (info.chrom = 5 AND info.pos IN ('[295169,295169]', '[1295228,1295228]', '[1295242,1295242]', '[1295250,1295250]')) OR
    (info.chrom = 15 AND info.pos IN ('[90631837,90631837]', '[90631838,90631838]', '[90631839,90631839]'))
),
t2 AS (
  SELECT DISTINCT ON (gt.aliquot_barcode,gt.variant_id)
    gt.aliquot_barcode,
    gt.variant_id,
    gt.case_barcode,
    gt.chrom,
    gt.pos,
    gt.alt,
    ad_ref,
    ad_alt,
    af,
    ssm2_pass_call,
    ssm2_saaf_none,
    ts.major_cn,
    ts.minor_cn,
    ts.clonal_cluster AS titan_cluster,
    ts.cellular_prevalence AS titan_ccf,
    pl.cluster_id AS pyclone_cluster,
    pl.cellular_prevalence AS pyclone_ccf,
    pl.cellular_prevalence_sd AS pyclone_ccf_sd,
    pl.variant_allele_frequency AS pyclone_vaf
  FROM variants.info
  INNER JOIN variants.geno gt ON gt.case_barcode = info.case_barcode AND gt.variant_id = info.variant_id
  LEFT JOIN variants.pyclone_loci pl ON pl.variant_id = info.variant_id AND pl.aliquot_barcode = gt.aliquot_barcode
  LEFT JOIN variants.titan_seg ts ON ts.aliquot_barcode = gt.aliquot_barcode AND ts.chrom = gt.chrom AND ts.pos && gt.pos
  WHERE
    info.filter = 'PASS'
)

SELECT aliquot_barcode, variant_id, case_barcode, chrom, pos, alt, ad_ref, ad_alt, af, ssm2_pass_call, ssm2_saaf_none, major_cn, minor_cn, titan_cluster, titan_ccf, pyclone_cluster, pyclone_ccf, pyclone_vaf FROM t1 WHERE sum_alt_case > 0
UNION
SELECT aliquot_barcode, variant_id, case_barcode, chrom, pos, alt, ad_ref, ad_alt, af, ssm2_pass_call, ssm2_saaf_none, major_cn, minor_cn, titan_cluster, titan_ccf, pyclone_cluster, pyclone_ccf, pyclone_vaf FROM t2