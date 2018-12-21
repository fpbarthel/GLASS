library(tidyverse)
library(DBI)
library(parallel)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

q <- "SELECT DISTINCT aliquot_analysis_type, case_barcode, aliquot_barcode FROM analysis.genotypes"
qres <- dbGetQuery(con, q)

sample_list <- lapply(split(qres, qres$aliquot_analysis_type, drop = TRUE), function(x) split(x$aliquot_barcode, x$case_barcode, drop = TRUE))

q <- "SELECT
      gt.aliquot_barcode,
      gt.aliquot_barcode || '__' || case_sex || '__' || gt.chrom || ':' || lower(gt.pos) || '-' || snvs.ref || '/' || gt.alt AS mutation_id,
      ref_count AS ref_counts,
      alt_count AS var_counts,
      (CASE WHEN case_sex = 'male' AND gt.chrom = 'X' THEN 1 ELSE 2 END) AS normal_cn,
      minor_cn,
      major_cn,
      (COUNT(*) OVER window_variant)::integer AS num_aliquots_total,
      (COUNT(CASE WHEN ref_count + alt_count >= 30 THEN 1 END) OVER window_variant)::integer AS num_aliquots_variant_30x,
      (COUNT(CASE WHEN ref_count + alt_count >= 50 THEN 1 END) OVER window_variant)::integer AS num_aliquots_variant_50x,
      gt.chrom AS chr,
      lower(gt.pos) AS \"start\",
      upper(gt.pos) - 1 AS \"end\",
      snvs.ref AS ref_allele,
      gt.alt AS alt_allele
      FROM analysis.genotypes gt
      LEFT JOIN analysis.pairs ps ON ps.tumor_barcode = gt.aliquot_barcode
      LEFT JOIN clinical.cases cs ON cs.case_barcode = gt.case_barcode
      LEFT JOIN analysis.titan_seg ts ON ts.pair_barcode = ps.pair_barcode AND ts.chrom = gt.chrom AND ts.pos && gt.pos
      LEFT JOIN analysis.snvs snvs ON snvs.chrom = gt.chrom AND snvs.pos = gt.pos AND snvs.alt = gt.alt
      WHERE gt.case_barcode = ? AND aliquot_analysis_type = ? AND (case_sex IS NOT NULL OR gt.chrom <> 'X')
      WINDOW window_variant AS (PARTITION BY gt.chrom, gt.pos, gt.alt)
      ORDER BY 9,8"

rs <- dbSendQuery(con, q)

## Loop over analysis types
for(analysis_type in names(sample_list)) {
  message(analysis_type)
  
  ## Loop over cases
  for(case_barcode in names(sample_list[[analysis_type]])) {
    
    message(case_barcode)
    dbBind(rs, list(case_barcode, analysis_type))
    qres <- dbFetch(rs)
    
    ## Loop over aliquots
    for(aliquot_barcode in sample_list[[analysis_type]][[case_barcode]]) {
      message(aliquot_barcode)
      tmp = aliquot_barcode
      
      df <- qres %>% 
        filter(aliquot_barcode == tmp) %>% 
        mutate(has_30x_in_all_subsamples = num_aliquots_variant_30x == num_aliquots_total,
               has_50x_in_all_subsamples = num_aliquots_variant_50x == num_aliquots_total) %>%
        select(-aliquot_barcode, -num_aliquots_variant_30x, -num_aliquots_variant_50x, -num_aliquots_total)
      
      fn = sprintf("results/pyclone/tsv/%s_%s_%s.tsv", analysis_type, case_barcode, aliquot_barcode)
      write.table(df, file = fn, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
  }
}
