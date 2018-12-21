
#library(devtools); install_github("im3sanger/dndscv")
library(dndscv)
library(tidyverse)
library(DBI)
library(ggthemes)
library(ggplot2)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

q <- "WITH selected_tumor_pairs AS
      (
        SELECT
          tumor_pair_barcode,
          tumor_barcode_a,
          tumor_barcode_b,
          row_number() OVER (PARTITION BY case_barcode ORDER BY surgical_interval_mo DESC, portion_a ASC, portion_b ASC, substring(tumor_pair_barcode from 27 for 3) ASC) AS priority
        FROM analysis.tumor_pairs ps
        LEFT JOIN analysis.blocklist b1 ON b1.aliquot_barcode = ps.tumor_barcode_a
        LEFT JOIN analysis.blocklist b2 ON b2.aliquot_barcode = ps.tumor_barcode_b
        LEFT JOIN analysis.mutation_freq mf1 ON mf1.aliquot_barcode = ps.tumor_barcode_a  -- join with mutation freq to remove hypermutators
        LEFT JOIN analysis.mutation_freq mf2 ON mf2.aliquot_barcode = ps.tumor_barcode_b
        WHERE
          comparison_type = 'longitudinal' AND
          sample_type_b <> 'M1' 													-- exclude metastatic samples here because this is outside the scope of our study
          AND b1.coverage_exclusion = 'allow' AND b2.coverage_exclusion = 'allow' %s
          --b1.cnv_exclusion = 'allow' AND b2.cnv_exclusion = 'allow' AND
          --b1.clinical_exclusion = 'allow' AND b2.clinical_exclusion = 'allow' 
      )
      SELECT DISTINCT -- remove duplicate entries
        gt.case_barcode,
        gt.chrom::character varying(2),
        lower(gt.pos) AS pos,
        snvs.ref,
        gt.alt AS mut,
        idh_codel_subtype
      FROM analysis.genotypes gt
      INNER JOIN selected_tumor_pairs stp ON stp.priority = 1 AND %s
      LEFT JOIN analysis.snvs ON snvs.chrom = gt.chrom AND snvs.pos = gt.pos AND snvs.alt = gt.alt
      LEFT JOIN analysis.variant_classifications vc ON gt.variant_classification = vc.variant_classification
      LEFT JOIN clinical.surgeries su ON su.sample_barcode = substring(gt.aliquot_barcode from 1 for 15)
      WHERE mutect2_call"

qres_all <- dbGetQuery(con, sprintf(q, "", "(stp.tumor_barcode_a = gt.aliquot_barcode OR stp.tumor_barcode_b = gt.aliquot_barcode)"))
dnds_all = dndscv(qres_all %>% select(-idh_codel_subtype), refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = 500)

result_list <- lapply(unique(qres$idh_codel_subtype), function(qres_subtype) {
  message("Computing dNdS for ", qres_subtype)
  qres_subset = qres_all %>% filter(idh_codel_subtype == qres_subtype) %>% dplyr::select(-idh_codel_subtype)
  dnds_subset = dndscv(qres_subset, refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = 500)
  cbind(subtype = qres_subtype, dnds_subset$globaldnds)
})

res <- data.table::rbindlist(result_list)

dbWriteTable(con, Id(schema="analysis",table="dndscv_global_by_subtype"), as.data.frame(res))


print(dnds_all$globaldnds)

sel_cv = dnds_all$sel_cv
print(head(sel_cv, 20), digits = 3)
print(sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")])

sitedn_all <- sitednds(dnds_all)

recursites = sitedn_all$recursites
print(head(recursites, 20), digits = 3)
print(recursites[recursites$qval < 0.05, c("gene","aachange","impact","freq","qval")])

## Push to database
dbWriteTable(con, Id(schema="analysis",table="dndscv_global"), dnds_all$globaldnds)
dbWriteTable(con, Id(schema="analysis",table="dndscv_gene"), dnds_all$sel_cv)
dbWriteTable(con, Id(schema="analysis",table="dndscv_sites"), sitedn_all$recursites)

## see if thers differenceds in genes per subtype
qres_subset = qres_all %>% filter(idh_codel_subtype == 'IDHmut_codel') %>% dplyr::select(-idh_codel_subtype)
dnds_subset = dndscv(qres_subset, refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = 500)
sel_cv = dnds_subset$sel_cv
print(head(sel_cv, 20), digits = 3)
print(sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")])

qres_subset = qres_all %>% filter(idh_codel_subtype == 'IDHmut_noncodel') %>% dplyr::select(-idh_codel_subtype)
dnds_subset = dndscv(qres_subset, refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = 500)
sel_cv = dnds_subset$sel_cv
print(head(sel_cv, 20), digits = 3)
print(sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")])

qres_subset = qres_all %>% filter(idh_codel_subtype == 'IDHwt_noncodel') %>% dplyr::select(-idh_codel_subtype)
dnds_subset = dndscv(qres_subset, refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = 500)
sel_cv = dnds_subset$sel_cv
print(head(sel_cv, 20), digits = 3)
print(sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")])

###
# Seperate hypermutants
###

qres_all_nonhypermut <- dbGetQuery(con, sprintf(q, "AND mf1.coverage_adj_mut_freq < 10 AND mf2.coverage_adj_mut_freq < 10", "(stp.tumor_barcode_a = gt.aliquot_barcode OR stp.tumor_barcode_b = gt.aliquot_barcode)"))
qres_all_yeshypermut <- dbGetQuery(con, sprintf(q, "AND mf1.coverage_adj_mut_freq >= 10 AND mf2.coverage_adj_mut_freq >= 10", "(stp.tumor_barcode_a = gt.aliquot_barcode OR stp.tumor_barcode_b = gt.aliquot_barcode)"))

dnds_all_n = dndscv(qres_all_nonhypermut, refdb = "hg19", outmats = TRUE)
dnds_all_y = dndscv(qres_all_yeshypermut, refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = Inf)

print(dnds_all_n$globaldnds)
print(dnds_all_y$globaldnds)

###
# RUN dNdS CV seperately for primary and recurrences
###

qres_pri <- dbGetQuery(con, sprintf(q, "", "(stp.tumor_barcode_a = gt.aliquot_barcode)"))
qres_rec <- dbGetQuery(con, sprintf(q, "", "(stp.tumor_barcode_b = gt.aliquot_barcode)"))

dnds_pri = dndscv(qres_pri, refdb = "hg19", outmats = TRUE)
dnds_rec = dndscv(qres_rec, refdb = "hg19", outmats = TRUE)

print(dnds_pri$globaldnds)
print(dnds_rec$globaldnds)

########################################################################################################################
## Seperate private and shared mutations

q <- "WITH
      /*
      Define a list of 'selected tumor pairs': select a single primary/recurrent pair per subject/case after removing any aliquots from the blocklist and sorting samples by surgical interval,
      so that pairs with the highest surgical interval are prioritized, and pairs with WGS are prioritized over WXS when both are available
      */
      selected_tumor_pairs AS
      (
        SELECT
          tumor_pair_barcode,
          row_number() OVER (PARTITION BY case_barcode ORDER BY surgical_interval_mo DESC, portion_a ASC, portion_b ASC, substring(tumor_pair_barcode from 27 for 3) ASC) AS priority
        FROM analysis.tumor_pairs ps
        LEFT JOIN analysis.blocklist b1 ON b1.aliquot_barcode = ps.tumor_barcode_a
        LEFT JOIN analysis.blocklist b2 ON b2.aliquot_barcode = ps.tumor_barcode_b
        WHERE
          comparison_type = 'longitudinal' AND
          sample_type_b <> 'M1' AND
          --b1.cnv_exclusion = 'allow' AND b2.cnv_exclusion = 'allow' AND
          b1.coverage_exclusion = 'allow' AND b2.coverage_exclusion = 'allow' --AND
          --b1.clinical_exclusion = 'allow' AND b2.clinical_exclusion = 'allow'
      )
      /*
      Aggregate counts over tumor pairs.
      Performs an INNER JOIN with selected_tumor_pairs from above so that only those pairs are included in the analysis.
      Restrict to events with coverage >= 15 in both A and B
      */
      SELECT
        gtc.case_barcode,
        gtc.chrom::character varying(2),
        lower(gtc.pos) AS pos,
        snvs.ref,
        gtc.alt AS mut,
        (CASE WHEN mutect2_call_a AND mutect2_call_b THEN 'shared' WHEN mutect2_call_a AND NOT mutect2_call_b THEN 'primary' WHEN mutect2_call_b AND NOT mutect2_call_a THEN 'recurrent' END) AS status
      FROM analysis.master_genotype_comparison gtc
      INNER JOIN selected_tumor_pairs stp ON stp.tumor_pair_barcode = gtc.tumor_pair_barcode AND stp.priority = 1
      LEFT JOIN analysis.snvs snvs ON snvs.chrom = gtc.chrom AND snvs.pos = gtc.pos AND snvs.alt = gtc.alt
      WHERE (mutect2_call_a OR mutect2_call_b) AND read_depth_a >= 15 AND read_depth_b >= 15;"

time <- system.time(qres <- dbGetQuery(con, q))

qres_shared = qres %>% filter(status == 'shared') %>% dplyr::select(-status)
qres_primary = qres %>% filter(status == 'primary') %>% dplyr::select(-status)
qres_recurrent = qres %>% filter(status == 'recurrent') %>% dplyr::select(-status)

dnds_shared = dndscv(qres_shared, refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = 500)
dnds_primary = dndscv(qres_primary, refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = 500)
dnds_recurrent = dndscv(qres_recurrent, refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = 500)

print(dnds_shared$globaldnds)
print(dnds_primary$globaldnds)
print(dnds_recurrent$globaldnds)

df <- rbind(cbind(status = "shared", dnds_shared$globaldnds),
            cbind(status = "primary", dnds_primary$globaldnds),
            cbind(status = "recurrent", dnds_recurrent$globaldnds))
rownames(df) <- NULL
dbWriteTable(con, Id(schema="analysis",table="dndscv_global_private_vs_shared"), df)

sel_cv = dnds_shared$sel_cv
print(head(sel_cv, 20), digits = 3)
print(sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")])

sel_cv = dnds_primary$sel_cv
print(head(sel_cv, 20), digits = 3)
print(sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")])

sel_cv = dnds_recurrent$sel_cv
print(head(sel_cv, 20), digits = 3)
print(sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")])

########################################################################################################################
## By idh-codel status

q <- "WITH
      /*
      Define a list of 'selected tumor pairs': select a single primary/recurrent pair per subject/case after removing any aliquots from the blocklist and sorting samples by surgical interval,
      so that pairs with the highest surgical interval are prioritized, and pairs with WGS are prioritized over WXS when both are available
      */
      selected_tumor_pairs AS
      (
        SELECT
          tumor_pair_barcode,
          row_number() OVER (PARTITION BY case_barcode ORDER BY surgical_interval_mo DESC, portion_a ASC, portion_b ASC, substring(tumor_pair_barcode from 27 for 3) ASC) AS priority
        FROM analysis.tumor_pairs ps
        LEFT JOIN analysis.blocklist b1 ON b1.aliquot_barcode = ps.tumor_barcode_a
        LEFT JOIN analysis.blocklist b2 ON b2.aliquot_barcode = ps.tumor_barcode_b
        WHERE
          comparison_type = 'longitudinal' AND
          sample_type_b <> 'M1' AND
          --b1.cnv_exclusion = 'allow' AND b2.cnv_exclusion = 'allow' AND
          b1.coverage_exclusion = 'allow' AND b2.coverage_exclusion = 'allow' --AND
          --b1.clinical_exclusion = 'allow' AND b2.clinical_exclusion = 'allow'
      )
      /*
      Aggregate counts over tumor pairs.
      Performs an INNER JOIN with selected_tumor_pairs from above so that only those pairs are included in the analysis.
      Restrict to events with coverage >= 15 in both A and B
      */
      SELECT
        gtc.case_barcode,
        gtc.chrom::character varying(2),
        lower(gtc.pos) AS pos,
        snvs.ref,
        gtc.alt AS mut,
        (CASE WHEN mutect2_call_a AND mutect2_call_b THEN 'shared' WHEN mutect2_call_a AND NOT mutect2_call_b THEN 'primary' WHEN mutect2_call_b AND NOT mutect2_call_a THEN 'recurrent' END) AS status,
        su.idh_codel_subtype
      FROM analysis.master_genotype_comparison gtc
      INNER JOIN selected_tumor_pairs stp ON stp.tumor_pair_barcode = gtc.tumor_pair_barcode AND stp.priority = 1
      LEFT JOIN analysis.snvs snvs ON snvs.chrom = gtc.chrom AND snvs.pos = gtc.pos AND snvs.alt = gtc.alt
      LEFT JOIN clinical.surgeries su ON su.sample_barcode = substring(gtc.tumor_pair_barcode from 1 for 15)
      WHERE (mutect2_call_a OR mutect2_call_b) AND read_depth_a >= 15 AND read_depth_b >= 15;"

time <- system.time(qres <- dbGetQuery(con, q))

result_list <- lapply(unique(qres$status), function(qres_status) {
  lapply(unique(qres$idh_codel_subtype), function(qres_subtype) {
    message("Computing dNdS for ", qres_status, " and ", qres_subtype)
    qres_subset = qres %>% filter(status == qres_status, idh_codel_subtype == qres_subtype) %>% dplyr::select(-status, -idh_codel_subtype)
    dnds_subset = dndscv(qres_subset, refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = 500)
    cbind(status = qres_status, subtype = qres_subtype, dnds_subset$globaldnds)
  })
})

res <- data.table::rbindlist(lapply(result_list, data.table::rbindlist))
dbWriteTable(con, Id(schema="analysis",table="dndscv_global_private_vs_shared_by_subtype"), as.data.frame(res))

qres_shared = qres %>% filter(status == 'shared') %>% dplyr::select(-status)
qres_primary = qres %>% filter(status == 'primary') %>% dplyr::select(-status)
qres_recurrent = qres %>% filter(status == 'recurrent') %>% dplyr::select(-status)

dnds_shared = dndscv(qres_shared, refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = 500)
dnds_primary = dndscv(qres_primary, refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = 500)
dnds_recurrent = dndscv(qres_recurrent, refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = 500)


print(dnds_primary$globaldnds)
print(dnds_recurrent$globaldnds)

df <- rbind(cbind(status = "shared", dnds_shared$globaldnds),
            cbind(status = "primary", dnds_primary$globaldnds),
            cbind(status = "recurrent", dnds_recurrent$globaldnds))
rownames(df) <- NULL
dbWriteTable(con, Id(schema="analysis",table="dndscv_global_private_vs_shared"), as.data.frame(df))

sel_cv = dnds_shared$sel_cv
print(head(sel_cv, 20), digits = 3)
print(sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")])

sel_cv = dnds_primary$sel_cv
print(head(sel_cv, 20), digits = 3)
print(sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")])

sel_cv = dnds_recurrent$sel_cv
print(head(sel_cv, 20), digits = 3)
print(sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")])

## END ##

## END ##