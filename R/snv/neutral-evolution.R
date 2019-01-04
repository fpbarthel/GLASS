#######################################################
# Test neutral evolution using package provided by Marc Williams (Nature Genetics paper).
# Updated: 2019.01.02
# Author: Kevin J.
#######################################################

# Load coverage files from Hoon and the WGS samples. Note that WXS coverage data for samples not processed by 
# Hoon are not available. For these samples, we can default to a read_depth of 100X OR omit.
# These files are not in the database, and will need to be uploaded to the github OR recalculate coverage.
wxs_coverage = "/Users/johnsk/Documents/Life-History/glass-analyses/data/glass_analysis-paired_bam_map_wes.variant_20181017.txt"
wgs_coverage =  "/Users/johnsk/Documents/Life-History/glass-analyses/data/GLASS-WG-sequencing-metrics.csv"

#######################################################
# Necessary packages:
library(tidyverse)
library(DBI)
library(neutralitytestr)

#######################################################
# Establish connection with Floris' database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

###################################################
# Step 1: Prep coverage, purity, and ploidy values
##################################################
# Load additional tables from the database.
pairs = dbReadTable(con,  Id(schema="analysis",table="pairs"))
aliquots = dbReadTable(con,  Id(schema="biospecimen",table="aliquots"))
surgeries = dbReadTable(con,  Id(schema="clinical",table="surgeries"))
titan_param = dbReadTable(con,  Id(schema="analysis",table="titan_params"))

## WXS coverage files.
wxs_coverage_unmerged = read_tsv(wxs_coverage)
wxs_cov_aliquots = wxs_coverage_unmerged %>% 
  select(tm_samplename, tm_bam.mean_wescov) %>% 
  left_join(aliquots, by=c("tm_samplename"="aliquot_id_legacy")) %>% 
  select(aliquot_barcode, avg_coverage = tm_bam.mean_wescov) %>% 
  filter(!grepl("-WGS-", aliquot_barcode))

## Only prepare WGS samples as we don't want merging with WXS.
wgs_coverage_unformatted = read_csv(wgs_coverage)
aliquots_wgs = aliquots %>% 
  filter(grepl("-WGS-", aliquot_barcode))

## Some aliquot_barcodes have changed since the avg_coverage was originally calculated. Remove renamed values.
wgs_cov_aliquots = wgs_coverage_unformatted %>% 
  mutate(sample_barcode = substr(aliquot_id, 1, 15)) %>% 
  left_join(aliquots_wgs, by="sample_barcode") %>% 
  select(aliquot_barcode, avg_coverage = MEDIAN_COVERAGE) %>% 
  filter(!is.na(aliquot_barcode))
  
## Combine all available WGS and WXS average coverage data.
avail_cov = bind_rows(wgs_cov_aliquots, wxs_cov_aliquots)
  
## Combine the avg. coverage, purity, and ploidy values. 
tumor_cov_titan = aliquots %>% 
  left_join(avail_cov, by="aliquot_barcode") %>% 
  distinct() %>% 
  left_join(pairs, by=c("aliquot_barcode"="tumor_barcode")) %>% 
  left_join(titan_param, by="pair_barcode") %>% 
  filter(!grepl("-NB-", aliquot_barcode))

###################################################
# Step 2: Query database to generate vaf data
###################################################
## DATABASE QUERY TO RETRIEVE OPTIMAL SAMPLE PAIRS for THIS analysis.
## We have prioritized WXS for this analysis (to have higher coverage at specific loci), but we may want to consider an WXS ONLY analysis.

###### NO FILTER #########
q_vaf = "WITH
        /*
Define a list of 'selected tumor pairs': select a single primary/recurrent pair per subject/case after removing any aliquots from the blocklist and sorting samples by surgical interval,
so that pairs with the highest surgical interval are prioritized, and pairs with WGS are prioritized over WXS when both are available
*/
selected_tumor_pairs AS
(
  SELECT
  tumor_pair_barcode,
  row_number() OVER (PARTITION BY case_barcode ORDER BY surgical_interval_mo DESC, portion_a ASC, portion_b ASC, substring(tumor_pair_barcode from 27 for 3) DESC) AS priority
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
(CASE WHEN mutect2_call_a AND mutect2_call_b THEN 'shared' WHEN mutect2_call_a AND NOT mutect2_call_b THEN 'primary' WHEN mutect2_call_b AND NOT mutect2_call_a THEN 'recurrent' END) AS status
FROM analysis.master_genotype_comparison gtc
INNER JOIN selected_tumor_pairs stp ON stp.tumor_pair_barcode = gtc.tumor_pair_barcode AND stp.priority = 1
LEFT JOIN analysis.snvs snvs ON snvs.chrom = gtc.chrom AND snvs.pos = gtc.pos AND snvs.alt = gtc.alt
WHERE (mutect2_call_a OR mutect2_call_b) AND (gtc.alt_count_a + gtc.ref_count_a) >= 15 AND (gtc.alt_count_b + gtc.ref_count_b) >= 15"


####### FILTER ########
q_vaf_filtered = "WITH
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
(CASE WHEN mutect2_call_a AND mutect2_call_b THEN 'shared' WHEN mutect2_call_a AND NOT mutect2_call_b THEN 'primary' WHEN mutect2_call_b AND NOT mutect2_call_a THEN 'recurrent' END) AS status
FROM analysis.master_genotype_comparison gtc
INNER JOIN selected_tumor_pairs stp ON stp.tumor_pair_barcode = gtc.tumor_pair_barcode AND stp.priority = 1
LEFT JOIN analysis.snvs snvs ON snvs.chrom = gtc.chrom AND snvs.pos = gtc.pos AND snvs.alt = gtc.alt
WHERE (mutect2_call_a OR mutect2_call_b) AND (gtc.alt_count_a + gtc.ref_count_a) >= 30 AND (gtc.alt_count_b + gtc.ref_count_b) >= 30 AND corrected_call_a = 'HET' AND corrected_call_b = 'HET'"


# Retrieve the relevant variant_allele_frequencies for the optimally selected GLASS pair.
glass_vaf <- dbGetQuery(con, q_vaf)

# If it is desired to have data filtered by coverage and ploidy.
glass_vaf <- dbGetQuery(con, q_vaf_filtered)


###################################################
# Step 3: Broad neutral vaf analyses. Tumor type and subtype.
###################################################
## Package example.
out <- neutralitytest(VAFselection, read_depth = 100.0, cellularity = 0.6, rho = 0.0, ploidy = 2)
plot_all(out)

## Tests across all samples. Appropriate?? Turns out, no. Not UNLESS mutation rates are assumed to be the same.
# PRIMARY
glass_vaf_primary = glass_vaf %>% 
  filter(status %in%c("primary")) 
out_primary = neutralitytest(glass_vaf_primary$vaf_a, fmin = 0.1, fmax = 0.25)
out_primary # return statistics
plot_all(out_primary)

# RECURRENCE
glass_vaf_recur = glass_vaf %>% 
  filter(status %in%c("recurrent")) 
out_recur = neutralitytest(glass_vaf_recur$vaf_b, fmin = 0.1, fmax = 0.25)
out_recur # return statistics
plot_all(out_recur)

# SHARED
glass_vaf_shared_a = glass_vaf %>% 
  filter(status %in%c("shared")) 
out_shared_a = neutralitytest(glass_vaf_shared_a$vaf_a, fmin = 0.1, fmax = 0.25)
out_shared_a # return statistics
plot_all(out_shared_a)
out_shared_b = neutralitytest(glass_vaf_shared_a$vaf_b, fmin = 0.1, fmax = 0.25)
out_shared_b
plot_all(out_shared_b)


###################################################
# Step 4: Per sample neutrality. See if "status" has an impact.
###################################################
### For each tumor calculate the neutraliytest().
## Inclusion of "read_depth" and "cellularity" impact the model.

# For loop method; revisit to clean up the code.
##### "PRIMARY ONLY" in "tumor_a".
glass_vaf_primary = glass_vaf %>% 
  filter(status %in% c("primary")) 

# Store results from all "PRIMARY ONLY" results:
neutral_results_primary = matrix(NA, nrow = length(unique(glass_vaf_primary$case_barcode)), ncol = 13)
colnames(neutral_results_primary) =  c("case_barcode", "aliquot_barcode", "read_depth", "cellularity", "mutation_rate", "model_rsq","model_pval", "area_value", "area_pval","meanDist_value", "meanDist_pval", "Dk_value", "Dk_pval")
neutral_results_primary = as.data.frame(neutral_results_primary)

# Set directory for output.
setwd("/Users/johnsk/Documents/Life-History/glass-analyses/figures/neutral-evolution/primary/")

## Generate results for "PRIMARY ONLY" in "vaf_a" primary sample.
for ( i in 1:length(unique(glass_vaf_primary$case_barcode))) { 
  # Create subsetted data for this case.
  Yi = glass_vaf_primary[glass_vaf_primary$case_barcode==unique(glass_vaf_primary$case_barcode)[i], ]
  xi = unique(Yi$tumor_barcode_a)
  
    # Store sample information.
  neutral_results_primary[i, "case_barcode"] = unique(glass_vaf_primary$case_barcode)[i]
  neutral_results_primary[i, "aliquot_barcode"] = xi
  neutral_results_primary[i, "read_depth"] =  tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "avg_coverage"]
  neutral_results_primary[i, "cellularity"] = tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "purity"]
  
  # The model can sometimes fail. Skip over samples that error out.
  fit_a =  try(neutralitytest(Yi$vaf_a, read_depth = tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "avg_coverage"], cellularity = tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "purity"]), silent = T) 
  if(!inherits(fit_a, "try-error")) {
  neutral_results_primary[i, "mutation_rate"] = fit_a$mutation.rate
  neutral_results_primary[i, "model_rsq"] = fit_a$model$rsq
  neutral_results_primary[i, "model_pval"] = fit_a$model$pval
  neutral_results_primary[i, "area_value"] = fit_a$area$metric
  neutral_results_primary[i, "area_pval"] = fit_a$area$pval
  neutral_results_primary[i, "meanDist_value"] = fit_a$meanDist$metric
  neutral_results_primary[i, "meanDist_pval"] = fit_a$meanDist$pval
  neutral_results_primary[i, "Dk_value"] = fit_a$Dk$metric
  neutral_results_primary[i, "Dk_pval"] = fit_a$Dk$pval
  
  # Print out model plots to directory of interest.
  pdf(paste(unique(glass_vaf_primary$case_barcode)[i],"-primary",".pdf", sep=""), width=8, height=6)
  print(plot_all(fit_a))
  dev.off()
  }
}

# Determine the number of samples that reject the null hypothesis of neutral evolution.
sum(is.na(neutral_results_primary$model_rsq))
sum(neutral_results_primary$model_rsq>0.98, na.rm = T)

######## PRIMARY | SHARED ###########
glass_vaf_pri_shared = glass_vaf %>% 
  filter(status %in% c("primary", "shared")) 

## Store results for Primary | Shared
neutral_results_pri_shared = matrix(NA, nrow = length(unique(glass_vaf_pri_shared$case_barcode)), ncol = 13)
colnames(neutral_results_pri_shared) =  c("case_barcode", "aliquot_barcode", "read_depth", "cellularity", "mutation_rate", "model_rsq","model_pval", "area_value", "area_pval","meanDist_value", "meanDist_pval", "Dk_value", "Dk_pval")
neutral_results_pri_shared = as.data.frame(neutral_results_pri_shared)

# Set directory for output.
setwd("/Users/johnsk/Documents/Life-History/glass-analyses/figures/neutral-evolution/primary-shared/")

## Generate results for "PRIMARY AND SHARED" in "vaf_a" primary sample.
for ( i in 1:length(unique(glass_vaf_pri_shared$case_barcode))) { 
  # Create subsetted data for this case.
  Yi = glass_vaf_pri_shared[glass_vaf_pri_shared$case_barcode==unique(glass_vaf_pri_shared$case_barcode)[i], ]
  xi = unique(Yi$tumor_barcode_a)
  
  # Store sample information.
  neutral_results_pri_shared[i, "case_barcode"] = unique(glass_vaf_pri_shared$case_barcode)[i]
  neutral_results_pri_shared[i, "aliquot_barcode"] = xi
  neutral_results_pri_shared[i, "read_depth"] =  tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "avg_coverage"]
  neutral_results_pri_shared[i, "cellularity"] = tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "purity"]
  
  # The model can sometimes fail. Skip over samples that error out.
  fit_a =  try(neutralitytest(Yi$vaf_a, read_depth = tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "avg_coverage"], cellularity = tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "purity"]), silent = T) 
  if(!inherits(fit_a, "try-error")) {
    neutral_results_pri_shared[i, "mutation_rate"] = fit_a$mutation.rate
    neutral_results_pri_shared[i, "model_rsq"] = fit_a$model$rsq
    neutral_results_pri_shared[i, "model_pval"] = fit_a$model$pval
    neutral_results_pri_shared[i, "area_value"] = fit_a$area$metric
    neutral_results_pri_shared[i, "area_pval"] = fit_a$area$pval
    neutral_results_pri_shared[i, "meanDist_value"] = fit_a$meanDist$metric
    neutral_results_pri_shared[i, "meanDist_pval"] = fit_a$meanDist$pval
    neutral_results_pri_shared[i, "Dk_value"] = fit_a$Dk$metric
    neutral_results_pri_shared[i, "Dk_pval"] = fit_a$Dk$pval
    
    # Print out model plots to directory of interest.
    pdf(paste(unique(glass_vaf_pri_shared$case_barcode)[i],"-shared_primary",".pdf", sep=""), width=8, height=6)
    print(plot_all(fit_a))
    dev.off()
  }
}
sum(is.na(neutral_results_pri_shared$model_rsq))
sum(neutral_results_pri_shared$model_rsq>0.98, na.rm = T)


#### SHARED | RECURRENT  #####
# Investigating tumor_b.
glass_vaf_recur_shared = glass_vaf %>% 
  filter(status %in% c("recurrent", "shared")) 

## Store results for Primary| Shared
neutral_results_recur_shared = matrix(NA, nrow = length(unique(glass_vaf_recur_shared$case_barcode)), ncol = 13)
colnames(neutral_results_recur_shared) =  c("case_barcode", "aliquot_barcode", "read_depth", "cellularity", "mutation_rate", "model_rsq","model_pval", "area_value", "area_pval","meanDist_value", "meanDist_pval", "Dk_value", "Dk_pval")
neutral_results_recur_shared = as.data.frame(neutral_results_recur_shared)

# Set directory for output.
setwd("/Users/johnsk/Documents/Life-History/glass-analyses/figures/neutral-evolution-unfiltered/recurrent_shared/")

## Generate results for "RECURRENT AND SHARED" in "vaf_b" recurrent sample.
for ( i in 1:length(unique(glass_vaf_recur_shared$case_barcode))) { 
  # Create subsetted data for this case.
  Yi = glass_vaf_recur_shared[glass_vaf_recur_shared$case_barcode==unique(glass_vaf_recur_shared$case_barcode)[i], ]
  xi = unique(Yi$tumor_barcode_a)
  
  # Store sample information.
  neutral_results_recur_shared[i, "case_barcode"] = unique(glass_vaf_recur_shared$case_barcode)[i]
  neutral_results_recur_shared[i, "aliquot_barcode"] = xi
  neutral_results_recur_shared[i, "read_depth"] =  tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "avg_coverage"]
  neutral_results_recur_shared[i, "cellularity"] = tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "purity"]
  
  # The model can sometimes fail. Skip over samples that error out.
  fit_a =  try(neutralitytest(Yi$vaf_b, read_depth = tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "avg_coverage"], cellularity = tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "purity"]), silent = T) 
  if(!inherits(fit_a, "try-error")) {
    neutral_results_recur_shared[i, "mutation_rate"] = fit_a$mutation.rate
    neutral_results_recur_shared[i, "model_rsq"] = fit_a$model$rsq
    neutral_results_recur_shared[i, "model_pval"] = fit_a$model$pval
    neutral_results_recur_shared[i, "area_value"] = fit_a$area$metric
    neutral_results_recur_shared[i, "area_pval"] = fit_a$area$pval
    neutral_results_recur_shared[i, "meanDist_value"] = fit_a$meanDist$metric
    neutral_results_recur_shared[i, "meanDist_pval"] = fit_a$meanDist$pval
    neutral_results_recur_shared[i, "Dk_value"] = fit_a$Dk$metric
    neutral_results_recur_shared[i, "Dk_pval"] = fit_a$Dk$pval
    
    # Print out model plots to directory of interest.
    pdf(paste(unique(glass_vaf_recur_shared$case_barcode)[i],"-shared_recur",".pdf", sep=""), width=8, height=6)
    print(plot_all(fit_a))
    dev.off()
  }
}
sum(is.na(neutral_results_recur_shared$model_rsq))
sum(neutral_results_recur_shared$model_rsq>0.98, na.rm = T)


####  RECURRENT  #####
glass_vaf_recur = glass_vaf %>% 
  filter(status %in% c("recurrent")) 

## Store results for RECURERNT
neutral_results_recur = matrix(NA, nrow = length(unique(glass_vaf_recur$case_barcode)), ncol = 13)
colnames(neutral_results_recur) =  c("case_barcode", "aliquot_barcode", "read_depth", "cellularity", "mutation_rate", "model_rsq","model_pval", "area_value", "area_pval","meanDist_value", "meanDist_pval", "Dk_value", "Dk_pval")
neutral_results_recur = as.data.frame(neutral_results_recur)

# Set directory for output.
setwd("/Users/johnsk/Documents/Life-History/glass-analyses/figures/neutral-evolution/recurrent/")

## Generate results for "RECURRENT" in "vaf_b" recurrent sample.
for ( i in 1:length(unique(glass_vaf_recur$case_barcode))) { 
  # Create subsetted data for this case.
  Yi = glass_vaf_recur[glass_vaf_recur$case_barcode==unique(glass_vaf_recur$case_barcode)[i], ]
  xi = unique(Yi$tumor_barcode_a)
  
  # Store sample information.
  neutral_results_recur[i, "case_barcode"] = unique(glass_vaf_recur$case_barcode)[i]
  neutral_results_recur[i, "aliquot_barcode"] = xi
  neutral_results_recur[i, "read_depth"] =  tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "avg_coverage"]
  neutral_results_recur[i, "cellularity"] = tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "purity"]
  
  # The model can sometimes fail. Skip over samples that error out.
  fit_a =  try(neutralitytest(Yi$vaf_b, read_depth = tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "avg_coverage"], cellularity = tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "purity"]), silent = T) 
  if(!inherits(fit_a, "try-error")) {
    neutral_results_recur[i, "mutation_rate"] = fit_a$mutation.rate
    neutral_results_recur[i, "model_rsq"] = fit_a$model$rsq
    neutral_results_recur[i, "model_pval"] = fit_a$model$pval
    neutral_results_recur[i, "area_value"] = fit_a$area$metric
    neutral_results_recur[i, "area_pval"] = fit_a$area$pval
    neutral_results_recur[i, "meanDist_value"] = fit_a$meanDist$metric
    neutral_results_recur[i, "meanDist_pval"] = fit_a$meanDist$pval
    neutral_results_recur[i, "Dk_value"] = fit_a$Dk$metric
    neutral_results_recur[i, "Dk_pval"] = fit_a$Dk$pval
    
    # Print out model plots to directory of interest.
    pdf(paste(unique(glass_vaf_recur$case_barcode)[i],"-recur",".pdf", sep=""), width=8, height=6)
    print(plot_all(fit_a))
    dev.off()
  }
}


###### ALL ###########
glass_vaf_all_a = glass_vaf %>% 
  filter(status %in% c("primary", "recurrent","shared")) 

## Store results for ALL
neutral_results_all_a = matrix(NA, nrow = length(unique(glass_vaf_all_a$case_barcode)), ncol = 13)
colnames(neutral_results_all_a) =  c("case_barcode", "aliquot_barcode", "read_depth", "cellularity", "mutation_rate", "model_rsq","model_pval", "area_value", "area_pval","meanDist_value", "meanDist_pval", "Dk_value", "Dk_pval")
neutral_results_all_a = as.data.frame(neutral_results_all_a)

# Set directory for output.
setwd("/Users/johnsk/Documents/Life-History/glass-analyses/figures/neutral-evolution/all_a/")

## Generate results for "ALL" in "vaf_a" (primary) and "vaf_b" (recurrent) sample.
for ( i in 1:length(unique(glass_vaf_all_a$case_barcode))) { 
  # Create subsetted data for this case.
  Yi = glass_vaf_all_a[glass_vaf_all_a$case_barcode==unique(glass_vaf_all_a$case_barcode)[i], ]
  xi = unique(Yi$tumor_barcode_a)
  
  # Store sample information.
  neutral_results_all_a[i, "case_barcode"] = unique(glass_vaf_all_a$case_barcode)[i]
  neutral_results_all_a[i, "aliquot_barcode"] = xi
  neutral_results_all_a[i, "read_depth"] =  tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "avg_coverage"]
  neutral_results_all_a[i, "cellularity"] = tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "purity"]
  
  # The model can sometimes fail. Skip over samples that error out.
  fit_a =  try(neutralitytest(Yi$vaf_a, read_depth = tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "avg_coverage"], cellularity = tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "purity"]), silent = T) 
  if(!inherits(fit_a, "try-error")) {
    neutral_results_all_a[i, "mutation_rate"] = fit_a$mutation.rate
    neutral_results_all_a[i, "model_rsq"] = fit_a$model$rsq
    neutral_results_all_a[i, "model_pval"] = fit_a$model$pval
    neutral_results_all_a[i, "area_value"] = fit_a$area$metric
    neutral_results_all_a[i, "area_pval"] = fit_a$area$pval
    neutral_results_all_a[i, "meanDist_value"] = fit_a$meanDist$metric
    neutral_results_all_a[i, "meanDist_pval"] = fit_a$meanDist$pval
    neutral_results_all_a[i, "Dk_value"] = fit_a$Dk$metric
    neutral_results_all_a[i, "Dk_pval"] = fit_a$Dk$pval
    
    # Print out model plots to directory of interest.
    pdf(paste(unique(glass_vaf_all_a$case_barcode)[i],"-all_a",".pdf", sep=""), width=8, height=6)
    print(plot_all(fit_a))
    dev.off()
  }
}

#### tumor_b ######
glass_vaf_all_b = glass_vaf %>% 
  filter(status %in% c("primary", "recurrent","shared")) 

## Store results for ALL
neutral_results_all_b = matrix(NA, nrow = length(unique(glass_vaf_all_b$case_barcode)), ncol = 13)
colnames(neutral_results_all_b) =  c("case_barcode", "aliquot_barcode", "read_depth", "cellularity", "mutation_rate", "model_rsq","model_pval", "area_value", "area_pval","meanDist_value", "meanDist_pval", "Dk_value", "Dk_pval")
neutral_results_all_b = as.data.frame(neutral_results_all_b)

# Set directory for output.
setwd("/Users/johnsk/Documents/Life-History/glass-analyses/figures/neutral-evolution/all_b/")

## Generate results for "ALL" in "vaf_a" (primary) and "vaf_b" (recurrent) sample.
for ( i in 1:length(unique(glass_vaf_all_b$case_barcode))) { 
  # Create subsetted data for this case.
  Yi = glass_vaf_all_b[glass_vaf_all_b$case_barcode==unique(glass_vaf_all_b$case_barcode)[i], ]
  xi = unique(Yi$tumor_barcode_a)
  
  # Store sample information.
  neutral_results_all_b[i, "case_barcode"] = unique(glass_vaf_all_b$case_barcode)[i]
  neutral_results_all_b[i, "aliquot_barcode"] = xi
  neutral_results_all_b[i, "read_depth"] =  tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "avg_coverage"]
  neutral_results_all_b[i, "cellularity"] = tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "purity"]
  
  # The model can sometimes fail. Skip over samples that error out.
  fit_a =  try(neutralitytest(Yi$vaf_b, read_depth = tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "avg_coverage"], cellularity = tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "purity"]), silent = T) 
  if(!inherits(fit_a, "try-error")) {
    neutral_results_all_b[i, "mutation_rate"] = fit_a$mutation.rate
    neutral_results_all_b[i, "model_rsq"] = fit_a$model$rsq
    neutral_results_all_b[i, "model_pval"] = fit_a$model$pval
    neutral_results_all_b[i, "area_value"] = fit_a$area$metric
    neutral_results_all_b[i, "area_pval"] = fit_a$area$pval
    neutral_results_all_b[i, "meanDist_value"] = fit_a$meanDist$metric
    neutral_results_all_b[i, "meanDist_pval"] = fit_a$meanDist$pval
    neutral_results_all_b[i, "Dk_value"] = fit_a$Dk$metric
    neutral_results_all_b[i, "Dk_pval"] = fit_a$Dk$pval
    
    # Print out model plots to directory of interest.
    pdf(paste(unique(glass_vaf_all_b$case_barcode)[i],"-all_b",".pdf", sep=""), width=8, height=6)
    print(plot_all(fit_a))
    dev.off()
  }
}


#### SHARED  #####
#### shared_a
glass_vaf_shared_a = glass_vaf %>% 
  filter(status %in% c("shared")) 

## Store results for SHARED
neutral_results_shared_a = matrix(NA, nrow = length(unique(glass_vaf_shared_a$case_barcode)), ncol = 13)
colnames(neutral_results_shared_a) =  c("case_barcode", "aliquot_barcode", "read_depth", "cellularity", "mutation_rate", "model_rsq","model_pval", "area_value", "area_pval","meanDist_value", "meanDist_pval", "Dk_value", "Dk_pval")
neutral_results_shared_a = as.data.frame(neutral_results_shared_a)

# Set directory for output.
setwd("/Users/johnsk/Documents/Life-History/glass-analyses/figures/neutral-evolution/shared_a/")

## Store results for "SHARED" in "vaf_a" primary sample.
for ( i in 1:length(unique(glass_vaf_shared_a$case_barcode))) { 
  # Create subsetted data for this case.
  Yi = glass_vaf_shared_a[glass_vaf_shared_a$case_barcode==unique(glass_vaf_shared_a$case_barcode)[i], ]
  xi = unique(Yi$tumor_barcode_a)
  
  # Store sample information.
  neutral_results_shared_a[i, "case_barcode"] = unique(glass_vaf_shared_a$case_barcode)[i]
  neutral_results_shared_a[i, "aliquot_barcode"] = xi
  neutral_results_shared_a[i, "read_depth"] =  tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "avg_coverage"]
  neutral_results_shared_a[i, "cellularity"] = tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "purity"]
  
  # The model can sometimes fail. Skip over samples that error out.
  fit_a =  try(neutralitytest(Yi$vaf_a, read_depth = tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "avg_coverage"], cellularity = tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "purity"]), silent = T) 
  if(!inherits(fit_a, "try-error")) {
    neutral_results_shared_a[i, "mutation_rate"] = fit_a$mutation.rate
    neutral_results_shared_a[i, "model_rsq"] = fit_a$model$rsq
    neutral_results_shared_a[i, "model_pval"] = fit_a$model$pval
    neutral_results_shared_a[i, "area_value"] = fit_a$area$metric
    neutral_results_shared_a[i, "area_pval"] = fit_a$area$pval
    neutral_results_shared_a[i, "meanDist_value"] = fit_a$meanDist$metric
    neutral_results_shared_a[i, "meanDist_pval"] = fit_a$meanDist$pval
    neutral_results_shared_a[i, "Dk_value"] = fit_a$Dk$metric
    neutral_results_shared_a[i, "Dk_pval"] = fit_a$Dk$pval
    
    # Print out model plots to directory of interest.
    pdf(paste(unique(glass_vaf_shared_a$case_barcode)[i],"-shared_a",".pdf", sep=""), width=8, height=6)
    print(plot_all(fit_a))
    dev.off()
  }
}


#### shared_b  ######
glass_vaf_shared_b = glass_vaf %>% 
  filter(status %in% c("shared")) 

## Store results for SHARED
neutral_results_shared_b = matrix(NA, nrow = length(unique(glass_vaf_shared_b$case_barcode)), ncol = 13)
colnames(neutral_results_shared_b) =  c("case_barcode", "aliquot_barcode", "read_depth", "cellularity", "mutation_rate", "model_rsq","model_pval", "area_value", "area_pval","meanDist_value", "meanDist_pval", "Dk_value", "Dk_pval")
neutral_results_shared_b = as.data.frame(neutral_results_shared_b)

# Set directory for output.
setwd("/Users/johnsk/Documents/Life-History/glass-analyses/figures/neutral-evolution/shared_b/")

## Store results for "SHARED" in "vaf_b" recurrent sample.
for ( i in 1:length(unique(glass_vaf_shared_b$case_barcode))) { 
  # Create subsetted data for this case.
  Yi = glass_vaf_shared_b[glass_vaf_shared_b$case_barcode==unique(glass_vaf_shared_b$case_barcode)[i], ]
  xi = unique(Yi$tumor_barcode_a)
  
  # Store sample information.
  neutral_results_shared_b[i, "case_barcode"] = unique(glass_vaf_shared_b$case_barcode)[i]
  neutral_results_shared_b[i, "aliquot_barcode"] = xi
  neutral_results_shared_b[i, "read_depth"] =  tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "avg_coverage"]
  neutral_results_shared_b[i, "cellularity"] = tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "purity"]
  
  # The model can sometimes fail. Skip over samples that error out.
  fit_a =  try(neutralitytest(Yi$vaf_b, read_depth = tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "avg_coverage"], cellularity = tumor_cov_titan[tumor_cov_titan$aliquot_barcode==xi, "purity"]), silent = T) 
  if(!inherits(fit_a, "try-error")) {
    neutral_results_shared_b[i, "mutation_rate"] = fit_a$mutation.rate
    neutral_results_shared_b[i, "model_rsq"] = fit_a$model$rsq
    neutral_results_shared_b[i, "model_pval"] = fit_a$model$pval
    neutral_results_shared_b[i, "area_value"] = fit_a$area$metric
    neutral_results_shared_b[i, "area_pval"] = fit_a$area$pval
    neutral_results_shared_b[i, "meanDist_value"] = fit_a$meanDist$metric
    neutral_results_shared_b[i, "meanDist_pval"] = fit_a$meanDist$pval
    neutral_results_shared_b[i, "Dk_value"] = fit_a$Dk$metric
    neutral_results_shared_b[i, "Dk_pval"] = fit_a$Dk$pval
    
    # Print out model plots to directory of interest.
    pdf(paste(unique(glass_vaf_shared_b$case_barcode)[i],"-shared_b",".pdf", sep=""), width=8, height=6)
    print(plot_all(fit_a))
    dev.off()
  }
}

######## Summary Tables #########
write.table(neutral_results_all_a, file = "/Users/johnsk/Documents/Life-History/glass-analyses/figures/neutral-evolution/tabular_results/neutral_results_all_a.txt", sep="\t", row.names = F, col.names = T, quote = F)
write.table(neutral_results_all_b, file = "/Users/johnsk/Documents/Life-History/glass-analyses/figures/neutral-evolution/tabular_results/neutral_results_all_b.txt", sep="\t", row.names = F, col.names = T, quote = F)
write.table(neutral_results_primary, file = "/Users/johnsk/Documents/Life-History/glass-analyses/figures/neutral-evolution/tabular_results/neutral_results_primary.txt", sep="\t", row.names = F, col.names = T, quote = F)
write.table(neutral_results_pri_shared, file = "/Users/johnsk/Documents/Life-History/glass-analyses/figures/neutral-evolution/tabular_results/neutral_results_pri_shared.txt", sep="\t", row.names = F, col.names = T, quote = F)
write.table(neutral_results_recur, file = "/Users/johnsk/Documents/Life-History/glass-analyses/figures/neutral-evolution/tabular_results/neutral_results_recur.txt", sep="\t", row.names = F, col.names = T, quote = F)
write.table(neutral_results_recur_shared, file = "/Users/johnsk/Documents/Life-History/glass-analyses/figures/neutral-evolution-unfiltered/tabular_results/neutral_results_shared_recur.txt", sep="\t", row.names = F, col.names = T, quote = F)
write.table(neutral_results_shared_a, file = "/Users/johnsk/Documents/Life-History/glass-analyses/figures/neutral-evolution/tabular_results/neutral_results_shared_a.txt", sep="\t", row.names = F, col.names = T, quote = F)
write.table(neutral_results_shared_b, file = "/Users/johnsk/Documents/Life-History/glass-analyses/figures/neutral-evolution/tabular_results/neutral_results_shared_b.txt", sep="\t", row.names = F, col.names = T, quote = F)



