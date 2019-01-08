#######################################################
# NeutralitytestR applied to coverage filtered and ploidy filtered variants (vaf thresholds 0.05)
# Updated: 2019.01.07
# Author: Kevin J.
#######################################################

# Working directory for the GLASS-analysis project. 
mybasedir = "/Volumes/verhaak-lab/GLASS-analysis/results/neutral_evolution/neutral-evolution-mutect2"
setwd(mybasedir)

#######################################################
# Necessary packages:
library(tidyverse)
library(DBI)
library(neutralitytestr)
library(ggExtra)

#######################################################
# Establish connection with Floris' database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

# Load additional tables from the database.
pairs = dbReadTable(con,  Id(schema="analysis",table="pairs"))
tumor_pairs = dbReadTable(con,  Id(schema="analysis", table="tumor_pairs"))
aliquots = dbReadTable(con,  Id(schema="biospecimen", table="aliquots"))
surgeries = dbReadTable(con,  Id(schema="clinical", table="surgeries"))
drivers = dbReadTable(con,  Id(schema="analysis", table="driver_status"))
titan_param = dbReadTable(con,  Id(schema="analysis",table="titan_params"))
mut_freq = dbReadTable(con,  Id(schema="analysis",table="mutation_freq"))

# Combine with pairs file to access "tumor_barcode".
titan_info = titan_param %>% 
  left_join(pairs, by="pair_barcode")


####### FILTER ########
# Floris updated the neutralitytestR input query to retrieve all possible tumor pairs AND designate mutation status.
neutrality_input_mutect2  = read_file("/Users/johnsk/Documents/Life-History/glass-analyses/scripts/neutrality-testr-input-mutect2.sql")

# Currently, there is a coverage filter requirement (30X in both tumor_a and tumor_b) PLUS a vaf threshold of 0.05.
glass_vaf_raw <- dbGetQuery(con, neutrality_input_mutect2)

# NOTE that there will be mutations with status `NA` for mutations because there was a Mutect2 call in one of the samples, BUT
# the mutation DID NOT reach a VAF > 0.05.

# Add some additional descriptors and filters.
glass_vaf = glass_vaf_raw %>% 
  mutate(tumor_a_type = substr(tumor_barcode_a, 14, 15),
         tumor_b_type = substr(tumor_barcode_b, 14, 15), 
         tumor_pair_barcode = paste(substr(tumor_barcode_a, 1, 18), substr(tumor_barcode_b, 14, 19), substr(tumor_barcode_a, 21, 23), sep="-"),
         comparison = paste(tumor_a_type, tumor_b_type, sep="-"),
         comparison_type = ifelse(tumor_a_type==tumor_b_type, "multisector", "longitudinal")) %>% 
  filter(!comparison%in%c("R1-M1", "R2-M1", "TP-M1")) %>% 
  filter(comparison_type == "longitudinal")


# First, output vaf_a x vaf_b plots using this 30X coverage in each tumor criteria.
setwd("/Volumes/verhaak-lab/GLASS-analysis/results/neutral_evolution/neutral-evolution-mutect2/vaf_compare")

## Generate pdf of "vaf_a" vs. "vaf_b" for all longitudinal samples.
## Not that the shared/private mutations will be highly dependent on selected threshold.
for ( i in 1:length(unique(glass_vaf$tumor_pair_barcode))) { 
  Yi = glass_vaf[glass_vaf$tumor_pair_barcode==unique(glass_vaf$tumor_pair_barcode)[i], ]
  sum_dat = Yi %>% group_by(status) %>% summarise(mutation_count = n())
  p1 = ggplot(Yi, aes(vaf_a, vaf_b, color = status)) + geom_point(alpha=0.4) +  coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +ggtitle(sprintf("%s", unique(glass_vaf$tumor_pair_barcode)[i])) +
    scale_color_discrete(labels= paste(sum_dat$status, " (n=", sum_dat$mutation_count, ")", sep = "")) + theme_bw()
  p2 = ggMarginal(p1, type = "histogram", groupFill = TRUE)
  pdf(paste(unique(glass_vaf$tumor_pair_barcode)[i],".pdf", sep=""), width=9, height=6)
  print(p2)
  dev.off()
}


###################################################
##### Data includes all tumor comparisons     #####
###################################################
## Inclusion of "ploidy" and "cellularity" variables.
#For each tumor calculate the neutraliytest().

##### Private in "tumor_a". Originally described as private, but P works better.
glass_vaf_primary = glass_vaf %>% 
  filter(status %in% c("P")) 

# Determine the number of mutations that are classified as primary per patient.
primary_mutation_counts = glass_vaf_primary %>% 
  filter(vaf_a > 0.1 & vaf_a < 0.25) %>% 
  group_by(case_barcode) %>% 
  summarize(mutation_count = n())

# Store results from all "PRIMARY ONLY" results:
neutral_results_primary = matrix(NA, nrow = length(unique(glass_vaf_primary$tumor_pair_barcode)), ncol = 14)
colnames(neutral_results_primary) =  c("tumor_pair_barcode", "tumor_barcode", "comparison", "ploidy", "cellularity", "mutation_rate", "model_rsq","model_pval", "area_value", "area_pval","meanDist_value", "meanDist_pval", "Dk_value", "Dk_pval")
neutral_results_primary = as.data.frame(neutral_results_primary)

# Set directory for output.
setwd("/Volumes/verhaak-lab/GLASS-analysis/results/neutral_evolution/neutral-evolution-mutect2/tumor_a_private/")

## Generate results for "tumor_a" using "vaf_a".
for ( i in 1:length(unique(glass_vaf_primary$tumor_pair_barcode))){ 
  # Create subsetted data for this case.
  Yi = glass_vaf_primary[glass_vaf_primary$tumor_pair_barcode==unique(glass_vaf_primary$tumor_pair_barcode)[i], ]
  xi = unique(Yi$tumor_barcode_a)
  zi = unique(Yi$tumor_pair_barcode)
  print(sprintf("Analyzing private mutations for: %s in tumor pair: %s",  xi, zi))
  # Store sample information.
  neutral_results_primary[i, "tumor_pair_barcode"] = unique(glass_vaf_primary$tumor_pair_barcode)[i]
  neutral_results_primary[i, "tumor_barcode"] = xi
  neutral_results_primary[i, "comparison"] = Yi[1, "comparison"]
  neutral_results_primary[i, "ploidy"] = titan_info[titan_info$tumor_barcode==xi, "ploidy"]
  neutral_results_primary[i, "cellularity"] = titan_info[titan_info$tumor_barcode==xi, "purity"]
  
  # The model can sometimes fail. Skip over samples that error out.
  fit_a =  try(neutralitytest(Yi$vaf_a, ploidy = titan_info[titan_info$tumor_barcode==xi, "ploidy"], cellularity = titan_info[titan_info$tumor_barcode==xi, "purity"]), silent = T) 
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
    pdf(paste(unique(glass_vaf_primary$tumor_pair_barcode)[i],"-private_a",".pdf", sep=""), width=9, height=6)
    print(plot_all(fit_a))
    dev.off()
  }
}


#### SHARED  #####
#### shared_a
glass_vaf_shared_a = glass_vaf %>% 
  filter(status %in% c("S")) 

# Determine the number of mutations that are classified as shared_a per patient.
shared_a_mutation_counts = glass_vaf_shared_a %>% 
  filter(vaf_a > 0.1 & vaf_a < 0.25) %>% 
  group_by(case_barcode) %>% 
  summarize(mutation_count = n())

## Store results for SHARED
neutral_results_shared_a = matrix(NA, nrow = length(unique(glass_vaf_shared_a$tumor_pair_barcode)), ncol = 14)
colnames(neutral_results_shared_a) =  c("tumor_pair_barcode", "tumor_barcode", "comparison", "ploidy", "cellularity", "mutation_rate", "model_rsq","model_pval", "area_value", "area_pval","meanDist_value", "meanDist_pval", "Dk_value", "Dk_pval")
neutral_results_shared_a = as.data.frame(neutral_results_shared_a)

# Set directory for output.
setwd("/Volumes/verhaak-lab/GLASS-analysis/results/neutral_evolution/neutral-evolution-mutect2/shared_a/")

## Generate results for "tumor_a" using "vaf_a".
for ( i in 1:length(unique(glass_vaf_shared_a$tumor_pair_barcode))){ 
  # Create subsetted data for this case.
  Yi = glass_vaf_shared_a[glass_vaf_shared_a$tumor_pair_barcode==unique(glass_vaf_shared_a$tumor_pair_barcode)[i], ]
  xi = unique(Yi$tumor_barcode_a)
  zi = unique(Yi$tumor_pair_barcode)
  print(sprintf("Analyzing private mutations for: %s in tumor pair: %s",  xi, zi))
  # Store sample information.
  neutral_results_shared_a[i, "tumor_pair_barcode"] = unique(glass_vaf_shared_a$tumor_pair_barcode)[i]
  neutral_results_shared_a[i, "tumor_barcode"] = xi
  neutral_results_shared_a[i, "comparison"] = Yi[1, "comparison"]
  neutral_results_shared_a[i, "ploidy"] = titan_info[titan_info$tumor_barcode==xi, "ploidy"]
  neutral_results_shared_a[i, "cellularity"] = titan_info[titan_info$tumor_barcode==xi, "purity"]
  
  # The model can sometimes fail. Skip over samples that error out.
  fit_a =  try(neutralitytest(Yi$vaf_a, ploidy = titan_info[titan_info$tumor_barcode==xi, "ploidy"], cellularity = titan_info[titan_info$tumor_barcode==xi, "purity"]), silent = T) 
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
    pdf(paste(unique(glass_vaf_shared_a$tumor_pair_barcode)[i],"-shared_a",".pdf", sep=""), width=9, height=6)
    print(plot_all(fit_a))
    dev.off()
  }
}



#### shared_b #####
glass_vaf_shared_b = glass_vaf %>% 
  filter(status %in% c("S")) 

# Determine the number of mutations that are classified as shared_b per patient.
shared_b_mutation_counts = glass_vaf_shared_b %>% 
  filter(vaf_b > 0.1 & vaf_b < 0.25) %>% 
  group_by(case_barcode) %>% 
  summarize(mutation_count = n())

## Store results for SHARED
neutral_results_shared_b = matrix(NA, nrow = length(unique(glass_vaf_shared_b$tumor_pair_barcode)), ncol = 14)
colnames(neutral_results_shared_b) =  c("tumor_pair_barcode", "tumor_barcode", "comparison", "ploidy", "cellularity", "mutation_rate", "model_rsq","model_pval", "area_value", "area_pval","meanDist_value", "meanDist_pval", "Dk_value", "Dk_pval")
neutral_results_shared_b = as.data.frame(neutral_results_shared_b)

# Set directory for output.
setwd("/Volumes/verhaak-lab/GLASS-analysis/results/neutral_evolution/neutral-evolution-mutect2/shared_b/")

## Generate results for "tumor_barcode_b" using "vaf_b".
for ( i in 1:length(unique(glass_vaf_shared_b$tumor_pair_barcode))){ 
  # Create subsetted data for this case.
  Yi = glass_vaf_shared_b[glass_vaf_shared_b$tumor_pair_barcode==unique(glass_vaf_shared_b$tumor_pair_barcode)[i], ]
  xi = unique(Yi$tumor_barcode_b)
  zi = unique(Yi$tumor_pair_barcode)
  print(sprintf("Analyzing private mutations for: %s in tumor pair: %s",  xi, zi))
  # Store sample information.
  neutral_results_shared_b[i, "tumor_pair_barcode"] = zi
  neutral_results_shared_b[i, "tumor_barcode"] = xi
  neutral_results_shared_b[i, "comparison"] = Yi[1, "comparison"]
  neutral_results_shared_b[i, "ploidy"] = titan_info[titan_info$tumor_barcode==xi, "ploidy"]
  neutral_results_shared_b[i, "cellularity"] = titan_info[titan_info$tumor_barcode==xi, "purity"]
  
  # The model can sometimes fail. Skip over samples that error out.
  fit_a =  try(neutralitytest(Yi$vaf_b, ploidy = titan_info[titan_info$tumor_barcode==xi, "ploidy"], cellularity = titan_info[titan_info$tumor_barcode==xi, "purity"]), silent = T) 
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
    pdf(paste(unique(glass_vaf_shared_b$tumor_pair_barcode)[i],"-shared_b",".pdf", sep=""), width=9, height=6)
    print(plot_all(fit_a))
    dev.off()
  }
}


####  RECURRENT  #####
glass_vaf_recur = glass_vaf %>% 
  filter(status %in% c("R")) 

# Determine the number of mutations that are classified as private to tumor_b.
recurrence_mutation_counts = glass_vaf_recur %>% 
  filter(vaf_b > 0.1 & vaf_b < 0.25) %>% 
  group_by(case_barcode) %>% 
  summarize(mutation_count = n())

## Store results for what's called "recurrence" although it really refers to private to tumor_b.
neutral_results_recur = matrix(NA, nrow = length(unique(glass_vaf_recur$tumor_pair_barcode)), ncol = 14)
colnames(neutral_results_recur) =  c("tumor_pair_barcode", "tumor_barcode", "comparison", "ploidy", "cellularity", "mutation_rate", "model_rsq","model_pval", "area_value", "area_pval","meanDist_value", "meanDist_pval", "Dk_value", "Dk_pval")
neutral_results_recur = as.data.frame(neutral_results_recur)

# Set directory for output.
setwd("/Volumes/verhaak-lab/GLASS-analysis/results/neutral_evolution/neutral-evolution-mutect2/tumor_b_private/")

## Generate results for "tumor_barcode_b" using "vaf_b".
for ( i in 1:length(unique(glass_vaf_recur$tumor_pair_barcode))){ 
  # Create subsetted data for this case.
  Yi = glass_vaf_recur[glass_vaf_recur$tumor_pair_barcode==unique(glass_vaf_recur$tumor_pair_barcode)[i], ]
  xi = unique(Yi$tumor_barcode_b)
  zi = unique(Yi$tumor_pair_barcode)
  print(sprintf("Analyzing private mutations for: %s in tumor pair: %s",  xi, zi))
  # Store sample information.
  neutral_results_recur[i, "tumor_pair_barcode"] = zi
  neutral_results_recur[i, "tumor_barcode"] = xi
  neutral_results_recur[i, "comparison"] = Yi[1, "comparison"]
  neutral_results_recur[i, "ploidy"] = titan_info[titan_info$tumor_barcode==xi, "ploidy"]
  neutral_results_recur[i, "cellularity"] = titan_info[titan_info$tumor_barcode==xi, "purity"]
  
  # The model can sometimes fail. Skip over samples that error out.
  fit_a =  try(neutralitytest(Yi$vaf_b, ploidy = titan_info[titan_info$tumor_barcode==xi, "ploidy"], cellularity = titan_info[titan_info$tumor_barcode==xi, "purity"]), silent = T) 
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
    pdf(paste(unique(glass_vaf_recur$tumor_pair_barcode)[i],"-recur",".pdf", sep=""), width=9, height=6)
    print(plot_all(fit_a))
    dev.off()
  }
}


######## Output: Summary Tables #########
write.table(neutral_results_primary, file = "/Volumes/verhaak-lab/GLASS-analysis/results/neutral_evolution/neutral-evolution-mutect2/tabular_results/neutral_results_primary.txt", sep="\t", row.names = F, col.names = T, quote = F)
write.table(neutral_results_recur, file = "/Volumes/verhaak-lab/GLASS-analysis/results/neutral_evolution/neutral-evolution-mutect2/tabular_results/neutral_results_recur.txt", sep="\t", row.names = F, col.names = T, quote = F)
write.table(neutral_results_shared_a, file = "/Volumes/verhaak-lab/GLASS-analysis/results/neutral_evolution/neutral-evolution-mutect2/tabular_results/neutral_results_shared_a.txt", sep="\t", row.names = F, col.names = T, quote = F)
write.table(neutral_results_shared_b, file = "/Volumes/verhaak-lab/GLASS-analysis/results/neutral_evolution/neutral-evolution-mutect2/tabular_results/neutral_results_shared_b.txt", sep="\t", row.names = F, col.names = T, quote = F)


#### Primary, Shared, and Recurrent mutations #####
# Mutations private to the primary tumor.
primary = read_tsv("/Volumes/verhaak-lab/GLASS-analysis/results/neutral_evolution/neutral-evolution-mutect2/tabular_results/neutral_results_primary.txt")
primary = primary %>% 
  mutate(case_barcode = substr(tumor_barcode, 1, 12)) %>% 
  left_join(primary_mutation_counts, by="case_barcode") %>% 
  mutate(analysis = "primary") %>% 
  filter(ploidy < 3, mutation_count >= 12) %>% 
  select(tumor_pair_barcode, primary_rsq = model_rsq, primary_pval = model_pval, analysis, mutation_count)

# shared_a = the shared mutations found in the primary tumor.
shared_a = read_tsv("/Volumes/verhaak-lab/GLASS-analysis/results/neutral_evolution/neutral-evolution-mutect2/tabular_results/neutral_results_shared_a.txt")
shared_a = shared_a %>% 
  mutate(case_barcode = substr(tumor_barcode, 1, 12)) %>% 
  left_join(shared_a_mutation_counts, by="case_barcode") %>% 
  mutate(analysis = "shared_a") %>% 
  filter(ploidy < 3, mutation_count >= 12) %>% 
  select(tumor_pair_barcode, shared_a_rsq = model_rsq, shared_a_pval = model_pval, analysis, mutation_count)

# shared_b = the shared mutations found in the recurrent tumor.
shared_b = read_tsv("/Volumes/verhaak-lab/GLASS-analysis/results/neutral_evolution/neutral-evolution-mutect2/tabular_results/neutral_results_shared_b.txt")
shared_b = shared_b %>% 
  mutate(case_barcode = substr(tumor_barcode, 1, 12)) %>% 
  left_join(shared_b_mutation_counts, by="case_barcode") %>% 
  mutate(analysis = "shared_b") %>% 
  filter(ploidy < 3, mutation_count >= 12) %>% 
  select(tumor_pair_barcode, shared_b_rsq = model_rsq, shared_b_pval = model_pval, analysis, mutation_count)

# Mutations private to the recurrent tumor.
recurrence = read_tsv("/Volumes/verhaak-lab/GLASS-analysis/results/neutral_evolution/neutral-evolution-mutect2/tabular_results/neutral_results_recur.txt")
recurrence = recurrence %>% 
  mutate(case_barcode = substr(tumor_barcode, 1, 12)) %>% 
  left_join(recurrence_mutation_counts, by="case_barcode") %>% 
  mutate(analysis = "recurrence") %>% 
  filter(ploidy < 3, mutation_count >= 12) %>% 
  select(tumor_pair_barcode, recur_rsq = model_rsq, recur_pval = model_pval, analysis, mutation_count)

# Derive case-level subtype so that each case is only represented once. Samples do not have subtype switch.
case_level_subtype = surgeries %>% 
  select(case_barcode, idh_codel_subtype) %>% 
  distinct() %>% 
  filter(!is.na(idh_codel_subtype)) 

# Combine mutations across the two time points.
all_samples = primary %>% 
  left_join(shared_a, by="tumor_pair_barcode") %>% 
  left_join(shared_b, by="tumor_pair_barcode") %>% 
  left_join(recurrence, by="tumor_pair_barcode") %>% 
  na.omit() %>% 
  mutate(primary_evolution = ifelse(primary_rsq> 0.98 & primary_pval > 0.05, "N", "S"), 
         shared_a_evolution = ifelse(shared_a_rsq> 0.98 & shared_a_pval > 0.05, "N", "S"),
         shared_b_evolution = ifelse(shared_b_rsq> 0.98 & shared_b_pval > 0.05, "N", "S"), 
         shared_evolution = ifelse(((shared_a_rsq+shared_b_rsq)/2)> 0.98, "N", "S"), 
         recurrence_evolution = ifelse(recur_rsq> 0.98 & recur_pval > 0.05, "N", "S"),
         evolution_mode = paste(primary_evolution, shared_evolution, recurrence_evolution, sep="-"))

# Assessment of missingness if na.omit is commented out.
sum(is.na(all_samples$primary_rsq))
sum(is.na(all_samples$shared_a_rsq))
sum(is.na(all_samples$shared_b_rsq))
sum(is.na(all_samples$recur_rsq))

# Quick look at neutrality estimates.
table(all_samples$primary_evolution)
table(all_samples$shared_a_evolution)
table(all_samples$shared_b_evolution)
table(all_samples$recurrence_evolution)


# Combine all of the R^2 information and gather into plot format.
plot_rsq = all_samples %>% 
  select(tumor_pair_barcode, primary_rsq, shared_a_rsq, shared_b_rsq, recur_rsq) %>% 
  gather("rsq_time", "rsq_value", c(primary_rsq, shared_a_rsq, shared_b_rsq, recur_rsq), -tumor_pair_barcode)
plot_rsq$rsq_time = factor(plot_rsq$rsq_time, levels = c("primary_rsq", "shared_a_rsq", "recur_rsq", "shared_b_rsq"))

# Manual input of the "primary" and "recurrence" state for plotting purposes.
plot_rsq$timepoint = NA
plot_rsq$timepoint = ifelse(plot_rsq$rsq_time == "primary_rsq", "primary", ifelse(plot_rsq$rsq_time=="shared_a_rsq", "primary", "recurrence"))

# Ladder plot for paired R-squared values based on neutrality testR.
ggplot(plot_rsq, aes(x = rsq_time, y = rsq_value, group = tumor_pair_barcode)) +
  geom_line(linetype="solid", size=1, alpha = 0.1) + 
  geom_point(color="black", size=2) + theme_bw() + ylab("R-squared values") + xlab("") + geom_hline(yintercept = 0.98, linetype="dotted") +
  facet_grid(~timepoint, scales = "free") + ggtitle("Mutect2, min 12 subclonal variants, Ploidy < 3 filter (n=108)")


# 1. Link tumor-pairs barcode with individual barcodes.
neutral_results_annot = all_samples %>% 
  left_join(tumor_pairs, by = "tumor_pair_barcode") %>% 
  left_join(case_level_subtype, by="case_barcode") %>% 
  select(-starts_with("analysis")) %>% 
  left_join(mut_freq, by=c("tumor_barcode_a"="aliquot_barcode")) %>% 
  left_join(mut_freq, by=c("tumor_barcode_b"="aliquot_barcode")) %>% 
  select(tumor_pair_barcode:idh_codel_subtype, mutation_count_a = mutation_count.x.x.x, tumor_a_mut_freq = coverage_adj_mut_freq.x ,mutation_count_b = mutation_count.y.y.y, tumor_b_mut_freq = coverage_adj_mut_freq.y) %>% 
  mutate(delta_mut_freq = tumor_b_mut_freq-tumor_a_mut_freq,
         hypermutant_recurrence = ifelse(tumor_b_mut_freq > 10.20, "1", "0"),
         sample_barcode_a = substr(tumor_barcode_a, 1, 15)) %>% 
  left_join(drivers, by= "tumor_pair_barcode") %>% 
  left_join(surgeries, by = c("sample_barcode_a"="sample_barcode")) 

# 68 unique cases
length(unique(neutral_results_annot$case_barcode))

# What about neutral evolution by subtype?
table(neutral_results_annot$evolution_mode, neutral_results_annot$idh_codel_subtype.y)
table(neutral_results_annot$evolution_mode, neutral_results_annot$hypermutant_recurrence)
table(neutral_results_annot$evolution_mode, neutral_results_annot$surgery_extent_of_resection)
# INTERESTING....
table(neutral_results_annot$evolution_mode, neutral_results_annot$treatment_tmz)
table(neutral_results_annot$evolution_mode, neutral_results_annot$driver_status)



################################
## Group-wise neutralitytestR ##
################################
# Set thresholds of > 0.5 purity AND < 3 ploidy.

# Primary.
glass_vaf_primary_filtered = glass_vaf_primary %>% 
  left_join(titan_info, by =c("tumor_barcode_a"="tumor_barcode")) %>% 
  left_join(titan_info, by =c("tumor_barcode_b"="tumor_barcode")) %>% 
  left_join(mut_freq, by=c("tumor_barcode_b"="aliquot_barcode")) %>% 
  mutate(hypermutant_recurrence = ifelse(coverage_adj_mut_freq > 10.2, "1", "0")) %>% 
  filter(purity.x > 0.5, purity.y > 0.5, ploidy.x < 3, ploidy.y < 3) 
  
# Fit neutral model to the data.
primary_hyper = glass_vaf_primary_filtered %>% filter(hypermutant_recurrence == 1)
primary_nonhyper = glass_vaf_primary_filtered %>% filter(hypermutant_recurrence == 0)
out_primary_1 = neutralitytest(primary_hyper$vaf_a, fmin = 0.1, fmax = 0.25)
out_primary_2 = neutralitytest(primary_nonhyper$vaf_a, fmin = 0.1, fmax = 0.25)
out_primary_1 # return statistics
out_primary_2
plot_all(out_primary_1)
plot_all(out_primary_2)


# Recurrence.
glass_vaf_recurrence_filtered = glass_vaf_recur %>% 
  left_join(titan_info, by =c("tumor_barcode_a"="tumor_barcode")) %>% 
  left_join(titan_info, by =c("tumor_barcode_b"="tumor_barcode")) %>%
  left_join(mut_freq, by=c("tumor_barcode_b"="aliquot_barcode")) %>% 
  mutate(hypermutant_recurrence = ifelse(coverage_adj_mut_freq > 10.2, "1", "0")) %>% 
  filter(purity.x > 0.5, purity.y > 0.5, ploidy.x < 3, ploidy.y < 3)

recur_hyper = glass_vaf_recurrence_filtered %>% filter(hypermutant_recurrence == 1)
recur_nonhyper = glass_vaf_recurrence_filtered %>% filter(hypermutant_recurrence == 0)
out_recur_1 = neutralitytest(recur_hyper$vaf_b, fmin = 0.1, fmax = 0.25)
out_recur_2 = neutralitytest(recur_nonhyper$vaf_b, fmin = 0.1, fmax = 0.25)
out_recur_1
out_recur_2 # return statistics
plot_all(out_recur_1)
plot_all(out_recur_2)

# SHARED
glass_vaf_shared_a_filtered = glass_vaf_shared_a %>% 
  left_join(titan_info, by =c("tumor_barcode_a"="tumor_barcode")) %>% 
  left_join(titan_info, by =c("tumor_barcode_b"="tumor_barcode")) %>% 
  left_join(mut_freq, by=c("tumor_barcode_b"="aliquot_barcode")) %>% 
  mutate(hypermutant_recurrence = ifelse(coverage_adj_mut_freq > 10.2, "1", "0")) %>% 
  filter(purity.x > 0.5, purity.y > 0.5, ploidy.x < 3, ploidy.y < 3)

shared_a_hyper = glass_vaf_shared_a_filtered %>% filter(hypermutant_recurrence == 1)
shared_a_nonhyper = glass_vaf_shared_a_filtered %>% filter(hypermutant_recurrence == 0)
out_shared_a_1 = neutralitytest(shared_a_hyper$vaf_b, fmin = 0.1, fmax = 0.25)
out_shared_a_2 = neutralitytest(shared_a_nonhyper$vaf_b, fmin = 0.1, fmax = 0.25)
out_shared_a_1
out_shared_a_2 # return statistics
plot_all(out_shared_a_1)
plot_all(out_shared_a_2)

# Shared_b
glass_vaf_shared_b_filtered = glass_vaf_shared_b %>% 
  left_join(titan_info, by =c("tumor_barcode_a"="tumor_barcode")) %>% 
  left_join(titan_info, by =c("tumor_barcode_b"="tumor_barcode")) %>%
  left_join(mut_freq, by=c("tumor_barcode_b"="aliquot_barcode")) %>% 
  mutate(hypermutant_recurrence = ifelse(coverage_adj_mut_freq > 10.2, "1", "0")) %>% 
  filter(purity.x > 0.5, purity.y > 0.5, ploidy.x < 3, ploidy.y < 3)

shared_b_hyper = glass_vaf_shared_b_filtered %>% filter(hypermutant_recurrence == 1)
shared_b_nonhyper = glass_vaf_shared_b_filtered %>% filter(hypermutant_recurrence == 0)
out_shared_b_1 = neutralitytest(shared_b_hyper$vaf_b, fmin = 0.1, fmax = 0.25)
out_shared_b_2 = neutralitytest(shared_b_nonhyper$vaf_b, fmin = 0.1, fmax = 0.25)
out_shared_b_1
out_shared_b_2 # return statistics
plot_all(out_shared_b_1)
plot_all(out_shared_b_2)

