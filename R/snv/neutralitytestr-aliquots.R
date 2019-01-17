##############################################
# NeutralitytestR applied to each aliquot
# Updated: 2019.01.16
# Author: Kevin J.
##################################################

# Working directory for this analysis in the GLASS-analysis project. 
mybasedir = "/Volumes/verhaak-lab/GLASS-analysis/"
setwd(mybasedir)

#######################################################
# Necessary packages:
library(tidyverse)
library(DBI)
library(neutralitytestr)
library(survminer)
library(survival)
library(ggExtra)
library(EnvStats)

#######################################################
# Establish connection with Floris' database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

# Load additional tables from the database.
pairs = dbReadTable(con,  Id(schema="analysis",table="pairs"))
tumor_pairs = dbReadTable(con,  Id(schema="analysis", table="tumor_pairs"))
aliquots = dbReadTable(con,  Id(schema="biospecimen", table="aliquots"))
surgeries = dbReadTable(con,  Id(schema="clinical", table="surgeries"))
cases = dbReadTable(con,  Id(schema="clinical", table="cases"))
titan_param = dbReadTable(con,  Id(schema="analysis",table="titan_params"))
mut_freq = dbReadTable(con,  Id(schema="analysis",table="mutation_freq"))
neutrality_aliquots = dbReadTable(con,  Id(schema="analysis",table="neutrality_aliquots"))
neutrality_tumor_pairs = dbReadTable(con,  Id(schema="analysis",table="neutrality_tumor_pairs"))

# These tables **MAY** change, especially the driver table.
clinal_tumor_pairs = dbGetQuery(con,"SELECT * FROM analysis.clinical_by_tumor_pair")
drivers = dbGetQuery(con, "SELECT * FROM analysis.driver_status")
silver_set = dbGetQuery(con, "SELECT * FROM analysis.silver_set")

# Combine with pairs file to access "tumor_barcode".
titan_info = titan_param %>% 
  left_join(pairs, by="pair_barcode")

#######################################################
# Traditional neutralitytestR analysis of all variants in each tumor
#######################################################
# To mirror what Georgette is doing with the SubclonalSelection script, 
# we will include all Mutect2 mutations within an aliquot.
neutrality_input_aliquot_mutect2 = read_file("sql/neutralitytestr-input-aliquot-level.sql")
glass_single_vaf <- dbGetQuery(con, neutrality_input_aliquot_mutect2)

# 514 distinct aliquots.
n_distinct(glass_single_vaf$aliquot_barcode)

# Determine the number of mutations that are considered subclonal in each tumor.
aliquot_mutation_counts = glass_single_vaf %>% 
  filter(vaf >= 0.1 & vaf <= 0.25) %>% 
  group_by(aliquot_barcode) %>% 
  summarize(subclonal_mut = n())

# Store results from per-level neutralitytestR.
aliquot_results = matrix(NA, nrow = length(unique(glass_single_vaf$tumor_pair_barcode)), ncol = 12)
colnames(aliquot_results) =  c("aliquot_barcode", "ploidy", "cellularity", "mutation_rate", "model_rsq","model_pval", "area_value", "area_pval","meanDist_value", "meanDist_pval", "Dk_value", "Dk_pval")
aliquot_results = as.data.frame(aliquot_results)

# Set directory for output.
setwd("/Volumes/verhaak-lab/GLASS-analysis/results/neutral_evolution/neutral-evolution-mutect2/aliquot-level/")

## Generate results for each tumor aliquot.
for ( i in 1:length(unique(glass_single_vaf$aliquot_barcode))){ 
  # Create subsetted data for this case.
  Yi = glass_single_vaf[glass_single_vaf$aliquot_barcode==unique(glass_single_vaf$aliquot_barcode)[i], ]
  xi = unique(Yi$aliquot_barcode)
  
  print(sprintf("Analyzing all mutations for: %s",  xi))
  # Store sample information.
  aliquot_results[i, "aliquot_barcode"] = xi
  aliquot_results[i, "ploidy"] = titan_info[titan_info$tumor_barcode==xi, "ploidy"]
  aliquot_results[i, "cellularity"] = titan_info[titan_info$tumor_barcode==xi, "purity"]
  
  # The model can sometimes fail. Skip over samples that error out.
  fit_a =  try(neutralitytest(Yi$vaf, ploidy = titan_info[titan_info$tumor_barcode==xi, "ploidy"], cellularity = titan_info[titan_info$tumor_barcode==xi, "purity"]), silent = T) 
  if(!inherits(fit_a, "try-error")) {
    aliquot_results[i, "mutation_rate"] = fit_a$mutation.rate
    aliquot_results[i, "model_rsq"] = fit_a$model$rsq
    aliquot_results[i, "model_pval"] = fit_a$model$pval
    aliquot_results[i, "area_value"] = fit_a$area$metric
    aliquot_results[i, "area_pval"] = fit_a$area$pval
    aliquot_results[i, "meanDist_value"] = fit_a$meanDist$metric
    aliquot_results[i, "meanDist_pval"] = fit_a$meanDist$pval
    aliquot_results[i, "Dk_value"] = fit_a$Dk$metric
    aliquot_results[i, "Dk_pval"] = fit_a$Dk$pval
    
    # Print out model plots to directory of interest.
    pdf(paste(unique(glass_single_vaf$aliquot_barcode)[i],"-_all",".pdf", sep=""), width=9, height=6)
    print(plot_all(fit_a))
    dev.off()
  }
}

# Combine data and examine the patterns.
aliquot_neutrality = aliquot_results %>% 
  left_join(aliquot_mutation_counts, by="aliquot_barcode") 

# Export to database:
dbWriteTable(con, Id(schema="analysis",table="neutrality_aliquots"), as.data.frame(aliquot_neutrality))