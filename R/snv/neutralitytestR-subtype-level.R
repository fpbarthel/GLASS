##################################################
# NeutralitytestR results at the subtype and mutation status-level (private/shared)
# Updated: 2019.01.29
# Author: Kevin J.
##################################################

# Working directory for this analysis in the GLASS-analysis project. 
mybasedir = "/Volumes/verhaak-lab/GLASS-analysis/"
setwd(mybasedir)

##################################################
# Necessary packages:
library(tidyverse)
library(DBI)
library(neutralitytestr)
library(survminer)
library(survival)
require(alluvial)
library(ggExtra)
library(EnvStats)
library(gridExtra)

##################################################
# Establish connection with database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

# Load additional tables from the database.
pairs = dbReadTable(con,  Id(schema="analysis", table="pairs"))
tumor_pairs = dbReadTable(con,  Id(schema="analysis", table="tumor_pairs"))
aliquots = dbReadTable(con,  Id(schema="biospecimen", table="aliquots"))
surgeries = dbReadTable(con,  Id(schema="clinical", table="surgeries"))
cases = dbReadTable(con,  Id(schema="clinical", table="cases"))
titan_param = dbReadTable(con,  Id(schema="analysis", table="titan_params"))
mut_freq = dbReadTable(con,  Id(schema="analysis", table="mutation_freq"))
aliquot_neutrality = dbReadTable(con,  Id(schema="analysis", table="neutrality_aliquots"))
clinal_tumor_pairs = dbReadTable(con,  Id(schema="clinical", table="clinical_by_tumor_pair"))  

# These tables **MAY** change, especially the driver table.
silver_set = dbGetQuery(con, "SELECT * FROM analysis.silver_set")
all_drivers = dbGetQuery(con, "SELECT * FROM analysis.driver_status")

# Reduce the query so that only subtype information is available.
silver_subtype = all_drivers %>% 
  inner_join(clinal_tumor_pairs, by="tumor_pair_barcode") %>% 
  select(tumor_pair_barcode, idh_codel_subtype, hypermutator_status)

# Combine with pairs file to access "tumor_barcode".
titan_info = titan_param %>% 
  left_join(pairs, by="pair_barcode")

#
neutrality_input_aliquot_mutect2 = read_file("sql/neutralitytestr-input-aliquot-level.sql")
glass_single_vaf <- dbGetQuery(con, neutrality_input_aliquot_mutect2)


################################
## Group-wise neutralitytestR ##
################################
####### FILTER ########
# Floris updated the neutralitytestR input query to retrieve all possible tumor pairs AND designate mutation status.
neutrality_input_mutect2  = read_file("sql/neutrality-testr-input-mutect2.sql")


# Currently, there is a coverage filter requirement (30X in both tumor_a and tumor_b) PLUS it is required that both tumors have Mutect2 calls for shared variants.
glass_vaf_raw <- dbGetQuery(con, neutrality_input_mutect2)

# Add some additional descriptors (longitudinal vs. multisector) and filters.
glass_vaf = glass_vaf_raw %>% 
  mutate(tumor_a_type = substr(tumor_barcode_a, 14, 15),
         tumor_b_type = substr(tumor_barcode_b, 14, 15), 
         # Need to be able to merge with the other data tables.
         tumor_pair_barcode = paste(substr(tumor_barcode_a, 1, 18), substr(tumor_barcode_b, 14, 19), substr(tumor_barcode_a, 21, 23), sep="-"),
         comparison = paste(tumor_a_type, tumor_b_type, sep="-"),
         comparison_type = ifelse(tumor_a_type==tumor_b_type, "multisector", "longitudinal")) %>% 
  # Remove North Sydney's metastasis from the data set.
  filter(!comparison%in%c("R1-M1", "R2-M1", "TP-M1")) %>% 
  # For now, only consider the longitudinal samples.
  filter(comparison_type == "longitudinal")

#######################################
## Silver set
## Separated into two different groups for Hypermutators
## and non-hypermutators
#######################################
## Set thresholds of > 0.5 purity AND < 3 ploidy for inclusion.
glass_vaf_subtype = glass_vaf %>% 
  select(tumor_pair_barcode, comparison, comparison_type, chrom:status) %>% 
  inner_join(silver_set, by="tumor_pair_barcode") %>% 
  left_join(silver_subtype, by ="tumor_pair_barcode") %>% 
  select(-tumor_barcode_a, -tumor_barcode_b) %>% 
  left_join(clinal_tumor_pairs, by ="case_barcode") %>%
  left_join(titan_info, by =c("tumor_barcode_a"="tumor_barcode")) %>% 
  left_join(titan_info, by =c("tumor_barcode_b"="tumor_barcode")) %>% 
  filter(purity.x > 0.5, purity.y > 0.5, ploidy.x < 3, ploidy.y < 3) 

## Create a 3 x 4 plot to represent the VAF across subtypes.
plot_vaf = glass_vaf_subtype %>% 
  filter(hypermutator_status==0) %>% 
  gather(mutation_subtype, vaf, c(vaf_a, vaf_b), -status) %>% 
  mutate(var_filter = paste(mutation_subtype, status, sep="_")) %>% 
  filter(var_filter %in% c("vaf_a_P", "vaf_a_S", "vaf_b_S", "vaf_b_R")) %>% 
  mutate(var_filter_revised = recode(var_filter, "vaf_a_P" = "private P", "vaf_a_S" = "shared P", "vaf_b_S" = "shared R", "vaf_b_R" = "private R"),
         idh_codel_subtype_revised = recode(idh_codel_subtype, "IDHwt_noncodel" = "IDHwt", "IDHmut_noncodel"= "IDHmut noncodel", "IDHmut_codel" = "IDHmut codel"))
plot_vaf$var_filter_revised = factor(plot_vaf$var_filter_revised, levels = c("private P", "shared P", "shared R", "private R"))

# Hypermutators removed. 167 samples (removed samples with low purity and high ploidy).
n_distinct(plot_vaf$case_barcode)
ggplot(plot_vaf, aes(x=vaf)) + geom_histogram(binwidth = 0.01) + theme_bw() + facet_grid(idh_codel_subtype_revised~var_filter_revised, scales="free") +
  xlab("Allele frequency (f)") + ylab("Number of mutations") + ggtitle("Non-hypermutators")


###########################
# Neutrality models
###########################
# Primary.
codel_primary_nonhyper = glass_vaf_subtype %>% 
  filter(status =="P", idh_codel_subtype == "IDHmut_codel", hypermutator_status==0) 
out_primary_1 = neutralitytest(codel_primary_nonhyper$vaf_a, fmin = 0.1, fmax = 0.25)
out_primary_1 # return statistics
noncodel_primary_nonhyper = glass_vaf_subtype %>% 
  filter(status =="P", idh_codel_subtype == "IDHmut_noncodel", hypermutator_status==0) 
out_primary_2 = neutralitytest(noncodel_primary_nonhyper$vaf_a, fmin = 0.1, fmax = 0.25)
out_primary_2 # return statistics
wt_primary_nonhyper = glass_vaf_subtype %>% 
  filter(status =="P", idh_codel_subtype == "IDHwt_noncodel", hypermutator_status==0) 
out_primary_3 = neutralitytest(wt_primary_nonhyper$vaf_a, fmin = 0.1, fmax = 0.25)
out_primary_3 # return statistics

# Recurrence.
codel_recur_nonhyper = glass_vaf_subtype %>% 
  filter(status =="R", idh_codel_subtype == "IDHmut_codel", hypermutator_status==0) 
out_recur_1 = neutralitytest(codel_recur_nonhyper$vaf_b, fmin = 0.1, fmax = 0.25)
out_recur_1 # return statistics
noncodel_recur_nonhyper = glass_vaf_subtype %>% 
  filter(status =="R", idh_codel_subtype == "IDHmut_noncodel", hypermutator_status==0) 
out_recur_2 = neutralitytest(noncodel_recur_nonhyper$vaf_b, fmin = 0.1, fmax = 0.25)
out_recur_2 # return statistics
wt_recur_nonhyper = glass_vaf_subtype %>% 
  filter(status =="R", idh_codel_subtype == "IDHwt_noncodel", hypermutator_status==0) 
out_recur_3 = neutralitytest(wt_recur_nonhyper$vaf_b, fmin = 0.1, fmax = 0.25)
out_recur_3 # return statistics

# SHARED - tumor_a
codel_shared_a_nonhyper = glass_vaf_subtype %>% 
  filter(status =="S", idh_codel_subtype == "IDHmut_codel", hypermutator_status==0) 
out_shared_a_1 = neutralitytest(codel_shared_a_nonhyper$vaf_a, fmin = 0.1, fmax = 0.25)
out_shared_a_1 # return statistics
noncodel_shared_a_nonhyper = glass_vaf_subtype %>% 
  filter(status =="S", idh_codel_subtype == "IDHmut_noncodel", hypermutator_status==0) 
out_shared_a_2 = neutralitytest(noncodel_shared_a_nonhyper$vaf_a, fmin = 0.1, fmax = 0.25)
out_shared_a_2 # return statistics
wt_shared_a_nonhyper = glass_vaf_subtype %>% 
  filter(status =="S", idh_codel_subtype == "IDHwt_noncodel", hypermutator_status==0) 
out_shared_a_3 = neutralitytest(wt_shared_a_nonhyper$vaf_a, fmin = 0.1, fmax = 0.25)
out_shared_a_3 # return statistics

# SHARED - tumor_b
codel_shared_b_nonhyper = glass_vaf_subtype %>% 
  filter(status =="S", idh_codel_subtype == "IDHmut_codel", hypermutator_status==0) 
out_shared_b_1 = neutralitytest(codel_shared_b_nonhyper$vaf_b, fmin = 0.1, fmax = 0.25)
out_shared_b_1 # return statistics
vnoncodel_shared_b_nonhyper = glass_vaf_subtype %>% 
  filter(status =="S", idh_codel_subtype == "IDHmut_noncodel", hypermutator_status==0) 
out_shared_b_2 = neutralitytest(noncodel_shared_a_nonhyper$vaf_b, fmin = 0.1, fmax = 0.25)
out_shared_b_2 # return statistics
wt_shared_b_nonhyper = glass_vaf_subtype %>% 
  filter(status =="S", idh_codel_subtype == "IDHwt_noncodel", hypermutator_status==0) 
out_shared_b_3 = neutralitytest(wt_shared_a_nonhyper$vaf_b, fmin = 0.1, fmax = 0.25)
out_shared_b_3 # return statistics

# Plot the cumulative distribution of the data as well as the best fit linear model line.
q1 = lsq_plot(out_primary_1) + ggtitle("IDHmut codel") + xlab("") 
q2 = lsq_plot(out_primary_2) + ggtitle("IDHmut noncodel") + ylab("") 
q3 = lsq_plot(out_primary_3) + ggtitle("IDHwt") + ylab("") + xlab("")
q4 = lsq_plot(out_recur_1)  + ggtitle("IDHmut codel") + xlab("")
q5 = lsq_plot(out_recur_2) + ggtitle("IDHmut noncodel") + ylab("")
q6 =lsq_plot(out_recur_3) + ggtitle("IDHwt") + ylab("") + xlab("")
q7 = lsq_plot(out_shared_a_1)  + ggtitle("IDHmut codel") + xlab("")
q8 = lsq_plot(out_shared_a_2)  + ggtitle("IDHmut noncodel")  + ylab("")
q9 = lsq_plot(out_shared_a_3) + ggtitle("IDHwt") + ylab("") + xlab("")
q10 = lsq_plot(out_shared_b_1)  + ggtitle("IDHmut noncodel") + xlab("") 
q11 = lsq_plot(out_shared_b_2)  + ggtitle("IDHmut noncodel") + ylab("") 
q12 = lsq_plot(out_shared_b_3) + ggtitle("IDHwt") + ylab("") + xlab("")

# Difficult to reproduce the faceted layout of the VAF plots.
grid.arrange(q1, q2, q3, ncol=3)
grid.arrange(q4, q5, q6, ncol=3)
grid.arrange(q7, q8, q9, ncol=3)
grid.arrange(q10, q11, q12, ncol=3)

# Plot the normalized cumulative distributions.
p1 = normalized_plot(out_primary_1) + ggtitle("IDHmut codel") + xlab("") 
p2 = normalized_plot(out_primary_2) + ggtitle("IDHmut noncodel") + ylab("") 
p3 = normalized_plot(out_primary_3) + ggtitle("IDHwt") + ylab("") + xlab("")
p4 = normalized_plot(out_recur_1)  + ggtitle("IDHmut codel") + xlab("")
p5 = normalized_plot(out_recur_2) + ggtitle("IDHmut noncodel") + ylab("")
p6 = normalized_plot(out_recur_3) + ggtitle("IDHwt") + ylab("") + xlab("")
p7 = normalized_plot(out_shared_a_1)  + ggtitle("IDHmut codel") + xlab("")
p8 = normalized_plot(out_shared_a_2)  + ggtitle("IDHmut noncodel")  + ylab("")
p9 = normalized_plot(out_shared_a_3) + ggtitle("IDHwt") + ylab("") + xlab("")
p10 = normalized_plot(out_shared_b_1) + ggtitle("IDHmut codel") + xlab("")
p11 = normalized_plot(out_shared_b_2)  +ggtitle("IDHmut noncodel") + ylab("") 
p12 = normalized_plot(out_shared_b_3) + ggtitle("IDHwt") + ylab("") + xlab("")

# Plot the normalized values here:
grid.arrange(p1, p2, p3, ncol=3)
grid.arrange(p4, p5, p6, ncol=3)
grid.arrange(p7, p8, p9, ncol=3)
grid.arrange(p10, p11, p12, ncol=3)

#######################
# Hypermutation
#######################
### Create a 3 x 4 plot to represent the VAF across subtypes.
plot_vaf_hyper = glass_vaf_subtype %>% 
  filter(hypermutator_status==1) %>% 
  gather(mutation_subtype, vaf, c(vaf_a, vaf_b), -status) %>% 
  mutate(var_filter = paste(mutation_subtype, status, sep="_")) %>% 
  filter(var_filter %in% c("vaf_a_P", "vaf_a_S", "vaf_b_S", "vaf_b_R")) %>% 
  mutate(var_filter_revised = recode(var_filter, "vaf_a_P" = "private P", "vaf_a_S" = "shared P", "vaf_b_S" = "shared R", "vaf_b_R" = "private R"),
         idh_codel_subtype_revised = recode(idh_codel_subtype, "IDHwt_noncodel" = "IDHwt", "IDHmut_noncodel"= "IDHmut noncodel", "IDHmut_codel" = "IDHmut codel"))
plot_vaf_hyper$var_filter_revised = factor(plot_vaf_hyper$var_filter_revised, levels = c("private P", "shared P", "shared R", "private R"))

# Hypermutators, n = 29 samples.
n_distinct(plot_vaf_hyper$case_barcode)
ggplot(plot_vaf_hyper, aes(x=vaf)) + geom_histogram(binwidth = 0.01) + theme_bw() + facet_wrap(idh_codel_subtype_revised~var_filter_revised, scales="free") +
  xlab("Allele frequency (f)") + ylab("Number of mutations") + ggtitle("Hypermutators")

###########################
# Neutrality models **HYPERMUTATORS**
###########################
# Primary.
codel_primary_hyper = glass_vaf_subtype %>% 
  filter(status =="P", idh_codel_subtype == "IDHmut_codel", hypermutator_status==1) 
out_primary_1 = neutralitytest(codel_primary_hyper$vaf_a, fmin = 0.1, fmax = 0.25)
out_primary_1 # return statistics
noncodel_primary_hyper = glass_vaf_subtype %>% 
  filter(status =="P", idh_codel_subtype == "IDHmut_noncodel", hypermutator_status==1) 
out_primary_2 = neutralitytest(noncodel_primary_hyper$vaf_a, fmin = 0.1, fmax = 0.25)
out_primary_2 # return statistics
wt_primary_hyper = glass_vaf_subtype %>% 
  filter(status =="P", idh_codel_subtype == "IDHwt_noncodel", hypermutator_status==1) 
out_primary_3 = neutralitytest(wt_primary_hyper$vaf_a, fmin = 0.1, fmax = 0.25)
out_primary_3 # return statistics

# Recurrence.
codel_recur_hyper = glass_vaf_subtype %>% 
  filter(status =="R", idh_codel_subtype == "IDHmut_codel", hypermutator_status==1) 
out_recur_1 = neutralitytest(codel_recur_hyper$vaf_b, fmin = 0.1, fmax = 0.25)
out_recur_1 # return statistics
noncodel_recur_hyper = glass_vaf_subtype %>% 
  filter(status =="R", idh_codel_subtype == "IDHmut_noncodel", hypermutator_status==1) 
out_recur_2 = neutralitytest(noncodel_recur_hyper$vaf_b, fmin = 0.1, fmax = 0.25)
out_recur_2 # return statistics
wt_recur_hyper = glass_vaf_subtype %>% 
  filter(status =="R", idh_codel_subtype == "IDHwt_noncodel", hypermutator_status==1) 
out_recur_3 = neutralitytest(wt_recur_hyper$vaf_b, fmin = 0.1, fmax = 0.25)
out_recur_3 # return statistics

# SHARED - tumor_a
codel_shared_a_hyper = glass_vaf_subtype %>% 
  filter(status =="S", idh_codel_subtype == "IDHmut_codel", hypermutator_status==1) 
out_shared_a_1 = neutralitytest(codel_shared_a_hyper$vaf_a, fmin = 0.1, fmax = 0.25)
out_shared_a_1 # return statistics
noncodel_shared_a_hyper = glass_vaf_subtype %>% 
  filter(status =="S", idh_codel_subtype == "IDHmut_noncodel", hypermutator_status==1) 
out_shared_a_2 = neutralitytest(noncodel_shared_a_hyper$vaf_a, fmin = 0.1, fmax = 0.25)
out_shared_a_2 # return statistics
wt_shared_a_hyper = glass_vaf_subtype %>% 
  filter(status =="S", idh_codel_subtype == "IDHwt_noncodel", hypermutator_status==1) 
out_shared_a_3 = neutralitytest(wt_shared_a_hyper$vaf_a, fmin = 0.1, fmax = 0.25)
out_shared_a_3 # return statistics

# SHARED - tumor_b
codel_shared_b_hyper = glass_vaf_subtype %>% 
  filter(status =="S", idh_codel_subtype == "IDHmut_codel", hypermutator_status==1) 
out_shared_b_1 = neutralitytest(codel_shared_b_hyper$vaf_b, fmin = 0.1, fmax = 0.25)
out_shared_b_1 # return statistics
vnoncodel_shared_b_hyper = glass_vaf_subtype %>% 
  filter(status =="S", idh_codel_subtype == "IDHmut_noncodel", hypermutator_status==1) 
out_shared_b_2 = neutralitytest(vnoncodel_shared_b_hyper$vaf_b, fmin = 0.1, fmax = 0.25)
out_shared_b_2 # return statistics
wt_shared_b_hyper = glass_vaf_subtype %>% 
  filter(status =="S", idh_codel_subtype == "IDHwt_noncodel", hypermutator_status==1) 
out_shared_b_3 = neutralitytest(wt_shared_a_hyper$vaf_b, fmin = 0.1, fmax = 0.25)
out_shared_b_3 # return statistics

# Plot the cumulative distribution of the data as well as the best fit linear model line.
q1 = lsq_plot(out_primary_1) + ggtitle("IDHmut codel") + xlab("") 
q2 = lsq_plot(out_primary_2) + ggtitle("IDHmut noncodel") + ylab("") 
q3 = lsq_plot(out_primary_3) + ggtitle("IDHwt") + ylab("") + xlab("")
q4 = lsq_plot(out_recur_1)  + ggtitle("IDHmut codel") + xlab("")
q5 = lsq_plot(out_recur_2) + ggtitle("IDHmut noncodel") + ylab("")
q6 =lsq_plot(out_recur_3) + ggtitle("IDHwt") + ylab("") + xlab("")
q7 = lsq_plot(out_shared_a_1)  + ggtitle("IDHmut codel") + xlab("")
q8 = lsq_plot(out_shared_a_2)  + ggtitle("IDHmut noncodel")  + ylab("")
q9 = lsq_plot(out_shared_a_3) + ggtitle("IDHwt") + ylab("") + xlab("")
q10 = lsq_plot(out_shared_b_1)  + ggtitle("IDHmut noncodel") + ylab("") 
q11 = lsq_plot(out_shared_b_2)  + ggtitle("IDHmut noncodel") + ylab("") 
q12 = lsq_plot(out_shared_b_3) + ggtitle("IDHwt") + ylab("") + xlab("")

# Difficult to reproduce the faceted layout of the VAF plots.
grid.arrange(q1, q2, q3, ncol=3)
grid.arrange(q4, q5, q6, ncol=3)
grid.arrange(q7, q8, q9, ncol=3)
grid.arrange(q10, q11, q12, ncol=3)

# Plot the normalized cumulative distributions.
p1 = normalized_plot(out_primary_1) + ggtitle("IDHmut codel") + xlab("") 
p2 = normalized_plot(out_primary_2) + ggtitle("IDHmut noncodel") + ylab("") 
p3 = normalized_plot(out_primary_3) + ggtitle("IDHwt") + ylab("") + xlab("")
p4 = normalized_plot(out_recur_1)  + ggtitle("IDHmut codel") + xlab("")
p5 = normalized_plot(out_recur_2) + ggtitle("IDHmut noncodel") + ylab("")
p6 = normalized_plot(out_recur_3) + ggtitle("IDHwt") + ylab("") + xlab("")
p7 = normalized_plot(out_shared_a_1)  + ggtitle("IDHmut codel") + xlab("")
p8 = normalized_plot(out_shared_a_2)  + ggtitle("IDHmut noncodel")  + ylab("")
p9 = normalized_plot(out_shared_a_3) + ggtitle("IDHwt") + ylab("") + xlab("")
p10 = normalized_plot(out_shared_b_1) + ggtitle("IDHmut codel") + xlab("")
p11 = normalized_plot(out_shared_b_2)  +ggtitle("IDHmut noncodel") + ylab("") 
p12 = normalized_plot(out_shared_b_3) + ggtitle("IDHwt") + ylab("") + xlab("")

# Plot the normalized values here:
grid.arrange(p1, p2, p3, ncol=3)
grid.arrange(p4, p5, p6, ncol=3)
grid.arrange(p7, p8, p9, ncol=3)
grid.arrange(p10, p11, p12, ncol=3)
