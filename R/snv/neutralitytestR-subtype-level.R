##################################################
# NeutralitytestR results at the subtype and mutation status-level (private/shared)
# Updated: 2019.05.06
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
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB2")

# Load additional tables.
silver_set = dbReadTable(con,  Id(schema="analysis", table="silver_set"))
gold_set = dbReadTable(con,  Id(schema="analysis", table="gold_set"))
subtypes = dbReadTable(con,  Id(schema="clinical", table="subtypes"))
variant_classifications = dbReadTable(con,  Id(schema="variants", table="variant_classifications"))

# Construct a table that provides subject-level information about clinical variables between two timepoints.
clinical_tumor_pairs_query = read_file("/sql/clinical/clinical-tumor-pairs-db2.sql")
clinical_tumor_pairs <- dbGetQuery(con, clinical_tumor_pairs_query)

################################
## Group-wise neutralitytestR ##
################################
# First, analyze the non-hypermutation in recurrent tumor.
# Currently, there is a coverage filter requirement (30X in both tumor_a and tumor_b) PLUS it is required that both tumors have Mutect2 calls for shared variants.
neutrality_input  = read_file("/sql/neutrality/neutralitytestr-subtype-resubmission.sql")
glass_vaf_cnv_unfiltered <- dbGetQuery(con, neutrality_input)

# Examine the VAFs for variants with 30X. Filter based on copy number status and hypermutation event.
vaf_cnv_silver = glass_vaf_cnv_unfiltered %>% 
  filter(tumor_pair_barcode %in% gold_set$tumor_pair_barcode) %>% 
  left_join(clinical_tumor_pairs, by=c("tumor_pair_barcode", "case_barcode", "tumor_barcode_a", "tumor_barcode_b")) %>% 
  left_join(subtypes, by="case_barcode") %>%
  left_join(variant_classifications, by=c("variant_classification"="variant_classification")) %>% 
  mutate(mut_class = ifelse(variant_classification_impact%in%c("HIGH", "MODERATE"), "nonsynonymous", ifelse(variant_classification_impact%in%c("LOW"), "synonymous", NA))) %>% 
  filter(fraction =="P" & cnv_call_a == 0 | fraction =="R" & cnv_call_b == 0 | fraction =="S" & cnv_call_a == 0 & cnv_call_b == 0) %>% 
  filter(hypermutator_status==0) 

# Put the data into a long format, separate by timepoint and fraction.
plot_vaf = vaf_cnv_silver %>% 
  select(tumor_pair_barcode, vaf_a = variant_allele_frequency_a, vaf_b = variant_allele_frequency_b, fraction, mut_class, idh_codel_subtype) %>% 
  gather(mutation_subtype, vaf, c(vaf_a, vaf_b), -fraction) %>% 
  mutate(var_filter = paste(mutation_subtype, fraction, sep="_")) %>% 
  filter(var_filter %in% c("vaf_a_P", "vaf_a_S", "vaf_b_S", "vaf_b_R")) %>% 
  mutate(var_filter_revised = recode(var_filter, "vaf_a_P" = "private P", "vaf_a_S" = "shared P", "vaf_b_S" = "shared R", "vaf_b_R" = "private R"),
         sample_barcode = substr(tumor_pair_barcode, 6, 12)) 
plot_vaf$var_filter_revised = factor(plot_vaf$var_filter_revised, levels = c("private P", "shared P", "shared R", "private R"))

# Hypermutators removed. 181 samples (probably removed samples with lots of CNVs that did not pass). [195 vs. 181 for silver_set vs. gold set].
n_distinct(plot_vaf$tumor_pair_barcode)
pdf(file = "/Users/johnsk/Documents/nonhypermutants-vaf-distribution-copy-neutral.pdf", height = 4, width = 8, bg = "transparent", useDingbats = FALSE)
ggplot(plot_vaf, aes(x=vaf)) + geom_histogram(binwidth = 0.01) + theme_bw() + facet_grid(idh_codel_subtype~var_filter_revised, scales="free") +
  xlab("Allele frequency (f)") + ylab("Number of mutations") + ggtitle("Non-hypermutators, copy neutral regions, all variant types (n=181)")
dev.off()

# Nonsynonymous mutations.
nonsyn_plot_vaf = plot_vaf %>% filter(mut_class == "nonsynonymous")
n_distinct(nonsyn_plot_vaf$tumor_pair_barcode)
pdf(file = "/Users/johnsk/Documents/nonhypermutants-vaf-distribution-copy-neutral-nonsynonymous.pdf", height = 4, width = 8, bg = "transparent", useDingbats = FALSE)
ggplot(nonsyn_plot_vaf, aes(x=vaf)) + geom_histogram(binwidth = 0.01) + theme_bw() + facet_grid(idh_codel_subtype~var_filter_revised, scales="free") +
  xlab("Allele frequency (f)") + ylab("Number of mutations") + ggtitle("Non-hypermutators, copy neutral regions, nonsynonymous mutations (n=173)")
dev.off()

# Synonymous mutations.
syn_plot_vaf = plot_vaf %>% filter(mut_class == "synonymous")
n_distinct(syn_plot_vaf$tumor_pair_barcode)
pdf(file = "/Users/johnsk/Documents/nonhypermutants-vaf-distribution-copy-neutral-synonymous.pdf", height = 4, width = 8, bg = "transparent", useDingbats = FALSE)
ggplot(syn_plot_vaf, aes(x=vaf)) + geom_histogram(binwidth = 0.01) + theme_bw() + facet_grid(idh_codel_subtype~var_filter_revised, scales="free") +
  xlab("Allele frequency (f)") + ylab("Number of mutations") + ggtitle("Non-hypermutators, copy neutral regions, synonymous mutations (n=168)")
dev.off()

## Going back to the non-subdivided data, inspect the cancer cell fraction for copy-neutral variants.
plot_ccf = vaf_cnv_silver %>% 
  select(tumor_pair_barcode, ccf_a = cellular_prevalence_a, ccf_b = cellular_prevalence_b, fraction, idh_codel_subtype) %>% 
  gather(mutation_subtype, ccf, c(ccf_a, ccf_b), -fraction) %>% 
  mutate(var_filter = paste(mutation_subtype, fraction, sep="_")) %>% 
  filter(var_filter %in% c("ccf_a_P", "ccf_a_S", "ccf_b_S", "ccf_b_R")) %>% 
  mutate(var_filter_revised = recode(var_filter, "ccf_a_P" = "private P", "ccf_a_S" = "shared P", "ccf_b_S" = "shared R", "ccf_b_R" = "private R"),
         sample_barcode = substr(tumor_pair_barcode, 6, 12))
plot_ccf$var_filter_revised = factor(plot_ccf$var_filter_revised, levels = c("private P", "shared P", "shared R", "private R"))

# Provides information about clonality (clonal mutations = CCF > 0.5). 
pdf(file = "/Users/johnsk/Documents/nonhypermutants-ccf-copy-neutral.pdf", height = 4, width = 8, bg = "transparent", useDingbats = FALSE)
ggplot(plot_ccf, aes(x=ccf)) + geom_histogram(binwidth = 0.01) + theme_bw() + facet_grid(idh_codel_subtype~var_filter_revised, scales="free") +
  xlab("Cancer Cell Fraction") + ylab("Number of mutations") + ggtitle("Non-hypermutators, copy neutral regions, all variant types (n=195)")
dev.off()

# Determine the proportion of ccf increases for shared variants.
plot_shared_ccf = vaf_cnv_silver %>% 
  filter(fraction == "S") %>% 
  mutate(ccf_diff = cellular_prevalence_b-cellular_prevalence_a,
         vaf_diff = variant_allele_frequency_b-variant_allele_frequency_a)
sum(plot_shared_ccf$ccf_diff > 0)/dim(plot_shared_ccf)[1] #  62% have increase in ccf

# Visualize the change in cancer cell fraction for shared variants.
pdf(file = "/Users/johnsk/Documents/nonhypermutants-shared-fraction-delta-ccf.pdf", height = 4, width = 8, bg = "transparent", useDingbats = FALSE)
ggplot(plot_shared_ccf, aes(x=fraction, y=ccf_diff)) + geom_violin() + theme_bw() + ylab("Delta CCF (recurrence-primary)") + xlab("Shared fraction") + facet_grid(~mut_class, scales="free")
dev.off()


###########################################
# Neutrality models for non-hypermutators
###########################################
# Primary.
codel_primary_nonhyper = vaf_cnv_silver %>% 
  filter(fraction =="P", idh_codel_subtype == "IDHmut-codel", hypermutator_status==0) 
out_primary_1 = neutralitytest(codel_primary_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
out_primary_1 # return statistics
noncodel_primary_nonhyper = vaf_cnv_silver %>% 
  filter(fraction =="P", idh_codel_subtype == "IDHmut-noncodel", hypermutator_status==0) 
out_primary_2 = neutralitytest(noncodel_primary_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
out_primary_2 # return statistics
wt_primary_nonhyper = vaf_cnv_silver %>% 
  filter(fraction =="P", idh_codel_subtype == "IDHwt", hypermutator_status==0) 
out_primary_3 = neutralitytest(wt_primary_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
out_primary_3 # return statistics

# Recurrence.
codel_recur_nonhyper = vaf_cnv_silver %>% 
  filter(fraction =="R", idh_codel_subtype == "IDHmut-codel", hypermutator_status==0) 
out_recur_1 = neutralitytest(codel_recur_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
out_recur_1 # return statistics
noncodel_recur_nonhyper = vaf_cnv_silver %>% 
  filter(fraction =="R", idh_codel_subtype == "IDHmut-noncodel", hypermutator_status==0) 
out_recur_2 = neutralitytest(noncodel_recur_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
out_recur_2 # return statistics
wt_recur_nonhyper = vaf_cnv_silver %>% 
  filter(fraction =="R", idh_codel_subtype == "IDHwt", hypermutator_status==0) 
out_recur_3 = neutralitytest(wt_recur_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
out_recur_3 # return statistics

# SHARED - tumor_a
codel_shared_a_nonhyper = vaf_cnv_silver %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHmut-codel", hypermutator_status==0) 
out_shared_a_1 = neutralitytest(codel_shared_a_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
out_shared_a_1 # return statistics
noncodel_shared_a_nonhyper = vaf_cnv_silver %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHmut-noncodel", hypermutator_status==0) 
out_shared_a_2 = neutralitytest(noncodel_shared_a_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
out_shared_a_2 # return statistics
wt_shared_a_nonhyper = vaf_cnv_silver %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHwt", hypermutator_status==0) 
out_shared_a_3 = neutralitytest(wt_shared_a_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
out_shared_a_3 # return statistics

# SHARED - tumor_b
codel_shared_b_nonhyper = vaf_cnv_silver %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHmut-codel", hypermutator_status==0) 
out_shared_b_1 = neutralitytest(codel_shared_b_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
out_shared_b_1 # return statistics
noncodel_shared_b_nonhyper = vaf_cnv_silver %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHmut-noncodel", hypermutator_status==0) 
out_shared_b_2 = neutralitytest(noncodel_shared_a_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
out_shared_b_2 # return statistics
wt_shared_b_nonhyper = vaf_cnv_silver %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHwt", hypermutator_status==0) 
out_shared_b_3 = neutralitytest(wt_shared_a_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
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

grid.arrange(q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, nrow=3, ncol=4)

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

##################################
## Hypermutator-specific analyses
###################################
# Select for hypermutation events.
vaf_cnv_silver_hyper = glass_vaf_cnv_unfiltered %>% 
  filter(tumor_pair_barcode %in% gold_set$tumor_pair_barcode) %>% 
  left_join(clinical_tumor_pairs, by=c("tumor_pair_barcode", "case_barcode", "tumor_barcode_a", "tumor_barcode_b")) %>% 
  left_join(subtypes, by="case_barcode") %>%
  left_join(variant_classifications, by=c("variant_classification"="variant_classification")) %>% 
  mutate(mut_class = ifelse(variant_classification_impact%in%c("HIGH", "MODERATE"), "nonsynonymous", ifelse(variant_classification_impact%in%c("LOW"), "synonymous", NA))) %>% 
  filter(fraction =="P" & cnv_call_a == 0 | fraction =="R" & cnv_call_b == 0 | fraction =="S" & cnv_call_a == 0 & cnv_call_b == 0) %>% 
  filter(hypermutator_status==1) 

# Create a vaf plot by timepoint and subtype.
plot_vaf_hyper = vaf_cnv_silver_hyper %>% 
  select(tumor_pair_barcode, vaf_a = variant_allele_frequency_a, vaf_b = variant_allele_frequency_b, fraction, idh_codel_subtype) %>% 
  gather(mutation_subtype, vaf, c(vaf_a, vaf_b), -fraction) %>% 
  mutate(var_filter = paste(mutation_subtype, fraction, sep="_")) %>% 
  filter(var_filter %in% c("vaf_a_P", "vaf_a_S", "vaf_b_S", "vaf_b_R")) %>% 
  mutate(var_filter_revised = recode(var_filter, "vaf_a_P" = "private P", "vaf_a_S" = "shared P", "vaf_b_S" = "shared R", "vaf_b_R" = "private R"),
         sample_barcode = substr(tumor_pair_barcode, 6, 12)) 
plot_vaf_hyper$var_filter_revised = factor(plot_vaf_hyper$var_filter_revised, levels = c("private P", "shared P", "shared R", "private R"))

# Hypermutators retained. 31 hypermutation events.
n_distinct(plot_vaf_hyper$tumor_pair_barcode)
pdf(file = "/Users/johnsk/Documents/hypermutants-vaf-distribution-copy neutral.pdf", height = 4, width = 8, bg = "transparent", useDingbats = FALSE)
ggplot(plot_vaf, aes(x=vaf)) + geom_histogram(binwidth = 0.01) + theme_bw() + facet_grid(idh_codel_subtype~var_filter_revised, scales="free") +
  xlab("plot_vaf_hyper frequency (f)") + ylab("Number of mutations") + ggtitle("Hypermutators, copy neutral regions (n=31)")
dev.off()

# Graph for IDHmut-noncodel subtype.
plot_vaf_hyper_noncodel = plot_vaf_hyper %>%  filter(idh_codel_subtype == "IDHmut-noncodel")
pdf(file = "/Users/johnsk/Documents/idh_noncodel_hypermutator_individual-vaf.pdf", height = 4, width = 8, bg = "transparent", useDingbats = FALSE)
ggplot(plot_vaf_hyper_noncodel, aes(x=vaf)) + geom_histogram(binwidth = 0.01) + theme_bw() + xlim(-0.05, 1.05)  + facet_grid(~sample_barcode, scales="free") +
  xlab("Allele frequency (f)") + ylab("Recurrence-only mutations (in R)") + theme(text = element_text(size=8), axis.text.x = element_text(angle=90, hjust=1)) + ggtitle("IDHmut-noncodel hypermutants, copy-neutral regions (n=15)")
dev.off()

# Plot the cancer cell fractions for the hypermutation cases.  
plot_ccf_hyper = vaf_cnv_silver_hyper %>%
  select(tumor_pair_barcode, ccf_a = cellular_prevalence_a, ccf_b = cellular_prevalence_b, fraction, idh_codel_subtype) %>% 
  gather(mutation_subtype, ccf, c(ccf_a, ccf_b), -fraction) %>% 
  mutate(var_filter = paste(mutation_subtype, fraction, sep="_")) %>% 
  filter(var_filter %in% c("ccf_a_P", "ccf_a_S", "ccf_b_S", "ccf_b_R")) %>% 
  mutate(var_filter_revised = recode(var_filter, "ccf_a_P" = "private P", "ccf_a_S" = "shared P", "ccf_b_S" = "shared R", "ccf_b_R" = "private R"),
         sample_barcode = substr(tumor_pair_barcode, 6, 12))
plot_ccf_hyper$var_filter_revised = factor(plot_ccf_hyper$var_filter_revised, levels = c("private P", "shared P", "shared R", "private R"))

# Hypermutators retained. 
n_distinct(plot_ccf_hyper$tumor_pair_barcode)
ggplot(plot_ccf_hyper, aes(x=ccf)) + geom_histogram(binwidth = 0.01) + theme_bw() + facet_grid(idh_codel_subtype~var_filter_revised, scales="free") +
  xlab("Cancer Cell Fraction") + ylab("Number of mutations") + ggtitle("Hypermutators (n=31)")

# Tabulate the number of mutations per sample that are considered clonal:
stacked_clonality_class = vaf_cnv_silver_hyper %>% 
  filter(fraction == "R") %>% 
  group_by(tumor_pair_barcode, clonality_b) %>% 
  summarise (n = n()) %>%
  mutate(freq = n / sum(n))
pdf(file = "/Users/johnsk/Documents/idh_noncodel_hypermutator_clonality.pdf", height = 4, width = 8, bg = "transparent", useDingbats = FALSE)
ggplot(data=stacked_clonality_class, aes(x=substr(tumor_pair_barcode, 6, 12), y=n, fill=clonality_b)) + 
  geom_bar(stat='identity') + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(fill="Hypermutation clonality") +
  xlab("Sample barcode") + ylab("Mutation count R  (Recurrence)") 
dev.off()

# Tabulate the proportion of non-hypermutation events that are synonymous vs. non-synonymous.
stacked_mut_class = vaf_cnv_silver %>% 
  filter(fraction == "R") %>% 
  group_by(tumor_pair_barcode, mut_class) %>% 
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  filter(mut_class == "synonymous")
hist(stacked_mut_class$freq)

stacked_mut_class_hyper = vaf_cnv_silver_hyper %>% 
  filter(fraction == "R") %>% 
  group_by(tumor_pair_barcode, mut_class) %>% 
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  filter(mut_class == "synonymous")
hist(stacked_mut_class_hyper$freq)

################################################
# Hypermutator-specific neutralitytestR models
################################################
# Primary.
codel_primary_nonhyper = vaf_cnv_silver_hyper %>% 
  filter(fraction =="P", idh_codel_subtype == "IDHmut-codel") 
out_primary_1 = neutralitytest(codel_primary_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
out_primary_1 # return statistics
noncodel_primary_nonhyper = vaf_cnv_silver_hyper %>% 
  filter(fraction =="P", idh_codel_subtype == "IDHmut-noncodel") 
out_primary_2 = neutralitytest(noncodel_primary_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
out_primary_2 # return statistics
wt_primary_nonhyper = vaf_cnv_silver_hyper %>% 
  filter(fraction =="P", idh_codel_subtype == "IDHwt") 
out_primary_3 = neutralitytest(wt_primary_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
out_primary_3 # return statistics

# Recurrence.
codel_recur_nonhyper = vaf_cnv_silver_hyper %>% 
  filter(fraction =="R", idh_codel_subtype == "IDHmut-codel") 
out_recur_1 = neutralitytest(codel_recur_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
out_recur_1 # return statistics
noncodel_recur_nonhyper = vaf_cnv_silver_hyper %>% 
  filter(fraction =="R", idh_codel_subtype == "IDHmut-noncodel") 
out_recur_2 = neutralitytest(noncodel_recur_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
out_recur_2 # return statistics
wt_recur_nonhyper = vaf_cnv_silver_hyper %>% 
  filter(fraction =="R", idh_codel_subtype == "IDHwt") 
out_recur_3 = neutralitytest(wt_recur_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
out_recur_3 # return statistics

# SHARED - tumor_a
codel_shared_a_nonhyper = vaf_cnv_silver_hyper %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHmut-codel") 
out_shared_a_1 = neutralitytest(codel_shared_a_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
out_shared_a_1 # return statistics
noncodel_shared_a_nonhyper = vaf_cnv_silver_hyper %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHmut-noncodel") 
out_shared_a_2 = neutralitytest(noncodel_shared_a_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
out_shared_a_2 # return statistics
wt_shared_a_nonhyper = vaf_cnv_silver_hyper %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHwt") 
out_shared_a_3 = neutralitytest(wt_shared_a_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
out_shared_a_3 # return statistics

# SHARED - tumor_b
codel_shared_b_nonhyper = vaf_cnv_silver_hyper %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHmut-codel") 
out_shared_b_1 = neutralitytest(codel_shared_b_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
out_shared_b_1 # return statistics
noncodel_shared_b_nonhyper = vaf_cnv_silver_hyper %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHmut-noncodel") 
out_shared_b_2 = neutralitytest(noncodel_shared_a_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
out_shared_b_2 # return statistics
wt_shared_b_nonhyper = vaf_cnv_silver_hyper %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHwt") 
out_shared_b_3 = neutralitytest(wt_shared_a_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
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
