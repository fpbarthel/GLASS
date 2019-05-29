##################################################
# NeutralitytestR results at the subtype and mutation status-level (private/shared)
# Updated: 2019.05.17
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
library(ggExtra)
library(EnvStats)
library(gridExtra)
library(RColorBrewer)

##################################################
# Establish connection with database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB2")

# Load additional tables.
silver_set = dbReadTable(con,  Id(schema="analysis", table="silver_set"))
gold_set = dbReadTable(con,  Id(schema="analysis", table="gold_set"))
subtypes = dbReadTable(con,  Id(schema="clinical", table="subtypes"))
variant_classifications = dbReadTable(con,  Id(schema="variants", table="variant_classifications"))

# Construct a table that provides subject-level information about clinical variables between two timepoints.
clinical_tumor_pairs_query = read_file("sql/clinical/clinical-tumor-pairs-db2.sql")
clinical_tumor_pairs <- dbGetQuery(con, clinical_tumor_pairs_query)

################################
## Group-wise neutralitytestR ##
################################
# First, analyze the non-hypermutation in recurrent tumor.
# Currently, there is a coverage filter requirement (30X in both tumor_a and tumor_b) PLUS it is required that both tumors have Mutect2 calls for shared variants.
neutrality_input  = read_file("sql/neutrality/neutralitytestr-subtype.sql")
glass_vaf_cnv_unfiltered <- dbGetQuery(con, neutrality_input)

# Examine the VAFs for variants with 30X. Filter based on copy number status and hypermutation event.
vaf_cnv_gold = glass_vaf_cnv_unfiltered %>% 
  filter(tumor_pair_barcode %in% gold_set$tumor_pair_barcode) %>% 
  left_join(clinical_tumor_pairs, by=c("tumor_pair_barcode", "case_barcode", "tumor_barcode_a", "tumor_barcode_b")) %>% 
  left_join(subtypes, by="case_barcode") %>%
  left_join(variant_classifications, by=c("variant_classification"="variant_classification")) %>% 
  mutate(mut_class = ifelse(variant_classification_impact%in%c("HIGH", "MODERATE"), "nonsynonymous", ifelse(variant_classification_impact%in%c("LOW"), "synonymous", NA))) %>% 
  filter(fraction =="P" & cnv_call_a == 0 | fraction =="R" & cnv_call_b == 0 | fraction =="S" & cnv_call_a == 0 & cnv_call_b == 0) %>% 
  filter(hypermutator_status==0) 

# Put the data into a long format, separate by timepoint and fraction.
plot_vaf = vaf_cnv_gold %>% 
  select(tumor_pair_barcode, vaf_a = variant_allele_frequency_a, vaf_b = variant_allele_frequency_b, fraction, mut_class, idh_codel_subtype) %>% 
  gather(mutation_subtype, vaf, c(vaf_a, vaf_b), -fraction) %>% 
  mutate(var_filter = paste(mutation_subtype, fraction, sep="_")) %>% 
  filter(var_filter %in% c("vaf_a_P", "vaf_a_S", "vaf_b_S", "vaf_b_R")) %>% 
  mutate(var_filter_revised = recode(var_filter, "vaf_a_P" = "private P", "vaf_a_S" = "shared P", "vaf_b_S" = "shared R", "vaf_b_R" = "private R"),
         sample_barcode = substr(tumor_pair_barcode, 6, 12)) 
plot_vaf$var_filter_revised = factor(plot_vaf$var_filter_revised, levels = c("private P", "shared P", "shared R", "private R"))

# Hypermutators removed. 181 samples (probably removed samples with lots of CNVs that did not pass). [195 vs. 181 for silver_set vs. gold set].
n_distinct(plot_vaf$tumor_pair_barcode)
# Remove the intronic variants from this visualization.
plot_vaf = plot_vaf %>% filter(!is.na(mut_class))

pdf(file = "/Users/johnsk/Documents/neutralitytestr-nonhyper-goldset-vaf.pdf", height = 6, width = 9, bg = "transparent", useDingbats = FALSE)
ggplot(plot_vaf, aes(x=vaf, fill=mut_class)) + geom_histogram(binwidth = 0.01) + theme_bw() + facet_grid(idh_codel_subtype~var_filter_revised, scales="free") +
  xlab("Allele frequency (f)") + ylab("Number of mutations") + scale_fill_manual(values = c("nonsynonymous" = "#b22222", "synonymous" = "#22B2B2"), name = "Mutation type\n(coding regions)") 
dev.off()


# Nonsynonymous mutations.
nonsyn_plot_vaf = plot_vaf %>% filter(mut_class == "nonsynonymous")
# Sample size may be smaller because of samples with few mutations.
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
plot_ccf = vaf_cnv_gold %>% 
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
plot_shared_ccf = vaf_cnv_gold %>% 
  filter(fraction == "S") %>% 
  mutate(ccf_diff = cellular_prevalence_b-cellular_prevalence_a,
         vaf_diff = variant_allele_frequency_b-variant_allele_frequency_a)
sum(plot_shared_ccf$ccf_diff > 0)/dim(plot_shared_ccf)[1] #  63% have increase in ccf

# Visualize the change in cancer cell fraction for shared variants.
pdf(file = "/Users/johnsk/Documents/nonhypermutants-shared-fraction-delta-ccf.pdf", height = 4, width = 8, bg = "transparent", useDingbats = FALSE)
ggplot(plot_shared_ccf, aes(x=fraction, y=ccf_diff)) + geom_violin() + theme_bw() + ylab("Delta CCF (recurrence-primary)") + xlab("Shared fraction") + facet_grid(~mut_class, scales="free")
dev.off()

################################
## Produce Figure 2B.
################################
# Samples have already been filtered that are not hypermutation events.
codel_primary_nonhyper = vaf_cnv_gold %>% 
  filter(fraction =="P", idh_codel_subtype == "IDHmut-codel", hypermutator_status==0) 
out_primary_1 = neutralitytest(codel_primary_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
out_primary_1 # return statistics
noncodel_primary_nonhyper = vaf_cnv_gold %>% 
  filter(fraction =="P", idh_codel_subtype == "IDHmut-noncodel", hypermutator_status==0) 
out_primary_2 = neutralitytest(noncodel_primary_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
out_primary_2 # return statistics
wt_primary_nonhyper = vaf_cnv_gold %>% 
  filter(fraction =="P", idh_codel_subtype == "IDHwt", hypermutator_status==0) 
out_primary_3 = neutralitytest(wt_primary_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
out_primary_3 # return statistics

# Recurrence.
codel_recur_nonhyper = vaf_cnv_gold %>% 
  filter(fraction =="R", idh_codel_subtype == "IDHmut-codel", hypermutator_status==0) 
out_recur_1 = neutralitytest(codel_recur_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
out_recur_1 # return statistics
noncodel_recur_nonhyper = vaf_cnv_gold %>% 
  filter(fraction =="R", idh_codel_subtype == "IDHmut-noncodel", hypermutator_status==0) 
out_recur_2 = neutralitytest(noncodel_recur_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
out_recur_2 # return statistics
wt_recur_nonhyper = vaf_cnv_gold %>% 
  filter(fraction =="R", idh_codel_subtype == "IDHwt", hypermutator_status==0) 
out_recur_3 = neutralitytest(wt_recur_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
out_recur_3 # return statistics

# SHARED - tumor_a
codel_shared_a_nonhyper = vaf_cnv_gold %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHmut-codel", hypermutator_status==0) 
out_shared_a_1 = neutralitytest(codel_shared_a_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
out_shared_a_1 # return statistics
noncodel_shared_a_nonhyper = vaf_cnv_gold %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHmut-noncodel", hypermutator_status==0) 
out_shared_a_2 = neutralitytest(noncodel_shared_a_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
out_shared_a_2 # return statistics
wt_shared_a_nonhyper = vaf_cnv_gold %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHwt", hypermutator_status==0) 
out_shared_a_3 = neutralitytest(wt_shared_a_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
out_shared_a_3 # return statistics

# SHARED - tumor_b
codel_shared_b_nonhyper = vaf_cnv_gold %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHmut-codel", hypermutator_status==0) 
out_shared_b_1 = neutralitytest(codel_shared_b_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
out_shared_b_1 # return statistics
noncodel_shared_b_nonhyper = vaf_cnv_gold %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHmut-noncodel", hypermutator_status==0) 
out_shared_b_2 = neutralitytest(noncodel_shared_a_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
out_shared_b_2 # return statistics
wt_shared_b_nonhyper = vaf_cnv_gold %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHwt", hypermutator_status==0) 
out_shared_b_3 = neutralitytest(wt_shared_a_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
out_shared_b_3 # return statistics

# Primary
codel_primary = out_primary_1$cumulativefrequency
codel_primary$idh_codel_subtype = "IDHmut-codel"
codel_primary$fraction = "Private P"
noncodel_primary = out_primary_2$cumulativefrequency
noncodel_primary$idh_codel_subtype = "IDHmut-noncodel"
noncodel_primary$fraction = "Private P"
wt_primary = out_primary_3$cumulativefrequency
wt_primary$idh_codel_subtype = "IDHwt"
wt_primary$fraction = "Private P"

# Recurrence
codel_recurrence = out_recur_1$cumulativefrequency
codel_recurrence$idh_codel_subtype = "IDHmut-codel"
codel_recurrence$fraction = "Private R"
noncodel_recurrence = out_recur_2$cumulativefrequency
noncodel_recurrence$idh_codel_subtype = "IDHmut-noncodel"
noncodel_recurrence$fraction = "Private R"
wt_recurrence = out_recur_3$cumulativefrequency
wt_recurrence$idh_codel_subtype = "IDHwt"
wt_recurrence$fraction = "Private R"

# Shared_a
codel_shared_a = out_shared_a_1$cumulativefrequency
codel_shared_a$idh_codel_subtype = "IDHmut-codel"
codel_shared_a$fraction = "Shared P"
noncodel_shared_a = out_shared_a_2$cumulativefrequency
noncodel_shared_a$idh_codel_subtype = "IDHmut-noncodel"
noncodel_shared_a$fraction = "Shared P"
wt_shared_a = out_shared_a_3$cumulativefrequency
wt_shared_a$idh_codel_subtype = "IDHwt"
wt_shared_a$fraction = "Shared P"

# Shared_b
codel_shared_b = out_shared_a_1$cumulativefrequency
codel_shared_b$idh_codel_subtype = "IDHmut-codel"
codel_shared_b$fraction = "Shared R"
noncodel_shared_b = out_shared_b_2$cumulativefrequency
noncodel_shared_b$idh_codel_subtype = "IDHmut-noncodel"
noncodel_shared_b$fraction = "Shared R"
wt_shared_b = out_shared_a_3$cumulativefrequency
wt_shared_b$idh_codel_subtype = "IDHwt"
wt_shared_b$fraction = "Shared R"

# Combine all fractions and subtypes into the same df.
all_samples = rbind(codel_primary, noncodel_primary, wt_primary,
                    codel_recurrence, noncodel_recurrence, wt_recurrence,
                    codel_shared_a, noncodel_shared_a, wt_shared_a,
                    codel_shared_b, noncodel_shared_b, wt_shared_b)
# Specify the order of the fractions.
all_samples$idh_codel_subtype <- as.factor(all_samples$idh_codel_subtype)
all_samples$fraction = factor(all_samples$fraction, levels = c("Private P", "Shared P", "Shared R", "Private R"))

# Create final figure 2B:
p = ggplot2::ggplot(data = all_samples, ggplot2::aes( x=inv_f, y=M_f, col = "1") ) +
  theme_bw() +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x + 0, se=FALSE)   +
  ggplot2::geom_point(ggplot2::aes(colour="2")) +
  ggplot2::scale_colour_manual(values = c("firebrick","black"),
                               labels = c("Best fit line", "Data"),
                               name = "") +
  ggplot2::xlab( "Inverse allelic frequency 1/f" ) +
  ggplot2::ylab( "Cumulative number \nof mutations M(f)" ) +
  ggplot2::scale_x_continuous( trans=scales::identity_trans(), breaks= 1 / c(max(all_samples$f), round((max(all_samples$f)-min(all_samples$f)) / 4, 2), min(all_samples$f)) - 1 / max(all_samples$f),
                               labels= paste("1/", c(max(all_samples$f), round((max(all_samples$f)-min(all_samples$f)) / 4, 2), min(all_samples$f)), sep="")) 

# Define the values to be plotted:
x_y_plot_values = all_samples %>% 
  group_by(idh_codel_subtype, fraction) %>% 
  summarise(max_Mf = max(M_f), 
            min_inv_f = 1/0.25)

# Create text data.frame to be plotted.
txt_df <- as.data.frame(matrix(ncol=6, nrow=12))
colnames(txt_df) <- c("rsq", "pvalue", "idh_codel_subtype", "fraction", "selection", "y")
txt_df$fraction <- as.factor(rep(unique(all_samples$fraction), 3))
txt_df$idh_codel_subtype <- as.factor(c(rep("IDHmut-codel", 4), rep("IDHmut-noncodel", 4), rep("IDHwt", 4)))

# Retrieve the R^2 values.
txt_df$rsq = c(formatC(out_primary_1$rsq$metric, digits = 2), formatC(out_recur_1$rsq$metric, digits = 2), formatC(out_shared_a_1$rsq$metric, digits = 2),
               formatC(out_shared_b_1$rsq$metric, digits = 2), formatC(out_primary_2$rsq$metric, digits = 2), formatC(out_recur_2$rsq$metric, digits = 2), formatC(out_shared_a_2$rsq$metric, digits = 2),
               formatC(out_shared_b_2$rsq$metric, digits = 2), formatC(out_primary_3$rsq$metric, digits = 2), formatC(out_recur_3$rsq$metric, digits = 2), formatC(out_shared_a_3$rsq$metric, digits = 2),
               formatC(out_shared_b_3$rsq$metric, digits = 2))

# Retrieve the p-value values.
txt_df$pvalue = c(formatC(out_primary_1$rsq$pval, format = "e", digits = 2), formatC(out_recur_1$rsq$pval, format = "e", digits = 2), formatC(out_shared_a_1$rsq$pval, format = "e", digits = 2),
                  formatC(out_shared_b_1$rsq$pval, format = "e", digits = 2), formatC(out_primary_2$rsq$pval, format = "e", digits = 2), formatC(out_recur_2$rsq$pval, format = "e", digits = 2), formatC(out_shared_a_2$rsq$pval, format = "e", digits = 2),
                  formatC(out_shared_b_2$rsq$pval, format = "e", digits = 2), formatC(out_primary_3$rsq$pval, format = "e", digits = 2), formatC(out_recur_3$rsq$pval, format = "e", digits = 2), formatC(out_shared_a_3$rsq$pval, format = "e", digits = 2),
                  formatC(out_shared_b_3$rsq$pval, format = "e", digits = 2))
txt_df$selection = ifelse(as.numeric(txt_df$pvalue) < 0.05, "selected", "neutral")

# Use geom_text to add R^2, P-values, and selection classification.
p2 = p + facet_wrap(idh_codel_subtype~fraction, scales = "free") 
txt_df$y = c(layer_scales(p2, 1, 1)$y$range$range[2], layer_scales(p2, 1, 4)$y$range$range[2], layer_scales(p2, 1, 2)$y$range$range[2], layer_scales(p2, 1, 3)$y$range$range[2], 
             layer_scales(p2, 2, 1)$y$range$range[2], layer_scales(p2, 2, 4)$y$range$range[2], layer_scales(p2, 2, 2)$y$range$range[2], layer_scales(p2, 2, 3)$y$range$range[2],
             layer_scales(p2, 3, 1)$y$range$range[2], layer_scales(p2, 3, 4)$y$range$range[2], layer_scales(p2, 3, 2)$y$range$range[2], layer_scales(p2, 3, 3)$y$range$range[2])
txt_df$pval = paste0("p = ", txt_df$pvalue, sep="")
txt_df$rsquared = paste0("R^2 = ", txt_df$rsq, sep="")

# Create final output for Figure 2b.
pdf(file = "/Users/johnsk/Documents/f2b-kcj.pdf", height = 6, width = 9, bg = "transparent", useDingbats = FALSE)
p + geom_text(data=txt_df, aes(x = 0.5, y= 1.1*y, label=rsquared, family="serif"), hjust=0) +
  geom_text(data=txt_df, aes(x = 0.5, y = y, label=pval, family="serif"), hjust=0) +
  geom_text(data=txt_df, aes(x = 0.5, y = 0.9*y, label=selection, family="serif"), hjust=0) +
  facet_wrap(idh_codel_subtype~fraction, scales = "free") 
dev.off()


###########################
##### Hypermutators     ###
###########################
vaf_cnv_gold_hyper = glass_vaf_cnv_unfiltered %>% 
  filter(tumor_pair_barcode %in% gold_set$tumor_pair_barcode) %>% 
  left_join(clinical_tumor_pairs, by=c("tumor_pair_barcode", "case_barcode", "tumor_barcode_a", "tumor_barcode_b")) %>% 
  left_join(subtypes, by="case_barcode") %>%
  left_join(variant_classifications, by=c("variant_classification"="variant_classification")) %>% 
  mutate(mut_class = ifelse(variant_classification_impact%in%c("HIGH", "MODERATE"), "nonsynonymous", ifelse(variant_classification_impact%in%c("LOW"), "synonymous", NA))) %>% 
  filter(fraction =="P" & cnv_call_a == 0 | fraction =="R" & cnv_call_b == 0 | fraction =="S" & cnv_call_a == 0 & cnv_call_b == 0) %>% 
  filter(hypermutator_status==1) 

# Samples have already been filtered that are not hypermutation events.
codel_primary_hyper = vaf_cnv_gold_hyper %>% 
  filter(fraction =="P", idh_codel_subtype == "IDHmut-codel") 
out_primary_1 = neutralitytest(codel_primary_hyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
out_primary_1 # return statistics
noncodel_primary_hyper = vaf_cnv_gold_hyper %>% 
  filter(fraction =="P", idh_codel_subtype == "IDHmut-noncodel") 
out_primary_2 = neutralitytest(noncodel_primary_hyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
out_primary_2 # return statistics
wt_primary_hyper = vaf_cnv_gold_hyper %>% 
  filter(fraction =="P", idh_codel_subtype == "IDHwt") 
out_primary_3 = neutralitytest(wt_primary_hyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
out_primary_3 # return statistics

# Recurrence.
codel_recur_hyper = vaf_cnv_gold_hyper %>% 
  filter(fraction =="R", idh_codel_subtype == "IDHmut-codel") 
out_recur_1 = neutralitytest(codel_recur_hyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
out_recur_1 # return statistics
noncodel_recur_hyper = vaf_cnv_gold_hyper %>% 
  filter(fraction =="R", idh_codel_subtype == "IDHmut-noncodel") 
out_recur_2 = neutralitytest(noncodel_recur_hyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
out_recur_2 # return statistics
wt_recur_hyper = vaf_cnv_gold_hyper %>% 
  filter(fraction =="R", idh_codel_subtype == "IDHwt") 
out_recur_3 = neutralitytest(wt_recur_hyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
out_recur_3 # return statistics

# SHARED - tumor_a
codel_shared_a_hyper = vaf_cnv_gold_hyper %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHmut-codel") 
out_shared_a_1 = neutralitytest(codel_shared_a_hyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
out_shared_a_1 # return statistics
noncodel_shared_a_hyper = vaf_cnv_gold_hyper %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHmut-noncodel") 
out_shared_a_2 = neutralitytest(noncodel_shared_a_hyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
out_shared_a_2 # return statistics
wt_shared_a_hyper = vaf_cnv_gold_hyper %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHwt") 
out_shared_a_3 = neutralitytest(wt_shared_a_hyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
out_shared_a_3 # return statistics

# SHARED - tumor_b
codel_shared_b_hyper = vaf_cnv_gold_hyper %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHmut-codel") 
out_shared_b_1 = neutralitytest(codel_shared_b_hyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
out_shared_b_1 # return statistics
noncodel_shared_b_hyper = vaf_cnv_gold_hyper %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHmut-noncodel") 
out_shared_b_2 = neutralitytest(noncodel_shared_a_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
out_shared_b_2 # return statistics
wt_shared_b_hyper = vaf_cnv_gold_hyper %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHwt") 
out_shared_b_3 = neutralitytest(wt_shared_b_hyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
out_shared_b_3 # return statistics

# Primary
codel_primary = out_primary_1$cumulativefrequency
codel_primary$idh_codel_subtype = "IDHmut-codel"
codel_primary$fraction = "Private P"
noncodel_primary = out_primary_2$cumulativefrequency
noncodel_primary$idh_codel_subtype = "IDHmut-noncodel"
noncodel_primary$fraction = "Private P"
wt_primary = out_primary_3$cumulativefrequency
wt_primary$idh_codel_subtype = "IDHwt"
wt_primary$fraction = "Private P"

# Recurrence
codel_recurrence = out_recur_1$cumulativefrequency
codel_recurrence$idh_codel_subtype = "IDHmut-codel"
codel_recurrence$fraction = "Private R"
noncodel_recurrence = out_recur_2$cumulativefrequency
noncodel_recurrence$idh_codel_subtype = "IDHmut-noncodel"
noncodel_recurrence$fraction = "Private R"
wt_recurrence = out_recur_3$cumulativefrequency
wt_recurrence$idh_codel_subtype = "IDHwt"
wt_recurrence$fraction = "Private R"

# Shared_a
codel_shared_a = out_shared_a_1$cumulativefrequency
codel_shared_a$idh_codel_subtype = "IDHmut-codel"
codel_shared_a$fraction = "Shared P"
noncodel_shared_a = out_shared_a_2$cumulativefrequency
noncodel_shared_a$idh_codel_subtype = "IDHmut-noncodel"
noncodel_shared_a$fraction = "Shared P"
wt_shared_a = out_shared_a_3$cumulativefrequency
wt_shared_a$idh_codel_subtype = "IDHwt"
wt_shared_a$fraction = "Shared P"

# Shared_b
codel_shared_b = out_shared_a_1$cumulativefrequency
codel_shared_b$idh_codel_subtype = "IDHmut-codel"
codel_shared_b$fraction = "Shared R"
noncodel_shared_b = out_shared_b_2$cumulativefrequency
noncodel_shared_b$idh_codel_subtype = "IDHmut-noncodel"
noncodel_shared_b$fraction = "Shared R"
wt_shared_b = out_shared_a_3$cumulativefrequency
wt_shared_b$idh_codel_subtype = "IDHwt"
wt_shared_b$fraction = "Shared R"

# Combine all fractions and subtypes into the same df.
all_samples = rbind(codel_primary, noncodel_primary, wt_primary,
                    codel_recurrence, noncodel_recurrence, wt_recurrence,
                    codel_shared_a, noncodel_shared_a, wt_shared_a,
                    codel_shared_b, noncodel_shared_b, wt_shared_b)
# Specify the order of the fractions.
all_samples$idh_codel_subtype <- as.factor(all_samples$idh_codel_subtype)
all_samples$fraction = factor(all_samples$fraction, levels = c("Private P", "Shared P", "Shared R", "Private R"))

# Create final figure 2B:
p = ggplot2::ggplot(data = all_samples, ggplot2::aes( x=inv_f, y=M_f, col = "1") ) +
  theme_bw() +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x + 0, se=FALSE)   +
  ggplot2::geom_point(ggplot2::aes(colour="2")) +
  ggplot2::scale_colour_manual(values = c("firebrick","black"),
                               labels = c("Best fit line", "Data"),
                               name = "") +
  ggplot2::xlab( "Inverse allelic frequency 1/f" ) +
  ggplot2::ylab( "Cumulative number \nof mutations M(f)" ) +
  ggplot2::scale_x_continuous( trans=scales::identity_trans(), breaks= 1 / c(max(all_samples$f), round((max(all_samples$f)-min(all_samples$f)) / 4, 2), min(all_samples$f)) - 1 / max(all_samples$f),
                               labels= paste("1/", c(max(all_samples$f), round((max(all_samples$f)-min(all_samples$f)) / 4, 2), min(all_samples$f)), sep="")) 

# Define the values to be plotted:
x_y_plot_values = all_samples %>% 
  group_by(idh_codel_subtype, fraction) %>% 
  summarise(max_Mf = max(M_f), 
            min_inv_f = 1/0.25)

# Create text data.frame to be plotted.
txt_df <- as.data.frame(matrix(ncol=6, nrow=12))
colnames(txt_df) <- c("rsq", "pvalue", "idh_codel_subtype", "fraction", "selection", "y")
txt_df$fraction <- as.factor(rep(unique(all_samples$fraction), 3))
txt_df$idh_codel_subtype <- as.factor(c(rep("IDHmut-codel", 4), rep("IDHmut-noncodel", 4), rep("IDHwt", 4)))

# Retrieve the R^2 values.
txt_df$rsq = c(formatC(out_primary_1$rsq$metric, digits = 2), formatC(out_recur_1$rsq$metric, digits = 2), formatC(out_shared_a_1$rsq$metric, digits = 2),
               formatC(out_shared_b_1$rsq$metric, digits = 2), formatC(out_primary_2$rsq$metric, digits = 2), formatC(out_recur_2$rsq$metric, digits = 2), formatC(out_shared_a_2$rsq$metric, digits = 2),
               formatC(out_shared_b_2$rsq$metric, digits = 2), formatC(out_primary_3$rsq$metric, digits = 2), formatC(out_recur_3$rsq$metric, digits = 2), formatC(out_shared_a_3$rsq$metric, digits = 2),
               formatC(out_shared_b_3$rsq$metric, digits = 2))

# Retrieve the p-value values.
txt_df$pvalue = c(formatC(out_primary_1$rsq$pval, format = "e", digits = 2), formatC(out_recur_1$rsq$pval, format = "e", digits = 2), formatC(out_shared_a_1$rsq$pval, format = "e", digits = 2),
                  formatC(out_shared_b_1$rsq$pval, format = "e", digits = 2), formatC(out_primary_2$rsq$pval, format = "e", digits = 2), formatC(out_recur_2$rsq$pval, format = "e", digits = 2), formatC(out_shared_a_2$rsq$pval, format = "e", digits = 2),
                  formatC(out_shared_b_2$rsq$pval, format = "e", digits = 2), formatC(out_primary_3$rsq$pval, format = "e", digits = 2), formatC(out_recur_3$rsq$pval, format = "e", digits = 2), formatC(out_shared_a_3$rsq$pval, format = "e", digits = 2),
                  formatC(out_shared_b_3$rsq$pval, format = "e", digits = 2))
txt_df$selection = ifelse(as.numeric(txt_df$pvalue) < 0.05, "selected", "neutral")

# Use geom_text to add R^2, P-values, and selection classification.
p2 = p + facet_wrap(idh_codel_subtype~fraction, scales = "free") 
txt_df$y = c(layer_scales(p2, 1, 1)$y$range$range[2], layer_scales(p2, 1, 4)$y$range$range[2], layer_scales(p2, 1, 2)$y$range$range[2], layer_scales(p2, 1, 3)$y$range$range[2], 
             layer_scales(p2, 2, 1)$y$range$range[2], layer_scales(p2, 2, 4)$y$range$range[2], layer_scales(p2, 2, 2)$y$range$range[2], layer_scales(p2, 2, 3)$y$range$range[2],
             layer_scales(p2, 3, 1)$y$range$range[2], layer_scales(p2, 3, 4)$y$range$range[2], layer_scales(p2, 3, 2)$y$range$range[2], layer_scales(p2, 3, 3)$y$range$range[2])
txt_df$pval = paste0("p = ", txt_df$pvalue, sep="")
txt_df$rsquared = paste0("R^2 = ", txt_df$rsq, sep="")

# Create final output for Figure 2b.
pdf(file = "/Users/johnsk/Documents/neutralitytestr-hypermutators-goldset.pdf", height = 6, width = 9, bg = "transparent", useDingbats = FALSE)
p + geom_text(data=txt_df, aes(x = 0.5, y= 1.1*y, label=rsquared, family="serif"), hjust=0) +
  geom_text(data=txt_df, aes(x = 0.5, y = y, label=pval, family="serif"), hjust=0) +
  geom_text(data=txt_df, aes(x = 0.5, y = 0.9*y, label=selection, family="serif"), hjust=0) +
  facet_wrap(idh_codel_subtype~fraction, scales = "free") 
dev.off()


plot_vaf_hyper = vaf_cnv_gold_hyper %>% 
  select(tumor_pair_barcode, vaf_a = variant_allele_frequency_a, vaf_b = variant_allele_frequency_b, fraction, mut_class, idh_codel_subtype) %>% 
  gather(mutation_subtype, vaf, c(vaf_a, vaf_b), -fraction) %>% 
  mutate(var_filter = paste(mutation_subtype, fraction, sep="_")) %>% 
  filter(var_filter %in% c("vaf_a_P", "vaf_a_S", "vaf_b_S", "vaf_b_R")) %>% 
  mutate(var_filter_revised = recode(var_filter, "vaf_a_P" = "private P", "vaf_a_S" = "shared P", "vaf_b_S" = "shared R", "vaf_b_R" = "private R"),
         sample_barcode = substr(tumor_pair_barcode, 6, 12)) 
plot_vaf_hyper$var_filter_revised = factor(plot_vaf_hyper$var_filter_revised, levels = c("private P", "shared P", "shared R", "private R"))

# Remove mutations in non-coding regions.
plot_vaf_hyper = plot_vaf_hyper %>% filter(!is.na(mut_class))
pdf(file = "/Users/johnsk/Documents/neutralitytestr-hyper-goldset-vaf.pdf", height = 6, width = 9, bg = "transparent", useDingbats = FALSE)
ggplot(plot_vaf_hyper, aes(x=vaf, fill=mut_class)) + geom_histogram(binwidth = 0.01)+ theme_bw() + facet_grid(idh_codel_subtype~var_filter_revised, scales="free") +
  xlab("Allele frequency (f)") + ylab("Number of mutations") + scale_fill_manual(values = c("nonsynonymous" = "#b22222", "synonymous" = "#22B2B2"), name = "Mutation type\n(coding regions)") 
dev.off()

# Create density plot.
pdf(file = "/Users/johnsk/Documents/neutralitytestr-hyper-goldset-density.pdf", height = 6, width = 9, bg = "transparent", useDingbats = FALSE)
ggplot(plot_vaf_hyper, aes(x=vaf, color=mut_class)) +geom_density() + theme_bw() + facet_grid(idh_codel_subtype~var_filter_revised, scales="free") +
  xlab("Allele frequency (f)") + ylab("Number of mutations") + scale_color_manual(values = c("nonsynonymous" = "#b22222", "synonymous" = "#22B2B2"), name = "Mutation type\n(coding regions)") 
dev.off()


##################################
## Hypermutator-specific analyses
###################################
# Select for hypermutation events.
vaf_cnv_gold_hyper = glass_vaf_cnv_unfiltered %>% 
  filter(tumor_pair_barcode %in% gold_set$tumor_pair_barcode) %>% 
  left_join(clinical_tumor_pairs, by=c("tumor_pair_barcode", "case_barcode", "tumor_barcode_a", "tumor_barcode_b")) %>% 
  left_join(subtypes, by="case_barcode") %>%
  left_join(variant_classifications, by=c("variant_classification"="variant_classification")) %>% 
  mutate(mut_class = ifelse(variant_classification_impact%in%c("HIGH", "MODERATE"), "nonsynonymous", ifelse(variant_classification_impact%in%c("LOW"), "synonymous", NA))) %>% 
  filter(fraction =="P" & cnv_call_a == 0 | fraction =="R" & cnv_call_b == 0 | fraction =="S" & cnv_call_a == 0 & cnv_call_b == 0) %>% 
  filter(hypermutator_status==1) 

# Create a vaf plot by timepoint and subtype.
plot_vaf_hyper = vaf_cnv_gold_hyper %>% 
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
plot_ccf_hyper = vaf_cnv_gold_hyper %>%
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
stacked_clonality_class = vaf_cnv_gold_hyper %>% 
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
stacked_mut_class = vaf_cnv_gold %>% 
  filter(fraction == "R") %>% 
  group_by(tumor_pair_barcode, mut_class) %>% 
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  filter(mut_class == "synonymous")
hist(stacked_mut_class$freq)

stacked_mut_class_hyper = vaf_cnv_gold_hyper %>% 
  filter(fraction == "R") %>% 
  group_by(tumor_pair_barcode, mut_class) %>% 
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  filter(mut_class == "synonymous")
hist(stacked_mut_class_hyper$freq)

###############################################
# NeutralitytestR for Figure 2B with histograms (associated R^2, P-values)
###############################################
# Color with the nonsynonymous and synonymous colors for R^2 and pvalue.

# Remove the intronic variants from this visualization.
# **nonsynonymous**
vaf_cnv_gold_nonsyn = vaf_cnv_gold %>% filter(mut_class == "nonsynonymous")

# Primary
codel_primary_nonhyper = vaf_cnv_gold_nonsyn %>% 
  filter(fraction =="P", idh_codel_subtype == "IDHmut-codel", hypermutator_status==0) 
out_primary_1 = neutralitytest(codel_primary_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
noncodel_primary_nonhyper = vaf_cnv_gold_nonsyn %>% 
  filter(fraction =="P", idh_codel_subtype == "IDHmut-noncodel", hypermutator_status==0) 
out_primary_2 = neutralitytest(noncodel_primary_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
wt_primary_nonhyper = vaf_cnv_gold_nonsyn %>% 
  filter(fraction =="P", idh_codel_subtype == "IDHwt", hypermutator_status==0) 
out_primary_3 = neutralitytest(wt_primary_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)

# Recurrence.
codel_recur_nonhyper = vaf_cnv_gold_nonsyn %>% 
  filter(fraction =="R", idh_codel_subtype == "IDHmut-codel", hypermutator_status==0) 
out_recur_1 = neutralitytest(codel_recur_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
noncodel_recur_nonhyper = vaf_cnv_gold_nonsyn %>% 
  filter(fraction =="R", idh_codel_subtype == "IDHmut-noncodel", hypermutator_status==0) 
out_recur_2 = neutralitytest(noncodel_recur_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
wt_recur_nonhyper = vaf_cnv_gold_nonsyn %>% 
  filter(fraction =="R", idh_codel_subtype == "IDHwt", hypermutator_status==0) 
out_recur_3 = neutralitytest(wt_recur_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)

# SHARED - tumor_a
codel_shared_a_nonhyper = vaf_cnv_gold_nonsyn %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHmut-codel", hypermutator_status==0) 
out_shared_a_1 = neutralitytest(codel_shared_a_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
noncodel_shared_a_nonhyper = vaf_cnv_gold_nonsyn %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHmut-noncodel", hypermutator_status==0) 
out_shared_a_2 = neutralitytest(noncodel_shared_a_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
wt_shared_a_nonhyper = vaf_cnv_gold_nonsyn %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHwt", hypermutator_status==0) 
out_shared_a_3 = neutralitytest(wt_shared_a_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)

# SHARED - tumor_b
codel_shared_b_nonhyper = vaf_cnv_gold_nonsyn %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHmut-codel", hypermutator_status==0) 
out_shared_b_1 = neutralitytest(codel_shared_b_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
noncodel_shared_b_nonhyper = vaf_cnv_gold_nonsyn %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHmut-noncodel", hypermutator_status==0) 
out_shared_b_2 = neutralitytest(noncodel_shared_a_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
wt_shared_b_nonhyper = vaf_cnv_gold_nonsyn %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHwt", hypermutator_status==0) 
out_shared_b_3 = neutralitytest(wt_shared_a_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)

# ** synonymous **
vaf_cnv_gold_syn = vaf_cnv_gold %>% filter(mut_class == "synonymous")

# Primary
codel_primary_nonhyper = vaf_cnv_gold_syn %>% 
  filter(fraction =="P", idh_codel_subtype == "IDHmut-codel", hypermutator_status==0) 
out_primary_1s = neutralitytest(codel_primary_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
noncodel_primary_nonhyper = vaf_cnv_gold_syn %>% 
  filter(fraction =="P", idh_codel_subtype == "IDHmut-noncodel", hypermutator_status==0) 
out_primary_2s = neutralitytest(noncodel_primary_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
wt_primary_nonhyper = vaf_cnv_gold_syn %>% 
  filter(fraction =="P", idh_codel_subtype == "IDHwt", hypermutator_status==0) 
out_primary_3s = neutralitytest(wt_primary_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)

# Recurrence.
codel_recur_nonhyper = vaf_cnv_gold_syn %>% 
  filter(fraction =="R", idh_codel_subtype == "IDHmut-codel", hypermutator_status==0) 
out_recur_1s = neutralitytest(codel_recur_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
noncodel_recur_nonhyper = vaf_cnv_gold_syn %>% 
  filter(fraction =="R", idh_codel_subtype == "IDHmut-noncodel", hypermutator_status==0) 
out_recur_2s = neutralitytest(noncodel_recur_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
wt_recur_nonhyper = vaf_cnv_gold_syn %>% 
  filter(fraction =="R", idh_codel_subtype == "IDHwt", hypermutator_status==0) 
out_recur_3s = neutralitytest(wt_recur_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)

# SHARED - tumor_a
codel_shared_a_nonhyper = vaf_cnv_gold_syn %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHmut-codel", hypermutator_status==0) 
out_shared_a_1s = neutralitytest(codel_shared_a_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
noncodel_shared_a_nonhyper = vaf_cnv_gold_syn %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHmut-noncodel", hypermutator_status==0) 
out_shared_a_2s = neutralitytest(noncodel_shared_a_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)
wt_shared_a_nonhyper = vaf_cnv_gold_syn %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHwt", hypermutator_status==0) 
out_shared_a_3s = neutralitytest(wt_shared_a_nonhyper$variant_allele_frequency_a, fmin = 0.1, fmax = 0.25)

# SHARED - tumor_b
codel_shared_b_nonhyper = vaf_cnv_gold_syn %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHmut-codel", hypermutator_status==0) 
out_shared_b_1s = neutralitytest(codel_shared_b_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
noncodel_shared_b_nonhyper = vaf_cnv_gold_syn %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHmut-noncodel", hypermutator_status==0) 
out_shared_b_2s = neutralitytest(noncodel_shared_a_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)
wt_shared_b_nonhyper = vaf_cnv_gold_syn %>% 
  filter(fraction =="S", idh_codel_subtype == "IDHwt", hypermutator_status==0) 
out_shared_b_3s = neutralitytest(wt_shared_a_nonhyper$variant_allele_frequency_b, fmin = 0.1, fmax = 0.25)


# Plotting.
plot_vaf = plot_vaf %>% 
  mutate(var_filter_revised = recode(var_filter_revised, "private P"= "private I", "shared P"= "shared I", "shared R"= "shared R", "private R" = "private R"))
levels(plot_vaf$var_filter_revised)
p = ggplot(plot_vaf, aes(x=vaf, fill=mut_class)) + geom_histogram(binwidth = 0.01) + theme_bw() +
  coord_cartesian(xlim = c(0,1)) +
  xlab("Allele frequency (f)") + ylab("Number of mutations") + scale_fill_manual(values = c("nonsynonymous" = "#b22222", "synonymous" = "#22B2B2"), name = "Mutation type\n(coding regions)") 


# Create text data.frame to be plotted.
txt_df <- as.data.frame(matrix(ncol=7, nrow=24))
colnames(txt_df) <- c("rsq", "pvalue", "idh_codel_subtype", "var_filter_revised", "selection", "mut_class", "y")
txt_df$var_filter_revised <- as.factor(rep(c("private I", "private R", "shared I", "shared R"), 6))
txt_df$idh_codel_subtype <- as.factor(c(rep("IDHmut-codel", 4), rep("IDHmut-noncodel", 4), rep("IDHwt", 4), 
                                        rep("IDHmut-codel", 4), rep("IDHmut-noncodel", 4), rep("IDHwt", 4)))
txt_df$mut_class <- c(rep("nonsynonymous", 12), rep("synonymous", 12))
txt_df$var_filter_revised = factor(txt_df$var_filter_revised, levels = c("private I", "shared I", "shared R", "private R"))

# Retrieve the R^2 values.
txt_df$rsq = c(formatC(out_primary_1$rsq$metric, digits = 2), formatC(out_recur_1$rsq$metric, digits = 2), formatC(out_shared_a_1$rsq$metric, digits = 2),
               formatC(out_shared_b_1$rsq$metric, digits = 2), formatC(out_primary_2$rsq$metric, digits = 2), formatC(out_recur_2$rsq$metric, digits = 2), formatC(out_shared_a_2$rsq$metric, digits = 2),
               formatC(out_shared_b_2$rsq$metric, digits = 2), formatC(out_primary_3$rsq$metric, digits = 2), formatC(out_recur_3$rsq$metric, digits = 2), formatC(out_shared_a_3$rsq$metric, digits = 2),
               formatC(out_shared_b_3$rsq$metric, digits = 2), formatC(out_primary_1s$rsq$metric, digits = 2), formatC(out_recur_1s$rsq$metric, digits = 2), formatC(out_shared_a_1s$rsq$metric, digits = 2),
               formatC(out_shared_b_1s$rsq$metric, digits = 2), formatC(out_primary_2s$rsq$metric, digits = 2), formatC(out_recur_2s$rsq$metric, digits = 2), formatC(out_shared_a_2s$rsq$metric, digits = 2),
               formatC(out_shared_b_2s$rsq$metric, digits = 2), formatC(out_primary_3s$rsq$metric, digits = 2), formatC(out_recur_3s$rsq$metric, digits = 2), formatC(out_shared_a_3s$rsq$metric, digits = 2),
               formatC(out_shared_b_3s$rsq$metric, digits = 2))

# Retrieve the p-value values.
txt_df$pvalue = c(formatC(out_primary_1$rsq$pval, format = "e", digits = 2), formatC(out_recur_1$rsq$pval, format = "e", digits = 2), formatC(out_shared_a_1$rsq$pval, format = "e", digits = 2),
                  formatC(out_shared_b_1$rsq$pval, format = "e", digits = 2), formatC(out_primary_2$rsq$pval, format = "e", digits = 2), formatC(out_recur_2$rsq$pval, format = "e", digits = 2), formatC(out_shared_a_2$rsq$pval, format = "e", digits = 2),
                  formatC(out_shared_b_2$rsq$pval, format = "e", digits = 2), formatC(out_primary_3$rsq$pval, format = "e", digits = 2), formatC(out_recur_3$rsq$pval, format = "e", digits = 2), formatC(out_shared_a_3$rsq$pval, format = "e", digits = 2),
                  formatC(out_shared_b_3$rsq$pval, format = "e", digits = 2), formatC(out_primary_1s$rsq$pval, format = "e", digits = 2), formatC(out_recur_1s$rsq$pval, format = "e", digits = 2), formatC(out_shared_a_1s$rsq$pval, format = "e", digits = 2),
                  formatC(out_shared_b_1s$rsq$pval, format = "e", digits = 2), formatC(out_primary_2s$rsq$pval, format = "e", digits = 2), formatC(out_recur_2s$rsq$pval, format = "e", digits = 2), formatC(out_shared_a_2s$rsq$pval, format = "e", digits = 2),
                  formatC(out_shared_b_2s$rsq$pval, format = "e", digits = 2), formatC(out_primary_3s$rsq$pval, format = "e", digits = 2), formatC(out_recur_3s$rsq$pval, format = "e", digits = 2), formatC(out_shared_a_3s$rsq$pval, format = "e", digits = 2),
                  formatC(out_shared_b_3s$rsq$pval, format = "e", digits = 2))
txt_df$selection = ifelse(as.numeric(txt_df$pvalue) < 0.05, "selected", "neutral")

# Use geom_text to add R^2, P-values, and selection classification.
p2 = p + facet_wrap(idh_codel_subtype~var_filter_revised, scales="free")
txt_df$y = rep(c(layer_scales(p2, 1, 1)$y$range$range[2], layer_scales(p2, 1, 4)$y$range$range[2], layer_scales(p2, 1, 2)$y$range$range[2], layer_scales(p2, 1, 3)$y$range$range[2], 
             layer_scales(p2, 2, 1)$y$range$range[2], layer_scales(p2, 2, 4)$y$range$range[2], layer_scales(p2, 2, 2)$y$range$range[2], layer_scales(p2, 2, 3)$y$range$range[2],
             layer_scales(p2, 3, 1)$y$range$range[2], layer_scales(p2, 3, 4)$y$range$range[2], layer_scales(p2, 3, 2)$y$range$range[2], layer_scales(p2, 3, 3)$y$range$range[2]), 2)
txt_df$pval = paste0("p = ", txt_df$pvalue, sep="")
txt_df$rsquared = paste0("R^2 = ", txt_df$rsq, sep="")
# Separate the text annotation into the mutation classes.
txt_df_ns <- txt_df[1:12, ]
txt_df_s <- txt_df[13:24, ]

pdf(file = "/Users/johnsk/Documents/a-neutralitytestr-nonhyper-goldset-vaf.pdf", height = 9, width = 12, bg = "transparent", useDingbats = FALSE)
p + geom_text(data=txt_df_s, aes(x = 0.5, y= 0.9*y, label=rsquared, family="serif", color="#b22222"), hjust=0) +
  geom_text(data=txt_df_s, aes(x = 0.5, y = 0.8*y, label=pval, family="serif", color="#b22222"), hjust=0) +
  geom_text(data=txt_df_ns, aes(x = 0.5, y= 1.1*y, label=rsquared, family="serif", color="#22B2B2"), hjust=0) +
  geom_text(data=txt_df_ns, aes(x = 0.5, y = y, label=pval, family="serif", color="#22B2B2"), hjust=0) + 
  guides(color = FALSE) +
  facet_wrap(idh_codel_subtype~var_filter_revised, scales = "free") 
dev.off()
