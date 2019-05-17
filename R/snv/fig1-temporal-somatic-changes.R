##################################################
# Determine association between time and somatic alteration acquisition.
# Updated: 2019.04.23
# Author: Kevin J.
##################################################

## OBJECTIVE: Address Reviewer #1's comment:
# Is there a correlation between time interval to recurrence and number of genetic alterations accrued in that time period? 
# We can test differences in mutational frequency/aneuploidy with surgical interval. 

# Working directory for this analysis in the GLASS-analysis project. 
mybasedir = "/Volumes/verhaak-lab/GLASS-analysis/"
setwd(mybasedir)

##################################################
# Necessary packages:
library(tidyverse)
library(DBI)
library(ggiraphExtra)

##################################################
# Establish connection with database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB2")

# Load in the necessary information regarding mutation frequency and all clinical information.
mutation_freq = dbGetQuery(con, "SELECT * FROM analysis.mut_freq")
silver_set = dbReadTable(con,  Id(schema="analysis", table="silver_set"))
gold_set = dbReadTable(con,  Id(schema="analysis", table="gold_set"))
tumor_pairs = dbReadTable(con,  Id(schema="analysis", table="tumor_pairs"))
gatk_aneuploidy = dbReadTable(con,  Id(schema="analysis", table="gatk_aneuploidy"))
subtypes = dbReadTable(con,  Id(schema="clinical", table="subtypes"))
cases = dbReadTable(con,  Id(schema="clinical", table="cases"))
surgeries = dbReadTable(con,  Id(schema="clinical", table="surgeries"))

# Clinical tumor pairs table including treatment information.
clinical_tumor_pairs_query = read_file("/Users/johnsk/Documents/Life-History/glass-analyses/sql/clinical-tumor-pairs-db2.sql")
clinical_tumor_pairs <- dbGetQuery(con, clinical_tumor_pairs_query)

###################################################################
# Time and CNV changes
###################################################################
# The proportion of different segments for all tumor pairs, also subdivided by amplification and deletion
gatk_seg_diff_prop = dbReadTable(con,  Id(schema="analysis", table="gatk_seg_diff_prop"))

# Combine with the `gold set` to know which samples have high-quality CNV data.
gatk_seg_diff_gold = gatk_seg_diff_prop %>% 
  filter(tumor_pair_barcode %in% gold_set$tumor_pair_barcode) %>% 
  left_join(clinical_tumor_pairs, by=c("tumor_pair_barcode", "case_barcode", "tumor_barcode_a", "tumor_barcode_b")) %>% 
  left_join(subtypes, by="case_barcode") %>% 
  filter(!is.na(surgical_interval))

# Define aneuploidy seg differences by IDH subtype.
gatk_seg_IDH_codel = gatk_seg_diff_gold %>% 
  filter(idh_codel_subtype == "IDHmut-codel", !is.na(surgical_interval))
gatk_seg_IDH_noncodel = gatk_seg_diff_gold %>% 
  filter(idh_codel_subtype == "IDHmut-noncodel", !is.na(surgical_interval))
gatk_seg_IDHwt = gatk_seg_diff_gold %>% 
  filter(idh_codel_subtype == "IDHwt", !is.na(surgical_interval))

# Histogram of surgical interval between samples with two separate timepoints.
ggplot(gatk_seg_diff_gold, aes(surgical_interval)) + geom_histogram() + theme_bw() + xlab("Surgical interval (months)") + ylab("Counts") + facet_wrap(~idh_codel_subtype) 
ggplot(gatk_seg_diff_gold, aes(prop_change)) + geom_histogram() + theme_bw() + xlab("Proportion of differential segments") + ylab("Counts") + facet_wrap(~idh_codel_subtype) 

# Transform the surgical interval so that the data is more normally distributed.
ggplot(gatk_seg_diff_gold, aes(log10(surgical_interval))) + geom_histogram() + theme_bw() + xlab("log10(Surgical interval (months))") + ylab("Counts") + facet_wrap(~idh_codel_subtype) 
# May not be necessary to log transform the proportion of copy number changes.
ggplot(gatk_seg_diff_gold, aes(log10(prop_change))) + geom_histogram() + theme_bw() + xlab("log10(Proportion of differential segments)") + ylab("Counts") + facet_wrap(~idh_codel_subtype) 

# TOTAL proportion change as it relates to time between surgeries.
pdf(file = "/Users/johnsk/Documents/Life-History/glass-analyses/figures/temporal-somatic-alterations/gold-time-total-prop-gatk-cnv-seg-diff.pdf", height = 4, width = 8, bg = "transparent", useDingbats = FALSE)
ggplot(gatk_seg_diff_gold, aes(x=log10(surgical_interval), y=prop_change, color= idh_codel_subtype)) +
  geom_point() + theme_bw() + xlab("log10(Surgical interval (months))") + ylab("Proportion of CNV \n segment differences") + labs("Tumor") + geom_smooth(method='lm', formula =y~x) +
  facet_wrap(~idh_codel_subtype) + labs(color = "Glioma subtype") + ggtitle("GLASS gold set (n = 222)")
dev.off()

# AMPLIFICATION proportion change.
ggplot(gatk_seg_diff_gold, aes(x=log10(surgical_interval), y=prop_amp, color= idh_codel_subtype)) +
  geom_point() + theme_bw() + xlab("") + ylab("") + labs("Tumor") + geom_smooth(method='lm', formula =y~x) +
  facet_wrap(~idh_codel_subtype, scale = "free")  
# DELETION proportion change.
ggplot(gatk_seg_diff_gold, aes(x=surgical_interval, y=prop_del, color= idh_codel_subtype)) +
  geom_point() + theme_bw() + xlab("") + ylab("") + geom_smooth(method='lm', formula =y~x) +
  facet_wrap(~idh_codel_subtype, scale = "free")  

# Build a multivariate model with hypermutation_status, subtype, and possibly treatment.
prop_change_fit = lm(prop_change~log10(surgical_interval) + idh_codel_subtype, data = gatk_seg_diff_gold)
summary(prop_change_fit)

# TOTAL Stratified by subtype.
IDHcodel_prop_change_fit = lm(prop_change~log10(surgical_interval), data = gatk_seg_IDH_codel)
summary(IDHcodel_prop_change_fit)
IDHnoncodel_prop_change_fit = lm(prop_change~log10(surgical_interval), data = gatk_seg_IDH_noncodel)
summary(IDHnoncodel_prop_change_fit)
# Significant change in the IDHwt tumros.
IDHwt_prop_change_fit = lm(prop_change~log10(surgical_interval), data = gatk_seg_IDHwt)
summary(IDHwt_prop_change_fit)

# AMP Stratified by subtype.
IDHcodel_amp_change_fit = lm(prop_amp~log10(surgical_interval), data = gatk_seg_IDH_codel)
summary(IDHcodel_amp_change_fit)
IDHnoncodel_amp_change_fit = lm(prop_amp~log10(surgical_interval), data = gatk_seg_IDH_noncodel)
summary(IDHnoncodel_amp_change_fit)
IDHwt_amp_change_fit = lm(prop_amp~log10(surgical_interval), data = gatk_seg_IDHwt)
summary(IDHwt_amp_change_fit)

# DEL Stratified by subtype.
IDHcodel_del_change_fit = lm(prop_change~log10(surgical_interval), data = gatk_seg_IDH_codel)
summary(IDHcodel_del_change_fit)
IDHnoncodel_del_change_fit = lm(prop_change~log10(surgical_interval), data = gatk_seg_IDH_noncodel)
summary(IDHnoncodel_del_change_fit)
# The association of "deletions" is the likely the LOSS of AMPLIFICATIONS.
IDHwt_del_change_fit = lm(prop_change~log10(surgical_interval), data = gatk_seg_IDHwt)
summary(IDHwt_del_change_fit)

###################################################################
# Time and SNV changes
###################################################################
# Database table the returns `mf_private_b` as the mutational frequency in RECURRENCE ONLY mutations.
tumor_mut_comparison_anno = dbGetQuery(con, "SELECT * FROM analysis.tumor_mut_comparison_anno")

# The data has already been prepared using the silver_set. 
# Histogram of surgical interval between samples with two separate timepoints.
tumor_pair_mf = tumor_mut_comparison_anno %>% 
  left_join(tumor_pairs, by=c("tumor_pair_barcode", "case_barcode", "tumor_barcode_a", "tumor_barcode_b"))

# How is mutation frequency specific to tumor_b distributed?  
ggplot(tumor_pair_mf, aes(mf_private_b)) + geom_histogram() + theme_bw() + xlab("Mutations per Mb") + ylab("Counts") + facet_wrap(~idh_codel_subtype) 

# Transform the `mf_private_b` so that the data is more normally distributed.
ggplot(tumor_pair_mf, aes(log10(mf_private_b))) + geom_histogram() + theme_bw() + xlab("log10(Mutations per Mb)") + ylab("Counts") + facet_wrap(~idh_codel_subtype) 

# Test whether mutation frequency is associated with surgical interval independent of subtype and hypermutation status.
tumor_pair_mf = tumor_pair_mf %>% 
  mutate(log10_mf_private_b = log10(mf_private_b),
         mf_ab = mf_private_a+mf_private_b,
         log10_mf_churn = log10(mf_ab),
         log10_surgical_int = log10(surgical_interval_mo),
         hypermutator_status = recode(hypermutator_status, `0` = "non-hypermutant", `1` = "hypermutant"))
tumor_pair_mf$hypermutator_status = factor(tumor_pair_mf$hypermutator_status, levels=c("non-hypermutant", "hypermutant"))

# Build a linear model using transformed values.
mf_change_fit = lm(log10_mf_private_b~log10_surgical_int + idh_codel_subtype + hypermutator_status, data = tumor_pair_mf)
summary(mf_change_fit)  ### Strong associations with hypermutation (obviously) and a positive association with surgical interval.
mf_churn_fit = lm(log10_mf_churn~log10_surgical_int + idh_codel_subtype + hypermutator_status, data = tumor_pair_mf)
summary(mf_churn_fit)  ### Strong associations with hypermutation (obviously) and a positive association with surgical interval.



# Visualized by subtype.
pdf(file = "/Users/johnsk/Documents/Life-History/glass-analyses/figures/temporal-somatic-alterations/silver-time-private-b-mut-frequency.pdf", height = 4, width = 8, bg = "transparent", useDingbats = FALSE)
ggplot(tumor_pair_mf, aes(x = log10_surgical_int, y = log10_mf_private_b, color=idh_codel_subtype)) +
  geom_point() + theme_bw() + xlab("log10(Surgical interval (months))")+ ylab("log10(Mutations per Mb \n specific to recurrence)") + geom_smooth(method='lm', formula =y~x) +
  facet_wrap(hypermutator_status~idh_codel_subtype) + labs(color = "Glioma subtype") + ggtitle("GLASS silver set (n = 257)")
dev.off()

# Investigate whether mutational "churn" makes a difference.
pdf(file = "/Users/johnsk/Documents/Life-History/glass-analyses/figures/temporal-somatic-alterations/silver-time-mut-frequency-churn.pdf", height = 4, width = 8, bg = "transparent", useDingbats = FALSE)
ggplot(tumor_pair_mf, aes(x = log10_surgical_int, y = log10_mf_churn, color=idh_codel_subtype)) +
  geom_point() + theme_bw() + xlab("log10(Surgical interval (months))")+ ylab("log10(Mutations per Mb \n 'churn')") + geom_smooth(method='lm', formula =y~x) +
  facet_wrap(hypermutator_status~idh_codel_subtype) + labs(color = "Glioma subtype") + ggtitle("GLASS silver set (n = 257)")
dev.off()


# Quickly visualize the associations by subtype and hypermutation status. Overall, data has an positive association.
ggPredict(mf_change_fit, interactive=TRUE)

# Does the relationship hold if you restrict to non-hypermutators?
tumor_pair_mf_nonhyper = tumor_pair_mf %>% 
  filter(hypermutator_status != "hypermutant")

# The positive association between these two attributes remains the same.
nohyper_mf_change_fit = lm(log10_mf_private_b~log10_surgical_int + idh_codel_subtype, data = tumor_pair_mf_nonhyper)
summary(nohyper_mf_change_fit)
# Plot regression lines across WHO subtypes.
ggPredict(nohyper_mf_change_fit, interactive=TRUE)

# Breakdown the data by WHO glioma subtypes.
tumor_pair_mf_IDH_codel = tumor_pair_mf %>% 
  filter(idh_codel_subtype == "IDHmut-codel", !is.na(log10_surgical_int))
tumor_pair_mf_noncodel = tumor_pair_mf %>% 
  filter(idh_codel_subtype == "IDHmut-noncodel", !is.na(log10_surgical_int))
tumor_pair_mf_IDHwt = tumor_pair_mf %>% 
  filter(idh_codel_subtype == "IDHwt", !is.na(log10_surgical_int))

# Subtype-specific model while adjusting for hypermutation event.
IDHcodel_mf_change_fit = lm(log10_mf_private_b~log10_surgical_int + hypermutator_status, data = tumor_pair_mf_IDH_codel)
summary(IDHcodel_mf_change_fit)
IDHnoncodel_mf_change_fit = lm(log10_mf_private_b~log10_surgical_int + hypermutator_status, data = tumor_pair_mf_noncodel)
summary(IDHnoncodel_mf_change_fit)
# The relationship seems to be driven by the larger number of IDHwt tumors.
IDHwt_mf_change_fit = lm(log10_mf_private_b~log10_surgical_int + hypermutator_status, data = tumor_pair_mf_IDHwt)
summary(IDHwt_mf_change_fit)

# Visualizations of the `mf_private_b` and surgical interval stratified by subtype.
ggplot(tumor_pair_mf_IDH_codel, aes(x = log10_surgical_int, y = log10_mf_private_b)) +
  geom_point() + theme_bw() + xlab("log10(surgical interval (months))")+ ylab("log10(Mutations per Mb)") + geom_smooth(method='lm', formula =y~x) +
  facet_wrap(~ hypermutator_status)  
ggplot(tumor_pair_mf_noncodel, aes(x = log10_surgical_int, y = log10_mf_private_b)) +
  geom_point() + theme_bw() + xlab("log10(surgical interval (months))")+ ylab("log10(Mutations per Mb)") + geom_smooth(method='lm', formula =y~x) +
  facet_wrap(~ hypermutator_status)  
ggplot(tumor_pair_mf_IDHwt, aes(x = log10_surgical_int, y = log10_mf_private_b)) +
  geom_point() + theme_bw() + xlab("log10(surgical interval (months))")+ ylab("log10(Mutations per Mb)") + geom_smooth(method='lm', formula =y~x) +
  facet_wrap(~ hypermutator_status)  

##########################################################
# Try subtype and hypermutator-specific correlation tests
#########################################################
# Further sub-divide based on subtype and hypermutation status.
tumor_pair_mf_IDH_codel_hyper = tumor_pair_mf_IDH_codel %>% 
  filter(hypermutator_status == "hypermutant")
tumor_pair_mf_IDH_codel_nonhyper = tumor_pair_mf_IDH_codel %>% 
  filter(hypermutator_status == "non-hypermutant")
tumor_pair_mf_noncodel_hyper = tumor_pair_mf_noncodel %>% 
  filter(hypermutator_status == "hypermutant")
tumor_pair_mf_noncodel_nonhyper = tumor_pair_mf_noncodel %>% 
  filter(hypermutator_status == "non-hypermutant")
tumor_pair_mf_IDHwt_hyper = tumor_pair_mf_IDHwt %>% 
  filter(hypermutator_status == "hypermutant")
tumor_pair_mf_IDHwt_nonhyper = tumor_pair_mf_IDHwt %>% 
  filter(hypermutator_status == "non-hypermutant")

# Non-parametric correlation tests based on subtype and surgical interval.
cor.test(tumor_pair_mf_IDH_codel_hyper$mf_recurrence, tumor_pair_mf_IDH_codel_hyper$surgical_interval_mo, method = "s")
cor.test(tumor_pair_mf_IDH_codel_nonhyper$mf_recurrence, tumor_pair_mf_IDH_codel_nonhyper$surgical_interval_mo, method = "s")
cor.test(tumor_pair_mf_noncodel_hyper$mf_recurrence, tumor_pair_mf_noncodel_hyper$surgical_interval_mo, method = "s")
cor.test(tumor_pair_mf_noncodel_nonhyper$mf_recurrence, tumor_pair_mf_noncodel_nonhyper$surgical_interval_mo, method = "s")
cor.test(tumor_pair_mf_IDHwt_hyper$mf_recurrence, tumor_pair_mf_IDHwt_hyper$surgical_interval_mo, method = "s")
# Again, it's only IDHwt non-hypermutators that really demonstrate a relationship with surgical interval.
cor.test(tumor_pair_mf_IDHwt_nonhyper$mf_recurrence, tumor_pair_mf_IDHwt_nonhyper$surgical_interval_mo, method = "s")


# Separate out subtypes based on hypermutation status.
gatk_seg_IDH_codel_hyper = gatk_seg_IDH_codel %>% 
  filter(hypermutator_status == 1)
gatk_seg_IDH_codel_nonhyper = gatk_seg_IDH_codel %>% 
  filter(hypermutator_status == 0)
gatk_seg_IDH_noncodel_hyper = gatk_seg_IDH_noncodel %>% 
  filter(hypermutator_status == 1)
gatk_seg_IDH_noncodel_nonhyper = gatk_seg_IDH_noncodel %>% 
  filter(hypermutator_status == 0)
gatk_seg_IDHwt_hyper = gatk_seg_IDHwt %>% 
  filter(hypermutator_status == 1)
gatk_seg_IDHwt_nonhyper = gatk_seg_IDHwt %>% 
  filter(hypermutator_status == 0)

# Non-parametric correlation between copy number changes and surgical interval.
cor.test(gatk_seg_IDH_codel_hyper$prop_change, gatk_seg_IDH_codel_hyper$surgical_interval, method = "s")
cor.test(gatk_seg_IDH_codel_nonhyper$prop_change, gatk_seg_IDH_codel_nonhyper$surgical_interval, method = "s")
cor.test(gatk_seg_IDH_noncodel_hyper$prop_change, gatk_seg_IDH_noncodel_hyper$surgical_interval, method = "s")
cor.test(gatk_seg_IDH_noncodel_hyper$prop_change, gatk_seg_IDH_noncodel_hyper$surgical_interval, method = "s")
cor.test(gatk_seg_IDHwt_hyper$prop_change, gatk_seg_IDHwt_hyper$surgical_interval, method = "s")
# Only IDHwt non-hypermutators demonstrate a relationship with 
cor.test(gatk_seg_IDHwt_nonhyper$prop_change, gatk_seg_IDHwt_nonhyper$surgical_interval, method = "s")

