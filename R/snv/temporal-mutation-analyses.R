#######################################################
# Analyse mutational frequency over time.
# Date: 2018.12.31 
# Author: Kevin J.
#######################################################

# Necessary packages:
library(tidyverse)
library(DBI)
library(gridExtra)
library(nlme)
library(EnvStats)

#######################################################
# Establish connection with the database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

# Retrieve the biospecimen_aliquots from the Database.
tumor_pairs = dbReadTable(con,  Id(schema="analysis",table="tumor_pairs"))
aliquots = dbReadTable(con,  Id(schema="biospecimen",table="aliquots"))
cases = dbReadTable(con,  Id(schema="clinical",table="cases"))
surgeries = dbReadTable(con,  Id(schema="clinical",table="surgeries"))
mut_freq = dbReadTable(con,  Id(schema="analysis",table="mutation_freq"))
tumor_mut_compare = dbReadTable(con,  Id(schema="analysis",table="tumor_mut_comparison_anno"))
silver_set = dbGetQuery(con, "SELECT * FROM analysis.silver_set")

## Define query for extracting mutational frequency data.
optimal_mut_freq = "WITH
selected_tumor_pairs AS
(
  SELECT
  tumor_pair_barcode,
  mf1.coverage_adj_mut_freq AS cov_adj_mut_freq_a,
  mf2.coverage_adj_mut_freq AS cov_adj_mut_freq_b,
  mf1.mutation_count AS mut_count_a,
  mf2.mutation_count AS mut_count_b,
  row_number() OVER (PARTITION BY case_barcode ORDER BY surgical_interval_mo DESC, portion_a ASC, portion_b ASC, substring(tumor_pair_barcode from 27 for 3) ASC) AS priority
  FROM analysis.tumor_pairs ps
  LEFT JOIN analysis.blocklist b1 ON b1.aliquot_barcode = ps.tumor_barcode_a
  LEFT JOIN analysis.blocklist b2 ON b2.aliquot_barcode = ps.tumor_barcode_b
  LEFT JOIN analysis.mutation_freq mf1 ON mf1.aliquot_barcode = ps.tumor_barcode_a  -- join with mutation freq to remove hypermutators
  LEFT JOIN analysis.mutation_freq mf2 ON mf2.aliquot_barcode = ps.tumor_barcode_b
  WHERE
  comparison_type = 'longitudinal' AND
  sample_type_b <> 'M1' AND                                                     -- exclude metastatic samples here because this is outside the scope of our study
  --b1.cnv_exclusion = 'allow' AND b2.cnv_exclusion = 'allow' AND
  b1.coverage_exclusion = 'allow' AND b2.coverage_exclusion = 'allow' AND
  b1.clinical_exclusion = 'allow' AND b2.clinical_exclusion = 'allow' --AND
  --mf1.coverage_adj_mut_freq < 10 AND mf2.coverage_adj_mut_freq < 10            -- filter hypermutators
)
SELECT * FROM selected_tumor_pairs stp
WHERE priority = 1"

## Retain all samples including the hypermutators.
mut_freq_pairs <- dbGetQuery(con, optimal_mut_freq)

## Create a flag for the hypermutated samples. The value `10.20` was determined by capturing the output of boxplot outliers (value = delta_mut_freq).
mut_freq_pairs = mut_freq_pairs %>% 
  # Compute the delta_mutation_frequency.
  mutate(delta_mut_freq = cov_adj_mut_freq_b-cov_adj_mut_freq_a) %>% 
  mutate(outlier_primary = ifelse(cov_adj_mut_freq_a > 10.20, "1", "0"),
    hypermutant_recurrence = ifelse(cov_adj_mut_freq_b > 10.20, "1", "0"))

# See whether any samples with `outlier_primary==1` are worth retaining for hypermutation analysis.
questionable_samples = mut_freq_pairs[which(mut_freq_pairs$outlier_primary==1), ]

# It seems like, "GLSS-CU-R007-TP-01-R1-01D-WXS" is a true hypermutator, retain it for the following analyses.
pairs_to_drop = questionable_samples$tumor_pair_barcode[2:9]
mut_freq_pairs_filtered = mut_freq_pairs %>%  # 193 SUBJECTS
  filter(!tumor_pair_barcode%in%pairs_to_drop)

# Create table with necessary covariate information for the subsequent linear models.
mut_freq_annot = mut_freq_pairs_filtered %>% 
  left_join(tumor_pairs, by="tumor_pair_barcode") %>% 
  left_join(cases, by="case_barcode") %>%
  left_join(aliquots, by=c("tumor_barcode_a" = "aliquot_barcode")) %>% 
  left_join(aliquots, by=c("tumor_barcode_b" = "aliquot_barcode")) %>% 
  left_join(surgeries, by=c("sample_barcode.x" = "sample_barcode")) %>% 
  left_join(surgeries, by=c("sample_barcode.y" = "sample_barcode")) %>% 
  mutate_if(bit64::is.integer64, as.double) %>% 
  select(tumor_pair_barcode:hypermutant_recurrence, tumor_barcode_a, tumor_barcode_b, case_age_diagnosis_years, surgical_interval_mo = surgical_interval_mo.x, idh_status = idh_status.y, idh_codel_subtype = idh_codel_subtype.y, 
         treatment_tmz_a = treatment_tmz.x, treatment_tmz_b = treatment_tmz.y, treatment_tmz_cycles_a = treatment_tmz_cycles.x, treatment_tmz_cycles_b = treatment_tmz_cycles.y,
         treatment_tmz_cycles_6 = treatment_tmz_cycles_6.x)

#####################
#       ALL         #
#####################
## Inspect the distributions for both the response and predictor variables.
## In this case, `delta_mut_freq` is not normally distributed. 
# Cannot log transform (generates many negative INF numbers) AND we want to separate the hypermutators from the nonhypermutators.
ggplot(mut_freq_annot, aes(x=delta_mut_freq)) + geom_histogram() + theme_bw() + xlab("recurrent mut. freq - primary mut. freq") + ggtitle("Delta mutation frequency (n=193 patients)")

## Inspect different transformations for the surigcal interval data.
# Non-normal time distribution.
hist(mut_freq_annot$surgical_interval_mo)
# Cube root that avoids complex numbers.
hist(sign(mut_freq_annot$surgical_interval_mo) * abs(mut_freq_annot$surgical_interval_mo)^(1/3))
# Square-root transformation.
hist(sqrt(mut_freq_annot$surgical_interval_mo)) 
# log transformation. ***SELECTED*** transformation, closest to normal distribution.
hist(log(mut_freq_annot$surgical_interval_mo))


#####################
#   HYPERMUTATORS   #
#####################
## Subset to the 31 hypermutators in the GLASS cohort.
hypermutator_annot = mut_freq_annot %>% 
  filter(hypermutant_recurrence == 1)

## Need to perform these analyses adjusting only for IDH status since the codel group is so small.
## IDHmut_codel IDHmut_noncodel  IDHwt_noncodel 
##       2              17              12 

## Plot variables to examine their distributions.
ggplot(hypermutator_annot, aes(x=delta_mut_freq)) + geom_histogram() + theme_bw() + xlab("recurrent mut. freq - primary mut. freq") + ggtitle("Hypermutators")
ggplot(hypermutator_annot, aes(x=log(delta_mut_freq))) + geom_histogram() + theme_bw() + xlab("log(delta_mut_freq)") + ggtitle("Hypermutators (n=31)")
ggplot(hypermutator_annot, aes(x=surgical_interval_mo)) + geom_histogram() + theme_bw() + xlab("surgical interval months") + ggtitle("Hypermutators")
ggplot(hypermutator_annot, aes(x=log(surgical_interval_mo))) + geom_histogram() + theme_bw() + xlab("log(surgical interval months)") + ggtitle("Hypermutators")

# We should remove the sample with the > 30 TMZ cycles (outlier). Brings the sample size down to 30.
ggplot(hypermutator_annot, aes(x=treatment_tmz_cycles_a)) + geom_histogram() + theme_bw() + xlab("tmz cycles") + ggtitle("Hypermutators")

# Filter out sample with the high number of TMZ cycles.
hypermutator_annot = hypermutator_annot %>% 
  filter(tumor_pair_barcode != "GLSS-MD-0030-TP-01-R1-01D-WGS")

## Various linear models with transformed and untransformed variables.
hyper_out_1 <- lm(log(delta_mut_freq) ~ log(treatment_tmz_cycles_a),
                 data = hypermutator_annot)
summary(hyper_out_1) 
hyper_out_2 <- lm(log(delta_mut_freq) ~ log(surgical_interval_mo)*log(treatment_tmz_cycles_a),
                  data = hypermutator_annot)
summary(hyper_out_2) # When you use an interaction term, the P-values trend toward significance.

## Assess whether there is a difference between the two models.
anova(hyper_out_1, hyper_out_2) # No, there is not.

### Conclusion: In hypermutators, there is a WEAK positive association between surgical interval and delta_mut_freq, WHILE
## there is a WEAK negative relationship between the interaction of surgical interval and treatment tmz cycles. See results from
## 'hyper_out_2'.

## Separate out the IDH wildtype and mutants. Separately run those models.
hypermutator_IDHwt = hypermutator_annot %>% 
  filter(idh_status == "IDHwt")
hypermutator_IDHmt = hypermutator_annot %>% 
  filter(idh_status == "IDHmut")

# Examine the relationship in just the patients with IDHwt tumors.
hyper_out_wt <- lm(log(delta_mut_freq) ~ log(surgical_interval_mo)*log(treatment_tmz_cycles_a),
                  data = hypermutator_IDHwt)
summary(hyper_out_wt) 
# Now for the IDH-mutant tumors.
hyper_out_mt <- lm(log(delta_mut_freq) ~ log(surgical_interval_mo)*log(treatment_tmz_cycles_a),
                   data = hypermutator_IDHmt)
summary(hyper_out_mt) 

# No relationship when separating the models out.

## Plot mutational frequency vs. surgical interval months.
ggplot(hypermutator_annot, aes(x=log(surgical_interval_mo), y=log(delta_mut_freq))) + geom_point() + geom_smooth(method="lm") +
  theme_bw() + xlab("log(surgical interval months)") + ylab("log(Delta Mutation Frequency)") + ggtitle("hypermutators (n=30)")
# TMZ cycles.
ggplot(hypermutator_annot, aes(x=treatment_tmz_cycles_a, y=log(delta_mut_freq))) + geom_point() + geom_smooth(method="lm") +
  theme_bw() + xlab("Number of TMZ cycles") + ylab("log(Delta Mutation Frequency)") + ggtitle("hypermutators (n=18)")

# Too many samples missing and one patient that received > 200 Grays. Do NOT radiation include in analysis.


############################
#   HYPERMUTATORS REMOVED  #
############################
# Variable defined as mutation frequency greater than 10.2.
nonhypermutator_annot = mut_freq_annot %>% 
  filter(hypermutant_recurrence == 0)

# Test whether there is a difference between the non-hypermutators and the hypermutators for surgical interval.
# In the section above I defined "hypermutator_IDHwt" and "hypermutator_IDHmt"
nonhypermutator_IDHwt = nonhypermutator_annot %>% 
  filter(idh_status == "IDHwt")
nonhypermutator_IDHmt = nonhypermutator_annot %>% 
  filter(idh_status == "IDHmut")

wilcox.test(nonhypermutator_IDHwt$surgical_interval_mo, hypermutator_IDHwt$surgical_interval_mo)
wilcox.test(nonhypermutator_IDHmt$surgical_interval_mo, hypermutator_IDHmt$surgical_interval_mo)

test = mut_freq_annot %>% 
  group_by(idh_codel_subtype, hypermutant_recurrence) %>% 
  tally()

q0 = ggplot(test, aes(x = idh_codel_subtype, y = n, fill = hypermutant_recurrence)) + 
  geom_bar(position = "fill",stat = "identity") + geom_text(label = "n=24", x=1, y=0) + geom_text(label = "n=69", x=2, y=0) + geom_text(label = "n=100", x=3, y=0) + theme_bw() + ylab("Proportion") + xlab("") +
  labs(fill='Hypermutator') 

q1 = ggplot(nonhypermutator_IDHwt, aes(x=surgical_interval_mo)) + geom_histogram() + theme_bw() + xlab("surgical interval (months)") + xlim(-5, 90) + ggtitle("IDHwt nonhypermutator")
q2 = ggplot(hypermutator_IDHwt, aes(x=surgical_interval_mo)) + geom_histogram() + geom_text(label = "P-value = 0.09",  x = 75, y = 1.5) + theme_bw() + xlab("surgical interval (months)") + xlim(-5, 90) + ggtitle("IDHwt hypermutator") 
grid.arrange(q1, q2, nrow = 2)

q3 = ggplot(nonhypermutator_IDHmt, aes(x=surgical_interval_mo)) + geom_histogram() + theme_bw() + xlab("surgical interval (months)") + xlim(-5, 160) + ggtitle("IDHmt nonhypermutator")
q4 = ggplot(hypermutator_IDHmt, aes(x=surgical_interval_mo)) + geom_histogram() + theme_bw()  + geom_text(label = "P-value = 0.46",  x = 115, y = 2.25) + xlim(-5, 160) + xlab("surgical interval (months)") + ggtitle("IDHmt hypermutator")
grid.arrange(q3, q4, nrow = 2)





# Create untransformed histogram of mutation frequency.
ggplot(nonhypermutator_annot, aes(x=delta_mut_freq)) + geom_histogram() + theme_bw() + xlab("recurrent mut. freq - primary mut. freq") + ggtitle("Delta mutation frequency (n=162 patients)")

# Examine the distributions of mutation frequency for tumor_a and tumor_b.
ggplot(nonhypermutator_annot, aes(x=cov_adj_mut_freq_a)) + geom_histogram() + theme_bw() + xlab("primary mut. freq") + ggtitle("(n=162 patients)")
ggplot(nonhypermutator_annot, aes(x=cov_adj_mut_freq_b)) + geom_histogram() + theme_bw() + xlab("recurrent mut. freq") + ggtitle("(n=162 patients)")
# Test for all subtypes combined.
wilcox.test(nonhypermutator_annot$cov_adj_mut_freq_a, nonhypermutator_annot$cov_adj_mut_freq_b, paired = TRUE, alternative = "two.sided")

# IDHwt
nonhypermutator_IDHwt = nonhypermutator_annot %>% 
  filter(idh_codel_subtype == "IDHwt_noncodel")
q1 = ggplot(nonhypermutator_IDHwt, aes(x=delta_mut_freq)) + geom_histogram() + theme_bw() + xlab("recurrent mut. freq - primary mut. freq") + ggtitle("Delta mutation frequency - IDHwt")
wilcox.test(nonhypermutator_IDHwt$cov_adj_mut_freq_a, nonhypermutator_IDHwt$cov_adj_mut_freq_b, paired = TRUE, alternative = "two.sided")
# IDHmut noncodel
nonhypermutator_IDHmut_noncodel = nonhypermutator_annot %>% 
  filter(idh_codel_subtype == "IDHmut_noncodel")
q2 = ggplot(nonhypermutator_IDHmut_noncodel, aes(x=delta_mut_freq)) + geom_histogram() + theme_bw() + xlab("recurrent mut. freq - primary mut. freq") + ggtitle("Delta mutation frequency - IDHmut noncodel")
wilcox.test(nonhypermutator_IDHmut_noncodel$cov_adj_mut_freq_a, nonhypermutator_IDHmut_noncodel$cov_adj_mut_freq_b, paired = TRUE, alternative = "two.sided")
# IDHmut codel
nonhypermutator_IDHmut_codel = nonhypermutator_annot %>% 
  filter(idh_codel_subtype == "IDHmut_codel")
q3 = ggplot(nonhypermutator_IDHmut_codel, aes(x=delta_mut_freq)) + geom_histogram() + theme_bw() + xlab("recurrent mut. freq - primary mut. freq") + ggtitle("Delta mutation frequency - IDHmut codel")
wilcox.test(nonhypermutator_IDHmut_codel$cov_adj_mut_freq_a, nonhypermutator_IDHmut_codel$cov_adj_mut_freq_b, paired = TRUE, alternative = "two.sided")
grid.arrange(q1, q2, q3, nrow = 3)


# How many samples have a higher mutation frequency at recurrence?
sum(nonhypermutator_annot$delta_mut_freq>0)

## EXCLUDING HYPERMUTATORS - table(nonhypermutator_annot$idh_codel_subtype)
# IDHmut_codel IDHmut_noncodel  IDHwt_noncodel 
#     22              52              88

## Plot variables to examine their distributions.
# Normally distributed delta_mut_freq in the non-hypermutators. No need to log transform.
ggplot(nonhypermutator_annot, aes(x=delta_mut_freq)) + geom_histogram() + theme_bw() + xlab("Delta mutation frequency") + ggtitle("Non-hypermutators (n=162)")
ggplot(nonhypermutator_annot, aes(x=log(surgical_interval_mo))) + geom_histogram() + theme_bw() + xlab("log(surgical interval months)") + ggtitle("Non-hypermutators (n=162)")
ggplot(nonhypermutator_annot, aes(x=log(treatment_tmz_cycles))) + geom_histogram() + theme_bw() + xlab("log(treatment cycles)") + ggtitle("Non-hypermutators (n=69)")

## Statistical models:
# No significant relationship with transformed surgical interval.
nonhyper_out_1 <- lm(delta_mut_freq ~ log(surgical_interval_mo),
           data = nonhypermutator_annot)
summary(nonhyper_out_1)

# Also, does not hold when using an interaction term for subtype. 
nonhyper_out_2 <- lm(delta_mut_freq ~ log(surgical_interval_mo) * idh_status,
            data = nonhypermutator_annot)
summary(nonhyper_out_2)
# Test if the model is better explained by including IDH status.
anova(nonhyper_out_1, nonhyper_out_2)

# No significant relationship with interval when it is transformed and treatment information is included.
nonhyper_out_3 <- lm(delta_mut_freq ~ log(surgical_interval_mo) * idh_status + log(treatment_tmz_cycles),
            data = nonhypermutator_annot)
summary(nonhyper_out_3)

## CONCLUSIONS: No relationship between non-hypermutators and surgical interval OR TMZ cycles.

## There could be an association if the models are separately run.
## Separate out the IDH wildtype and mutants.
nonhypermutator_IDHwt = nonhypermutator_annot %>% 
  filter(idh_status == "IDHwt")
nonhypermutator_IDHmt = nonhypermutator_annot %>% 
  filter(idh_status == "IDHmut")

# Restrict analysis by IDH mutation status.
nonhyper_idh_wt_out <- lm(delta_mut_freq ~ log(surgical_interval_mo)*log(treatment_tmz_cycles_a),
                     data = nonhypermutator_IDHwt)
# A stronger relationship with TMZ cycle, but still does NOT reach significance.
summary(nonhyper_idh_wt_out) 

##  Plot mutational frequency vs. surgical interval months.
ggplot(nonhypermutator_IDHwt, aes(x=log(surgical_interval_mo), y=log(treatment_tmz_cycles_a))) + geom_point() + geom_smooth(method="lm") +
  theme_bw() + xlab("log(surgical interval months)") + ylab("log(TMZ cycles)") + ggtitle("IDHwt nonhypermutator (n=53)")
ggplot(nonhypermutator_IDHwt, aes(x=log(surgical_interval_mo), y=delta_mut_freq)) + geom_point() + geom_smooth(method="lm") +
  theme_bw() + xlab("log(TMZ cycles)") + ylab("log(Delta Mutation Frequency)") + ggtitle("IDHwt nonhypermutator (n=53)")

# Run the same analyses for the IDHmut non-hypermutated tumors.
nonhyper_idh_mt_out <- lm(delta_mut_freq ~ log(surgical_interval_mo)*log(treatment_tmz_cycles_a),
                          data = nonhypermutator_IDHmt)
# No association between time or treatment.
summary(nonhyper_idh_mt_out)

###############################
# Mutation Freq. classified by
# private_a, shared, private_b
###############################
# This query calculates the mutation frequency for private to tumor_a (uses tumor_a coverage), shared (uses minimum coverage of tumor_a and tumor_b), 
# and private to tumor_b (uses tumor_b coverage).
mut_freq_tumor_type <- dbGetQuery(con, read_file("sql/mutation_freq_private_shared.sql"))
clindata <- dbGetQuery(con, "SELECT DISTINCT case_barcode, idh_codel_subtype FROM clinical.surgeries WHERE idh_codel_subtype IS NOT NULL")

# Remove the hypermutators from this analysis.
mut_freq_tumor_type <- mut_freq_tumor_type %>% left_join(clindata) %>% filter(mf_a < 10.2, mf_b < 10.2) %>% 
  left_join(cases, by="case_barcode") %>% 
  mutate(case_age_diagnosis_months = case_age_diagnosis_years*12,
         case_age_recurrence_months = case_age_diagnosis_months+surgical_interval_mo)

# Plot the log transformations. These transformed data should be approximately normal.
p1 = ggplot(mut_freq_tumor_type, aes(x=log(mf_shared))) + geom_histogram() + theme_bw() + xlab("log(mut_freq shared)") + xlim(-5, 2.5) + ggtitle("Mutation Frequency - Shared (n=171)")
p2 = ggplot(mut_freq_tumor_type, aes(x=log(mf_private_a))) + geom_histogram() + theme_bw() + xlab("log(mut_freq private_a)") + xlim(-5, 2.5) + ggtitle("Mutation Frequency - Private A (n=171)")
p3 = ggplot(mut_freq_tumor_type, aes(x=log(mf_private_b))) + geom_histogram() + theme_bw() + xlab("log(mut_freq private_b)") + xlim(-5, 2.5) + ggtitle("Mutation Frequency - Private B (n=171)")
grid.arrange(p1, p2, p3, nrow = 3)

# Again, for this set of data we have a non-normal surgical interval value.
ggplot(mut_freq_tumor_type, aes(x=surgical_interval_mo)) + geom_histogram() + theme_bw() + xlab("") 
# Log transformation of surgical interval yields a more gaussian distribution.
ggplot(mut_freq_tumor_type, aes(x=log(surgical_interval_mo))) + geom_histogram() + theme_bw() + xlab("log(surgical interval months)") 

# For the 'mf_shared' analysis, it would be helpful to remove those samples with ZERO shared mutations otherwise it errors (log(0) = -Inf).
mut_freq_tumor_type_sub = mut_freq_tumor_type %>% 
  filter(intersection_ab > 0)

# Transforming the two numeric dependent and independent variables to be normally distrbuted.
mf_shared_out_1 <- lm(log(mf_shared) ~ log(surgical_interval_mo),
                              data = mut_freq_tumor_type_sub)
summary(mf_shared_out_1) # Without considering subtype there is a very strong NEGATIVE association.
mf_privateA_out_2 <- lm(log(mf_private_a) ~ log(surgical_interval_mo),
                       data = mut_freq_tumor_type)
summary(mf_privateA_out_2) # No association.
mf_privateB_out_3a <- lm(log(mf_private_b) ~ log(surgical_interval_mo),
                       data = mut_freq_tumor_type)
summary(mf_privateB_out_3a) # Without considering subtype there is a strong POSITIVE association.
mf_privateB_out_3b <- lm(log(mf_private_b) ~ log(surgical_interval_mo)*idh_codel_subtype,
                        data = mut_freq_tumor_type)
summary(mf_privateB_out_3b) # If subtype is included as an effect modifier this relationship is no longer significant.

## Treat each subtype analysis as if it were a completely separate hypothesis.
## RATIONALE: Easier to interpret for subtype-specific inferences.

#### IDHwt
mut_freq_IDHwt = mut_freq_tumor_type %>% 
  filter(idh_codel_subtype=="IDHwt_noncodel") %>% 
  filter(intersection_ab > 0)

# Summary for each mutation frequency type. 
## SURGICAL INTERVAL.
idh_wt_out_1 <- lm(log(mf_private_a) ~ log(surgical_interval_mo),
                   data = mut_freq_IDHwt)
summary(idh_wt_out_1) # POSITIVE association with surgical interval.
idh_wt_out_2 <- lm(log(mf_private_b) ~ log(surgical_interval_mo),
                   data = mut_freq_IDHwt)
summary(idh_wt_out_2) # POSITIVE association with surgical interval.
idh_wt_out_3 <- lm(log(mf_shared) ~ log(surgical_interval_mo),
                   data = mut_freq_IDHwt)
summary(idh_wt_out_3) # NO association with surgical interval.

## Age at Diagnosis.
idh_wt_out_1 <- lm(log(mf_private_a) ~ case_age_diagnosis_months,
                   data = mut_freq_IDHwt)
summary(idh_wt_out_1) # No association with age at diagnosis.
idh_wt_out_2 <- lm(log(mf_private_b) ~ case_age_diagnosis_months,
                   data = mut_freq_IDHwt)
summary(idh_wt_out_2) # No association with age at diagnosis.
idh_wt_out_3 <- lm(log(mf_shared) ~ case_age_diagnosis_months,
                   data = mut_freq_IDHwt)
summary(idh_wt_out_3) # Statistically significant WEAK association with age at diagnosis.

## Age at Recurrence.
idh_wt_out_1 <- lm(log(mf_private_a) ~ case_age_recurrence_months,
                   data = mut_freq_IDHwt)
summary(idh_wt_out_1) # No association with age at recurrence
idh_wt_out_2 <- lm(log(mf_private_b) ~ case_age_recurrence_months,
                   data = mut_freq_IDHwt)
summary(idh_wt_out_2) # No association with age at recurrence
idh_wt_out_3 <- lm(log(mf_shared) ~ case_age_recurrence_months,
                   data = mut_freq_IDHwt)
summary(idh_wt_out_3) # Statistically significant WEAK association with age at recurrence.

#### IDHmut NON-codels
mut_freq_IDHmt = mut_freq_tumor_type %>% 
  filter(idh_codel_subtype=="IDHmut_noncodel")

## SURGICAL INTERVAL.
idh_mt_out_1 <- lm(log(mf_private_a) ~ log(surgical_interval_mo),
                 data = mut_freq_IDHmt)
summary(idh_mt_out_1) # No association.
idh_mt_out_2 <- lm(log(mf_private_b) ~ log(surgical_interval_mo),
                 data = mut_freq_IDHmt)
summary(idh_mt_out_2) # Weak POSITIVE association.
idh_mt_out_3 <- lm(log(mf_shared) ~ log(surgical_interval_mo),
                   data = mut_freq_IDHmt)
summary(idh_mt_out_3) # No association.

## Age at Diagnosis.
idh_mt_out_1 <- lm(log(mf_private_a) ~ case_age_diagnosis_months,
                   data = mut_freq_IDHmt)
summary(idh_mt_out_1) # No association with age at diagnosis.
idh_mt_out_2 <- lm(log(mf_private_b) ~ case_age_diagnosis_months,
                   data = mut_freq_IDHmt)
summary(idh_mt_out_2) # No association with age at diagnosis.
idh_mt_out_3 <- lm(log(mf_shared) ~ case_age_diagnosis_months,
                   data = mut_freq_IDHmt)
summary(idh_mt_out_3) # Statistically significant WEAK association with age at diagnosis.

## Age at Recurrence.
idh_mt_out_1 <- lm(log(mf_private_a) ~ case_age_recurrence_months,
                   data = mut_freq_IDHmt)
summary(idh_mt_out_1) # No association with age at recurrence
idh_mt_out_2 <- lm(log(mf_private_b) ~ case_age_recurrence_months,
                   data = mut_freq_IDHmt)
summary(idh_mt_out_2) # No association with age at recurrence
idh_mt_out_3 <- lm(log(mf_shared) ~ case_age_recurrence_months,
                   data = mut_freq_IDHmt)
summary(idh_mt_out_3) # Statistically significant WEAK association with age at recurrence.

## IDHmut CODELs. 
mut_freq_codel = mut_freq_tumor_type %>% 
  filter(idh_codel_subtype=="IDHmut_codel")

## SURGICAL INTERVAL
# Summary of the linear models. For all three models, there was a weak negative association between mf_* and surgical interval.
idh_codel_out_1 <- lm(log(mf_private_a) ~ log(surgical_interval_mo),
                      data = mut_freq_codel)
summary(idh_codel_out_1) # No association.
idh_codel_out_2 <- lm(log(mf_private_b) ~ log(surgical_interval_mo),
                 data = mut_freq_codel)
summary(idh_codel_out_2) # No association.
idh_codel_out_3 <- lm(log(mf_shared) ~ log(surgical_interval_mo),
                      data = mut_freq_codel)
summary(idh_codel_out_3) # No associaiton.

## Age at Diagnosis.
idh_codel_out_1 <- lm(log(mf_private_a) ~ case_age_diagnosis_months,
                   data = mut_freq_codel)
summary(idh_codel_out_1) # Significant positive association with age at diagnosis.
idh_codel_out_2 <- lm(log(mf_private_b) ~ case_age_diagnosis_months,
                   data = mut_freq_codel)
summary(idh_codel_out_2) # No association with age at diagnosis.
idh_codel_out_3 <- lm(log(mf_shared) ~ case_age_diagnosis_months,
                   data = mut_freq_codel)
summary(idh_codel_out_3) # Significant positive association with age at diagnosis.

## Age at Recurrence.
idh_codel_out_1 <- lm(log(mf_private_a) ~ case_age_recurrence_months,
                   data = mut_freq_codel)
summary(idh_codel_out_1) # Significant positive association with age at recurrence
idh_codel_out_2 <- lm(log(mf_private_b) ~ case_age_recurrence_months,
                   data = mut_freq_codel)
summary(idh_codel_out_2) # No association with age at recurrence
idh_codel_out_3 <- lm(log(mf_shared) ~ case_age_recurrence_months,
                   data = mut_freq_codel)
summary(idh_codel_out_3) # Statistically significant positive association with age at recurrence.


# Figures for each mutation type. SHARED.
ggplot(mut_freq_tumor_type_sub, aes(x = log(surgical_interval_mo), y = log(mf_shared))) + 
  geom_point() + ggtitle("shared mutations - mut. freq") +
  geom_smooth(method="lm") +
  facet_wrap(~idh_codel_subtype, scales = "free")
# PRIVATE to tumor A.
ggplot(mut_freq_tumor_type, aes(x = log(surgical_interval_mo), y = log(mf_private_a))) + 
  geom_point() + ggtitle("private_a mutations - mut. freq") +
  geom_smooth(method="lm") +
  facet_wrap(~idh_codel_subtype, scales = "free")
# PRIVATE to tumor B.
ggplot(mut_freq_tumor_type, aes(x = log(surgical_interval_mo), y = log(mf_private_b))) + 
  geom_point() + ggtitle("private_b mutations - mut. freq") +
  geom_smooth(method="lm") +
  facet_wrap(~idh_codel_subtype, scales = "free")



##########################################
# Mutational signatures and time analysis
##########################################
# Goal: associate mutational signatures scores with TMZ cycles, TMZ (yes/no),
# subject age at first Dx, and surgical interval.
mutsig = dbReadTable(con,  Id(schema="analysis",table="mutsig"))

# Limit this analysis to mutational signatures applied to all mutations in a specific aliquot.
mut_sig_meta = mutsig %>% 
  left_join(aliquots, by=c("tumor_pair_barcode" = "aliquot_barcode")) %>% 
  left_join(surgeries, by="sample_barcode") %>%
  left_join(cases, by="case_barcode") %>% 
  mutate(sample_type = substr(tumor_pair_barcode, 14, 15))

# Some visualizations of all tumors and their relationship with age at first Dx.
ggplot(mut_sig_meta, aes(x=case_age_diagnosis_years, y = relative_contribution)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~signature)

# Tabulate the largest signature contributors by subtype.
mut_sig_wxs_summary = mut_sig_meta %>% 
  group_by(signature, idh_codel_subtype) %>% 
  summarise(total_contribution = sum(relative_contribution),
            sample_number = n())

# Stacked barplot by sample. May reduce low numbers to one group to improve calriyty.
ggplot(mut_sig_meta, aes(x=tumor_pair_barcode, y=relative_contribution, fill=signature)) + geom_bar(stat="identity") +
  theme(axis.text.x=element_blank()) + labs(color="Signature") + xlab("") + facet_wrap(~idh_codel_subtype, scales = "free")

# Just examining mutational signature 1.
mut_sig1 = mut_sig_meta %>% 
  filter(signature=="Signature.1") 
ggplot(mut_sig1, aes(x=case_age_diagnosis_years, y = relative_contribution)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~surgery_number + idh_codel_subtype)

# Use the mutational signatures applied to private_a, shared, and private_b.
mutsig_private_vs_shared = dbReadTable(con,  Id(schema="analysis", table="mutsig_private_vs_shared"))

# Collect all of the metadata for this analysis.
mut_sig_type_meta = mutsig_private_vs_shared %>% 
  left_join(tumor_pairs, by="tumor_pair_barcode") %>% 
  left_join(cases, by="case_barcode") %>% 
  left_join(aliquots, by=c("tumor_barcode_a" = "aliquot_barcode")) %>% 
  left_join(surgeries, by="sample_barcode") %>% 
  select(tumor_pair_barcode:surgical_interval_mo.x, case_age_diagnosis_years, idh_status, idh_codel_subtype, treatment_tmz:treatment_tmz_cycles_6)

# In general, signature 1 always seems dominant. Signature 11 lights up at recurrence.
mut_sig_primary = mut_sig_type_meta %>% 
  filter(mutation_status == 'primary') %>% 
  group_by(signature, idh_codel_subtype) %>% 
  summarise(total_contribution = sum(relative_contribution),
            sample_number = n())
mut_sig_shared = mut_sig_type_meta %>% 
  filter(mutation_status == 'shared') %>% 
  group_by(signature, idh_codel_subtype) %>% 
  summarise(total_contribution = sum(relative_contribution),
            sample_number = n())
mut_sig_recurrent = mut_sig_type_meta %>% 
  filter(mutation_status == 'recurrent') %>% 
  group_by(signature, idh_codel_subtype) %>% 
  summarise(total_contribution = sum(relative_contribution),
            sample_number = n())

## Focus on signature 1 first.
mut_sig1_primary = mut_sig_type_meta %>% 
  filter(mutation_status == "primary") %>% 
  filter(signature =="Signature.1")
mut_sig1_shared = mut_sig_type_meta %>% 
  filter(mutation_status == "shared") %>% 
  filter(signature =="Signature.1")
mut_sig1_recurrent = mut_sig_type_meta %>% 
  filter(mutation_status == "recurrent") %>% 
  filter(signature =="Signature.1")
ggplot(mut_sig1_primary, aes(x=case_age_diagnosis_years, y = relative_contribution)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~idh_codel_subtype) + ggtitle("Primary")
ggplot(mut_sig1_shared, aes(x=case_age_diagnosis_years, y = relative_contribution)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~idh_codel_subtype) + ggtitle("Shared")
ggplot(mut_sig1_recurrent, aes(x=case_age_diagnosis_years, y = relative_contribution)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~idh_codel_subtype) + ggtitle("Recurrent")

## Signature 11.
mut_sig11_primary = mut_sig_type_meta %>% 
  filter(mutation_status == "primary") %>% 
  filter(signature =="Signature.11")
mut_sig11_shared = mut_sig_type_meta %>% 
  filter(mutation_status == "shared") %>% 
  filter(signature =="Signature.11")
mut_sig11_recurrent = mut_sig_type_meta %>% 
  filter(mutation_status == "recurrent") %>% 
  filter(signature =="Signature.11")

# There is NOT sufficient information for the 'treatment_tmz_cycles' variable.
ggplot(mut_sig11_primary, aes(x=log(treatment_tmz_cycles), y = relative_contribution)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~idh_codel_subtype) + ggtitle("Primary")
ggplot(mut_sig11_shared, aes(x=log(treatment_tmz_cycles), y = relative_contribution)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~idh_codel_subtype) + ggtitle("Shared")
ggplot(mut_sig11_recurrent, aes(x=log(treatment_tmz_cycles), y = relative_contribution)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~idh_codel_subtype) + ggtitle("Recurrent")

ggplot(mut_sig11_primary, aes(x=treatment_tmz, y = relative_contribution)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~idh_codel_subtype) + ggtitle("Primary")
ggplot(mut_sig11_shared, aes(x=treatment_tmz, y = relative_contribution)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~idh_codel_subtype) + ggtitle("Shared")
ggplot(mut_sig11_recurrent, aes(x=treatment_tmz, y = relative_contribution)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~idh_codel_subtype) + ggtitle("Recurrent")


#######################################
# OLD APPROACH
# LINEAR MIXED EFFECTS MODEL
# RANDOM EFFECT = PATIENT
#######################################
# Combine to retrieve all necessary data.
mut_freq_annot = mut_freq_pairs %>% 
  left_join(aliquots, by="aliquot_barcode") %>% 
  left_join(surgeries, by="sample_barcode") %>% 
  mutate_if(bit64::is.integer64, as.double) %>% 
  filter(coverage_adj_mut_freq < 10.2) # Remove all hypermutators.

# Visualize important clinical variables: surgical_interval_mo, case_age_diagnosis_years, idh_codel_subtype, and coverage_adj_mut_freq 
mut_freq_codel = mut_freq_annot %>% 
  filter(idh_codel_subtype == "IDHmut_codel")
mut_freq_noncodel = mut_freq_annot %>% 
  filter(idh_codel_subtype == "IDHmut_noncodel")
mut_freq_wt = mut_freq_annot %>% 
  filter(idh_codel_subtype == "IDHwt_noncodel")

# Analyzing both primary and recurrence in this analysis while adjusting for Age @ Dx.
fit_codel = lme(coverage_adj_mut_freq ~ surgical_interval_mo + case_age_diagnosis_years, 
                  random = ~1 | case_barcode, 
                  data = mut_freq_codel,
                  na.action = na.omit)
summary(fit_codel)
# Model for the non-codels.
fit_noncodel = lme(coverage_adj_mut_freq ~ surgical_interval_mo + case_age_diagnosis_years, 
                random = ~1 | case_barcode, 
                data = mut_freq_noncodel,
                na.action = na.omit)
summary(fit_noncodel)
# Model for the IDhwt tumors.
fit_wt = lme(coverage_adj_mut_freq ~ surgical_interval_mo + case_age_diagnosis_years, 
                   random = ~1 | case_barcode, 
                   data = mut_freq_wt,
                   na.action = na.omit)
summary(fit_wt)

# What happens if you include all of the samples in one model and adjust for IDH status confounding.
fit_all = lme(coverage_adj_mut_freq ~ surgical_interval_mo + case_age_diagnosis_years + as.factor(idh_status), 
             random = ~1 | case_barcode, 
             data = mut_freq_annot,
             na.action = na.omit)
summary(fit_all)
