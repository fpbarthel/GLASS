##############################################
# Analyze results from neutralitytestR applied to each aliquot
# Updated: 2019.01.26
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
require(alluvial)
library(ggExtra)
library(EnvStats)

#######################################################
# Establish connection with GLASS database.
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

# Define silver set drivers to perform association tests with neutrality estimates.
silver_drivers = silver_set %>% 
  inner_join(all_drivers, by="tumor_pair_barcode") %>% 
  inner_join(clinal_tumor_pairs, by="tumor_pair_barcode")
silver_drivers$any_driver_stability = ifelse(silver_drivers$snv_driver_stability=="Driver unstable" | silver_drivers$cnv_driver_stability=="Driver unstable", "Driver unstable", "Driver stable")

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

# Determine the number of mutations that are considered subclonal (VAF: 0.1-0.25) in each tumor.
aliquot_mutation_counts = glass_single_vaf %>% 
  filter(vaf >= 0.1 & vaf <= 0.25) %>% 
  group_by(aliquot_barcode) %>% 
  summarize(subclonal_mut = n())

# Although all aliquot results are present here, BUT need to filter low subclonal mutation count [purity < 0.4 AND ploidy > 3].
# Combine the neutrality at the aliquot with specific tumor_pairs.
per_sample_neutrality = tumor_pairs %>% 
  left_join(aliquot_neutrality, by=c("tumor_barcode_a"="aliquot_barcode")) %>% 
  left_join(aliquot_neutrality, by=c("tumor_barcode_b"="aliquot_barcode")) %>% 
  filter(subclonal_mut.x >= 30, subclonal_mut.y >= 30, cellularity.x > 0.4, cellularity.y > 0.4, ploidy.x < 3, ploidy.y < 3) %>% 
  mutate(primary_evolution = ifelse(area_pval.x < 0.05, "S", "N"), 
         recurrence_evolution = ifelse(area_pval.y < 0.05, "S", "N"),
         evolution_mode = paste(primary_evolution, recurrence_evolution, sep="-")) %>% 
  left_join(all_drivers, by = "tumor_pair_barcode") %>% 
  left_join(clinal_tumor_pairs, by = "tumor_pair_barcode") 

# Create a standard. of care treatment group for TMZ (at least 6 cycles).
per_sample_neutrality$tmz_std_care <- ifelse(per_sample_neutrality$received_tmz_sum_cycles >= 6, "1", "0")

# Subset to silver set so that a patient is only represented once at a single time point. N = 98 samples.
silver_neutrality = per_sample_neutrality %>% 
  inner_join(silver_set, by="tumor_pair_barcode") %>% 
  left_join(cases, by=c("case_barcode.x"="case_barcode")) %>% 
  mutate(treatment = ifelse(received_tmz == 1 | received_rt == 1, "YES", "NO"),
         treatment = ifelse(is.na(treatment), "Unknown", treatment))

# Write out table for Roel.
colnames(silver_neutrality)
neutral_out = silver_neutrality %>% 
  select(tumor_pair_barcode, ploidy_p = ploidy.x,  cellularity_p = cellularity.x, ploidy_r = ploidy.y,  cellularity_r = cellularity.y,
         idh_codel_subtype, primary_evolution, recurrence_evolution, evolution_mode)
write.table(neutral_out, "/Users/johnsk/Documents/neutrality_silver_set_lower_confidence.txt", sep="\t", row.names = F, col.names = T, quote = F)

## Color neutral-neutral or selection-selection.
# Use subtype, driver stability, and evolution pattern.
neutral_dat = silver_neutrality %>% 
  select(tumor_pair_barcode, idh_codel_subtype, hypermutator_status, received_tmz, received_rt, treatment, primary_evolution, recurrence_evolution) %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, "IDHmut_codel" = "IDHmut codel", "IDHmut_noncodel" = "IDHmut noncodel", "IDHwt_noncodel" = "IDHwt"),
         primary_evolution = recode(primary_evolution, "N" = "Neutral", "S" = "Selected"),
         recurrence_evolution = recode(recurrence_evolution, "N" = "Neutral", "S" = "Selected"))  %>% 
  group_by(idh_codel_subtype, primary_evolution, treatment, recurrence_evolution) %>% 
  summarise(Freq = n())

# Adjust the alluvial plot so that different colors map to different selection strength.
# Color scheme for evolution trajectory.
#2FB3CA - IDHwt
#F1564F - IDHmut noncodel
#F69654 - IDHmut codel
pal = c("#2FB3CA", "#F1564F", "#F69654", "lightgray")
neutral_dat$col_type = ifelse(neutral_dat$idh_codel_subtype == "IDHwt", pal[1], ifelse(neutral_dat$idh_codel_subtype == "IDHmut noncodel", pal[2], pal[3]))
neutral_dat$col_type = ifelse(neutral_dat$treatment == "Unknown", pal[4], neutral_dat$col_type)
alluvial(neutral_dat[,1:4], freq=neutral_dat$Freq, border=NA,
#       hide = neutral_dat$Freq < quantile(neutral_dat$Freq, .50),
         col= neutral_dat$col_type,
         cw = 0.25, axis_labels = c("IDH subtype", "Primary tumor", "Treatment", "Recurrent tumor"))
# Create legend for alluvial plot.
plot.new()
legend("center", title="IDH 1p/19q subtype",
       c("IDHwt","IDHmut noncodel","IDHmut codel", "Treatment unknown"), fill=pal, cex=1.5, bty = "n")

## Create revised survival variables for Kaplan-Meier curves and log-rank tests.
# The status indicator, normally is 0=alive, 1=dead.
silver_neutrality$patient_vital = ifelse(silver_neutrality$case_vital_status=="alive", 0, 1)

# Create a variable for post-recurrence survival (also presented in months).
silver_neutrality$patient_post_recur_surv = silver_neutrality$case_overall_survival_mo-silver_neutrality$surgical_interval
# Treat all "neutrality" and all "selection" as groups by themselves
silver_neutrality$binary_mode = ifelse(silver_neutrality$evolution_mode=="N-N", "Neutral-Neutral", "N-S|S-N|S-S")
silver_neutrality$s_binary_mode = ifelse(silver_neutrality$evolution_mode=="S-S", "Selection-Selection", "N-S|S-N|N-N")

# Subset to IDHwt, IDHmut noncodel, and IDHmut codel groups for silver set.
silver_neutrality_IDHwt = silver_neutrality %>% 
  filter(idh_codel_subtype %in% c("IDHwt_noncodel"))
silver_neutrality_IDHmut_noncodel = silver_neutrality %>% 
  filter(idh_codel_subtype %in% c("IDHmut_noncodel"))
silver_neutrality_IDHmut_codel = silver_neutrality %>% 
  filter(idh_codel_subtype %in% c("IDHmut_codel"))

# What is the median surgical interval for each subtype x mode of recurrence evolution? Seems to be a difference for IDHwt.
silver_neutrality %>% 
  group_by(recurrence_evolution, idh_codel_subtype) %>% 
  summarise(med_int = median(surgical_interval, na.rm = T),
            sample_size = n())

# Evolution at primary is marginally associated with idh_subtype.
fisher.test(table(silver_neutrality$idh_codel_subtype, silver_neutrality$primary_evolution))
# Evolution at recurrence is marginally associated with idh_subtype.
fisher.test(table(silver_neutrality$idh_codel_subtype, silver_neutrality$recurrence_evolution))
# Evolution transition between tumor_a and tumor_b is not strongly associated with subtype.
fisher.test(table(silver_neutrality$idh_codel_subtype, silver_neutrality$evolution_mode))
# Interesting, the IDHmut noncodels had the highest rate of neutral-neutral.
fisher.test(table(silver_neutrality$idh_codel_subtype, silver_neutrality$binary_mode))
# It seems again that IDHmut noncodels have fewer S-S events proportionally compared with IDHmut codels and IDHwt.
fisher.test(table(silver_neutrality$idh_codel_subtype, silver_neutrality$s_binary_mode))

# Test by method of **EVOLUTION AT PRIMARY** since this could be explained by other clinical variables.
## IDH MUT CODELs.
wilcox.test(silver_neutrality_IDHmut_codel$surgical_interval~silver_neutrality_IDHmut_codel$primary_evolution)
fisher.test(table(silver_neutrality_IDHmut_codel$primary_evolution, silver_neutrality_IDHmut_codel$snv_driver_stability))
fisher.test(table(silver_neutrality_IDHmut_codel$primary_evolution, silver_neutrality_IDHmut_codel$cnv_driver_stability))
fisher.test(table(silver_neutrality_IDHmut_codel$primary_evolution, silver_neutrality_IDHmut_codel$recurrence_location))
fisher.test(table(silver_neutrality_IDHmut_codel$primary_evolution, silver_neutrality_IDHmut_codel$received_tmz))
fisher.test(table(silver_neutrality_IDHmut_codel$primary_evolution, silver_neutrality_IDHmut_codel$tmz_std_care))
fisher.test(table(silver_neutrality_IDHmut_codel$primary_evolution, silver_neutrality_IDHmut_codel$received_rt))
fisher.test(table(silver_neutrality_IDHmut_codel$primary_evolution, silver_neutrality_IDHmut_codel$treatment))
fisher.test(table(silver_neutrality_IDHmut_codel$primary_evolution, silver_neutrality_IDHmut_codel$hypermutator_status))
fisher.test(table(silver_neutrality_IDHmut_codel$primary_evolution, silver_neutrality_IDHmut_codel$grade_change))

# IDH MUT NON-CODELs.
wilcox.test(silver_neutrality_IDHmut_noncodel$surgical_interval~silver_neutrality_IDHmut_noncodel$primary_evolution)
fisher.test(table(silver_neutrality_IDHmut_noncodel$primary_evolution, silver_neutrality_IDHmut_noncodel$snv_driver_stability))
# Interesting, weak association between CNV driver stability and primary evolution.
fisher.test(table(silver_neutrality_IDHmut_noncodel$primary_evolution, silver_neutrality_IDHmut_noncodel$cnv_driver_stability))
fisher.test(table(silver_neutrality_IDHmut_noncodel$primary_evolution, silver_neutrality_IDHmut_noncodel$recurrence_location))
fisher.test(table(silver_neutrality_IDHmut_noncodel$primary_evolution, silver_neutrality_IDHmut_noncodel$received_tmz))
fisher.test(table(silver_neutrality_IDHmut_noncodel$primary_evolution, silver_neutrality_IDHmut_noncodel$tmz_std_care))
fisher.test(table(silver_neutrality_IDHmut_noncodel$primary_evolution, silver_neutrality_IDHmut_noncodel$received_rt))
fisher.test(table(silver_neutrality_IDHmut_noncodel$primary_evolution, silver_neutrality_IDHmut_noncodel$treatment))
fisher.test(table(silver_neutrality_IDHmut_noncodel$primary_evolution, silver_neutrality_IDHmut_noncodel$hypermutator_status))
fisher.test(table(silver_neutrality_IDHmut_noncodel$primary_evolution, silver_neutrality_IDHmut_noncodel$grade_change))


# IDH WT - neutrality is association with surgical interval.
wilcox.test(silver_neutrality_IDHwt$surgical_interval~silver_neutrality_IDHwt$primary_evolution)
fisher.test(table(silver_neutrality_IDHwt$primary_evolution, silver_neutrality_IDHwt$snv_driver_stability))
fisher.test(table(silver_neutrality_IDHwt$primary_evolution, silver_neutrality_IDHwt$cnv_driver_stability))
fisher.test(table(silver_neutrality_IDHwt$primary_evolution, silver_neutrality_IDHwt$recurrence_location))
fisher.test(table(silver_neutrality_IDHwt$primary_evolution, silver_neutrality_IDHwt$received_tmz))
fisher.test(table(silver_neutrality_IDHwt$primary_evolution, silver_neutrality_IDHwt$tmz_std_care))
fisher.test(table(silver_neutrality_IDHwt$primary_evolution, silver_neutrality_IDHwt$received_rt))
fisher.test(table(silver_neutrality_IDHwt$primary_evolution, silver_neutrality_IDHwt$treatment))
fisher.test(table(silver_neutrality_IDHwt$primary_evolution, silver_neutrality_IDHwt$hypermutator_status))
fisher.test(table(silver_neutrality_IDHwt$primary_evolution, silver_neutrality_IDHwt$grade_change))
wilcox.test(silver_neutrality_IDHwt$snv_driver_count~silver_neutrality_IDHwt$primary_evolution)
wilcox.test(silver_neutrality_IDHwt$snv_driver_count_private_a~silver_neutrality_IDHwt$primary_evolution)


# Test by method of **EVOLUTION AT RECURRENCE** since this could be explained by other clinical variables.
## IDH MUT CODELs.
wilcox.test(silver_neutrality_IDHmut_codel$surgical_interval~silver_neutrality_IDHmut_codel$recurrence_evolution)
fisher.test(table(silver_neutrality_IDHmut_codel$recurrence_evolution, silver_neutrality_IDHmut_codel$snv_driver_stability))
fisher.test(table(silver_neutrality_IDHmut_codel$recurrence_evolution, silver_neutrality_IDHmut_codel$cnv_driver_stability))
fisher.test(table(silver_neutrality_IDHmut_codel$recurrence_evolution, silver_neutrality_IDHmut_codel$recurrence_location))
fisher.test(table(silver_neutrality_IDHmut_codel$recurrence_evolution, silver_neutrality_IDHmut_codel$received_tmz))
fisher.test(table(silver_neutrality_IDHmut_codel$recurrence_evolution, silver_neutrality_IDHmut_codel$tmz_std_care))
fisher.test(table(silver_neutrality_IDHmut_codel$recurrence_evolution, silver_neutrality_IDHmut_codel$received_rt))
fisher.test(table(silver_neutrality_IDHmut_codel$recurrence_evolution, silver_neutrality_IDHmut_codel$treatment))
fisher.test(table(silver_neutrality_IDHmut_codel$recurrence_evolution, silver_neutrality_IDHmut_codel$hypermutator_status))
fisher.test(table(silver_neutrality_IDHmut_codel$recurrence_evolution, silver_neutrality_IDHmut_codel$grade_change))

# IDH MUT NON-CODELs.
wilcox.test(silver_neutrality_IDHmut_noncodel$surgical_interval~silver_neutrality_IDHmut_noncodel$recurrence_evolution)
fisher.test(table(silver_neutrality_IDHmut_noncodel$recurrence_evolution, silver_neutrality_IDHmut_noncodel$snv_driver_stability))
# Interesting, weak association between CNV driver stability and recurrent evolution.
fisher.test(table(silver_neutrality_IDHmut_noncodel$recurrence_evolution, silver_neutrality_IDHmut_noncodel$cnv_driver_stability))
fisher.test(table(silver_neutrality_IDHmut_noncodel$recurrence_evolution, silver_neutrality_IDHmut_noncodel$recurrence_location))
fisher.test(table(silver_neutrality_IDHmut_noncodel$recurrence_evolution, silver_neutrality_IDHmut_noncodel$received_tmz))
fisher.test(table(silver_neutrality_IDHmut_noncodel$recurrence_evolution, silver_neutrality_IDHmut_noncodel$tmz_std_care))
fisher.test(table(silver_neutrality_IDHmut_noncodel$recurrence_evolution, silver_neutrality_IDHmut_noncodel$received_rt))
fisher.test(table(silver_neutrality_IDHmut_noncodel$recurrence_evolution, silver_neutrality_IDHmut_noncodel$treatment))
fisher.test(table(silver_neutrality_IDHmut_noncodel$recurrence_evolution, silver_neutrality_IDHmut_noncodel$hypermutator_status))
fisher.test(table(silver_neutrality_IDHmut_noncodel$recurrence_evolution, silver_neutrality_IDHmut_noncodel$grade_change))


# IDH WT - neutrality is association with surgical interval.
wilcox.test(silver_neutrality_IDHwt$surgical_interval~silver_neutrality_IDHwt$recurrence_evolution)
fisher.test(table(silver_neutrality_IDHwt$recurrence_evolution, silver_neutrality_IDHwt$snv_driver_stability))
fisher.test(table(silver_neutrality_IDHwt$recurrence_evolution, silver_neutrality_IDHwt$cnv_driver_stability))
fisher.test(table(silver_neutrality_IDHwt$recurrence_evolution, silver_neutrality_IDHwt$recurrence_location))
fisher.test(table(silver_neutrality_IDHwt$recurrence_evolution, silver_neutrality_IDHwt$received_tmz))
fisher.test(table(silver_neutrality_IDHwt$recurrence_evolution, silver_neutrality_IDHwt$tmz_std_care))
fisher.test(table(silver_neutrality_IDHwt$recurrence_evolution, silver_neutrality_IDHwt$received_rt))
fisher.test(table(silver_neutrality_IDHwt$recurrence_evolution, silver_neutrality_IDHwt$treatment))
fisher.test(table(silver_neutrality_IDHwt$recurrence_evolution, silver_neutrality_IDHwt$hypermutator_status))
fisher.test(table(silver_neutrality_IDHwt$recurrence_evolution, silver_neutrality_IDHwt$grade_change))
wilcox.test(silver_neutrality_IDHwt$snv_driver_count~silver_neutrality_IDHwt$recurrence_evolution)
wilcox.test(silver_neutrality_IDHwt$snv_driver_count_private_b~silver_neutrality_IDHwt$recurrence_evolution)

## EVOLUTION MODE ##
# ALL SAMPLES.
# What's the average surgical interval? Longest for IDHwt N -> N.
silver_neutrality %>% 
  group_by(idh_codel_subtype, evolution_mode) %>% 
  summarise(med_int = median(surgical_interval, na.rm = T),
            sample_size = n())

# Since there are more than two groups, I now use the kruskal-wallis rank sum test. Overall, significantly associated with interval.
kruskal.test(silver_neutrality_IDHmut_codel$surgical_interval, as.factor(silver_neutrality_IDHmut_codel$evolution_mode))
kruskal.test(silver_neutrality_IDHmut_noncodel$surgical_interval, as.factor(silver_neutrality_IDHmut_noncodel$evolution_mode))
# Again, evo mode for IDHwt tumors is associated with surgical interval (not progression). Test whether this would be true for surgeries 1-2 only.
# For group N-N only
kruskal.test(silver_neutrality_IDHwt$surgical_interval, as.factor(silver_neutrality_IDHwt$evolution_mode))
kruskal.test(silver_neutrality_IDHwt$surgical_interval, as.factor(silver_neutrality_IDHwt$binary_mode))


# IDHwt. Largest group of samples, which is why we would expect an association over the others.
fisher.test(table(silver_neutrality_IDHwt$evolution_mode, silver_neutrality_IDHwt$snv_driver_stability))
fisher.test(table(silver_neutrality_IDHwt$evolution_mode, silver_neutrality_IDHwt$cnv_driver_stability))
fisher.test(table(silver_neutrality_IDHwt$evolution_mode, silver_neutrality_IDHwt$recurrence_location))
fisher.test(table(silver_neutrality_IDHwt$evolution_mode, silver_neutrality_IDHwt$received_tmz))
fisher.test(table(silver_neutrality_IDHwt$evolution_mode, silver_neutrality_IDHwt$tmz_std_care))
fisher.test(table(silver_neutrality_IDHwt$evolution_mode, silver_neutrality_IDHwt$received_rt))
fisher.test(table(silver_neutrality_IDHwt$evolution_mode, silver_neutrality_IDHwt$treatment))
fisher.test(table(silver_neutrality_IDHwt$evolution_mode, silver_neutrality_IDHwt$hypermutator_status))
fisher.test(table(silver_neutrality_IDHwt$evolution_mode, silver_neutrality_IDHwt$grade_change))
# Perform test with binary mode for Neutral-Neutral:
fisher.test(table(silver_neutrality_IDHwt$binary_mode, silver_neutrality_IDHwt$snv_driver_stability))
fisher.test(table(silver_neutrality_IDHwt$binary_mode, silver_neutrality_IDHwt$cnv_driver_stability))
fisher.test(table(silver_neutrality_IDHwt$binary_mode, silver_neutrality_IDHwt$recurrence_location))
fisher.test(table(silver_neutrality_IDHwt$binary_mode, silver_neutrality_IDHwt$received_tmz))
fisher.test(table(silver_neutrality_IDHwt$binary_mode, silver_neutrality_IDHwt$tmz_std_care))
fisher.test(table(silver_neutrality_IDHwt$binary_mode, silver_neutrality_IDHwt$received_rt))
fisher.test(table(silver_neutrality_IDHwt$binary_mode, silver_neutrality_IDHwt$hypermutator_status))
# There's very little grade change in GBM, it may not be relevant.
fisher.test(table(silver_neutrality_IDHwt$binary_mode, silver_neutrality_IDHwt$grade_change))
wilcox.test(silver_neutrality_IDHwt$snv_driver_count~silver_neutrality_IDHwt$binary_mode)
wilcox.test(silver_neutrality_IDHwt$snv_driver_count_private_b~silver_neutrality_IDHwt$binary_mode)

# Perform test with binary mode for Selection-Selection:
fisher.test(table(silver_neutrality_IDHwt$s_binary_mode, silver_neutrality_IDHwt$snv_driver_stability))
fisher.test(table(silver_neutrality_IDHwt$s_binary_mode, silver_neutrality_IDHwt$recurrence_location))
fisher.test(table(silver_neutrality_IDHwt$s_binary_mode, silver_neutrality_IDHwt$received_tmz))
fisher.test(table(silver_neutrality_IDHwt$s_binary_mode, silver_neutrality_IDHwt$tmz_std_care))
fisher.test(table(silver_neutrality_IDHwt$s_binary_mode, silver_neutrality_IDHwt$received_rt))
fisher.test(table(silver_neutrality_IDHwt$s_binary_mode, silver_neutrality_IDHwt$treatment))
fisher.test(table(silver_neutrality_IDHwt$s_binary_mode, silver_neutrality_IDHwt$hypermutator_status))
fisher.test(table(silver_neutrality_IDHwt$s_binary_mode, silver_neutrality_IDHwt$grade_change))
wilcox.test(silver_neutrality_IDHwt$snv_driver_count~silver_neutrality_IDHwt$s_binary_mode)
wilcox.test(as.numeric(silver_neutrality_IDHwt$snv_driver_count_private_b)~silver_neutrality_IDHwt$s_binary_mode)


## SURVIVAL ANALYSES ##
# By subtype: IDH MUT CODEL.
fit_codel_primary <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ primary_evolution,
                         data = silver_neutrality_IDHmut_codel)
ggsurvplot(fit_codel_primary, data = silver_neutrality_IDHmut_codel, risk.table = TRUE, pval= TRUE)
fit_codel_recur <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ recurrence_evolution,
                             data = silver_neutrality_IDHmut_codel)
ggsurvplot(fit_codel_recur, data = silver_neutrality_IDHmut_codel, risk.table = TRUE, pval= TRUE)

# IDH MUT NONCODEL
fit_noncodel_primary <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ primary_evolution,
                                data = silver_neutrality_IDHmut_noncodel)
ggsurvplot(fit_noncodel_primary, data = silver_neutrality_IDHmut_noncodel, risk.table = TRUE, pval= TRUE)
fit_noncodel_recur <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ recurrence_evolution,
                              data = silver_neutrality_IDHmut_noncodel)
ggsurvplot(fit_noncodel_recur, data = silver_neutrality_IDHmut_noncodel, risk.table = TRUE, pval= TRUE)

#### IDH WT.
# *** Neutrality at primary is associated with OS survival.
fit_wt_primary <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ primary_evolution,
                              data = silver_neutrality_IDHwt)
ggsurvplot(fit_wt_primary, data = silver_neutrality_IDHwt, risk.table = TRUE, pval= TRUE)
fit_wt_recur <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ recurrence_evolution,
                          data = silver_neutrality_IDHwt)
ggsurvplot(fit_wt_recur, data = silver_neutrality_IDHwt, risk.table = TRUE, pval= TRUE)
fit_wt_post_recur <- survfit(Surv(patient_post_recur_surv, patient_vital) ~ recurrence_evolution,
                        data = silver_neutrality_IDHwt)
ggsurvplot(fit_wt_post_recur, data = silver_neutrality_IDHwt, risk.table = TRUE, pval= TRUE)

# Mode of evolution:
# IDHmut CODEL.
fit_codel_mode <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ evolution_mode,
                           data = silver_neutrality_IDHmut_codel)
ggsurvplot(fit_codel_mode, data = silver_neutrality_IDHmut_codel, risk.table = TRUE, pval= TRUE)
fit_codel_mode_post <- survfit(Surv(patient_post_recur_surv, patient_vital) ~ evolution_mode,
                          data = silver_neutrality_IDHmut_codel)
ggsurvplot(fit_codel_mode_post, data = silver_neutrality_IDHmut_codel, risk.table = TRUE, pval= TRUE)
fit_codel_bi_mode <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ binary_mode,
                             data = silver_neutrality_IDHmut_codel)
ggsurvplot(fit_codel_bi_mode, data = silver_neutrality_IDHmut_codel, risk.table = TRUE, pval= TRUE)
fit_codel_binary_post <- survfit(Surv(patient_post_recur_surv, patient_vital) ~ binary_mode,
                               data = silver_neutrality_IDHmut_codel)
ggsurvplot(fit_codel_binary_post, data = silver_neutrality_IDHmut_codel, risk.table = TRUE, pval= TRUE)

# IDHmut NONCODEL.
fit_noncodel_mode <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ evolution_mode,
                              data = silver_neutrality_IDHmut_noncodel)
ggsurvplot(fit_noncodel_mode, data = silver_neutrality_IDHmut_noncodel, risk.table = TRUE, pval= TRUE)
fit_noncodel_mode_post <- survfit(Surv(patient_post_recur_surv, patient_vital) ~ evolution_mode,
                             data = silver_neutrality_IDHmut_noncodel)
ggsurvplot(fit_noncodel_mode_post, data = silver_neutrality_IDHmut_noncodel, risk.table = TRUE, pval= TRUE)
fit_noncodel_bi_mode <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ binary_mode,
                                data = silver_neutrality_IDHmut_noncodel)
ggsurvplot(fit_noncodel_bi_mode, data = silver_neutrality_IDHmut_noncodel, risk.table = TRUE, pval= TRUE)

# IDHwt.
fit_wt_mode <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ evolution_mode,
                        data = silver_neutrality_IDHwt)
ggsurvplot(fit_wt_mode, data = silver_neutrality_IDHwt, risk.table = TRUE, pval= TRUE)
fit_wt_mode_post <- survfit(Surv(patient_post_recur_surv, patient_vital) ~ evolution_mode,
                       data = silver_neutrality_IDHwt)
ggsurvplot(fit_wt_mode_post, data = silver_neutrality_IDHwt, risk.table = TRUE, pval= TRUE)

fit_wt_bi_mode <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ binary_mode,
                       data = silver_neutrality_IDHwt)
ggsurvplot(fit_wt_bi_mode, data = silver_neutrality_IDHwt, risk.table = TRUE, pval= TRUE)
fit_wt_bi_mode_post <- survfit(Surv(patient_post_recur_surv, patient_vital) ~ binary_mode,
                          data = silver_neutrality_IDHwt)
ggsurvplot(fit_wt_bi_mode_post, data = silver_neutrality_IDHwt, risk.table = TRUE, pval= TRUE)
# **
fit_wt_bi_mode <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ s_binary_mode,
                          data = silver_neutrality_IDHwt)
ggsurvplot(fit_wt_bi_mode, data = silver_neutrality_IDHwt, risk.table = TRUE, pval= TRUE)
fit_wt_bi_mode <- survfit(Surv(patient_post_recur_surv, patient_vital) ~ s_binary_mode,
                          data = silver_neutrality_IDHwt)
ggsurvplot(fit_wt_bi_mode, data = silver_neutrality_IDHwt, risk.table = TRUE, pval= TRUE)

## Filter to only have two-groups:
silver_neutrality_IDHwt_2groups = silver_neutrality_IDHwt %>% 
  filter(evolution_mode%in%c("N-N", "S-S"))
fit_wt_bi_mode <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ evolution_mode,
                          data = silver_neutrality_IDHwt_2groups)
ggsurvplot(fit_wt_bi_mode, data = silver_neutrality_IDHwt_2groups, risk.table = TRUE, pval= TRUE)
fit_wt_bi_mode_post <- survfit(Surv(patient_post_recur_surv, patient_vital) ~ evolution_mode,
                          data = silver_neutrality_IDHwt_2groups)
ggsurvplot(fit_wt_bi_mode_post, data = silver_neutrality_IDHwt_2groups, risk.table = TRUE, pval= TRUE)


####### Cox proportional hazards model #####
res_cox_idh_wt <- coxph(Surv(case_overall_survival_mo, patient_vital)~ case_age_diagnosis_years + primary_evolution, data = silver_neutrality_IDHwt)
summary(res_cox_idh_wt)
res_cox_idh_wt <- coxph(Surv(case_overall_survival_mo, patient_vital)~ case_age_diagnosis_years + recurrence_evolution, data = silver_neutrality_IDHwt)
summary(res_cox_idh_wt)
res_cox_idh_wt <- coxph(Surv(case_overall_survival_mo, patient_vital)~ case_age_diagnosis_years + binary_mode, data = silver_neutrality_IDHwt)
summary(res_cox_idh_wt)
res_cox_idh_wt <- coxph(Surv(case_overall_survival_mo, patient_vital)~ case_age_diagnosis_years + binary_mode, data = silver_neutrality_IDHwt)
summary(res_cox_idh_wt)
# PRIMARY evolution is not significant when included primary evolution as a predictor in the cox model.
res_cox_drivers <- coxph(Surv(case_overall_survival_mo, patient_vital)~ case_age_diagnosis_years + idh_codel_subtype + primary_evolution, data = silver_neutrality)
summary(res_cox_drivers)


#### PRIMARY WITH SURGICAL INTERVAL.
fit_codel_interval <- survfit(Surv(surgical_interval, rep(1, length(silver_neutrality_IDHmut_codel$surgical_interval))) ~ primary_evolution,
                              data = silver_neutrality_IDHmut_codel)
ggsurvplot(fit_codel_interval, data = silver_neutrality_IDHmut_codel, risk.table = TRUE, pval= TRUE, pval.method = TRUE, ylab = "Proportion Progression Free")
fit_noncodel_interval <- survfit(Surv(surgical_interval, rep(1, length(silver_neutrality_IDHmut_noncodel$surgical_interval))) ~ primary_evolution,
                                 data = silver_neutrality_IDHmut_noncodel)
ggsurvplot(fit_noncodel_interval, data = silver_neutrality_IDHmut_noncodel, risk.table = TRUE, pval= TRUE, pval.method = TRUE, ylab = "Proportion Progression Free")
fit_IDHwt_interval <- survfit(Surv(surgical_interval, rep(1, length(silver_neutrality_IDHwt$surgical_interval))) ~ primary_evolution,
                              data = silver_neutrality_IDHwt)
ggsurvplot(fit_IDHwt_interval, data = silver_neutrality_IDHwt, risk.table = TRUE, pval= TRUE, pval.method = TRUE, ylab = "Proportion Progression Free")


#### RECURRENCE WITH SURGICAL INTERVAL.
fit_codel_interval <- survfit(Surv(surgical_interval, rep(1, length(silver_neutrality_IDHmut_codel$surgical_interval))) ~ recurrence_evolution,
                              data = silver_neutrality_IDHmut_codel)
ggsurvplot(fit_codel_interval, data = silver_neutrality_IDHmut_codel, risk.table = TRUE, pval= TRUE, pval.method = TRUE, ylab = "Proportion Progression Free")
fit_noncodel_interval <- survfit(Surv(surgical_interval, rep(1, length(silver_neutrality_IDHmut_noncodel$surgical_interval))) ~ recurrence_evolution,
                                 data = silver_neutrality_IDHmut_noncodel)
ggsurvplot(fit_noncodel_interval, data = silver_neutrality_IDHmut_noncodel, risk.table = TRUE, pval= TRUE, pval.method = TRUE, ylab = "Proportion Progression Free")
#### ***Neutrality at recurrence was associated with a longer surgical interval.
fit_IDHwt_interval <- survfit(Surv(surgical_interval, rep(1, length(silver_neutrality_IDHwt$surgical_interval))) ~ recurrence_evolution,
                              data = silver_neutrality_IDHwt)
ggsurvplot(fit_IDHwt_interval, data = silver_neutrality_IDHwt, risk.table = TRUE, pval= TRUE, pval.method = TRUE, ylab = "Proportion Progression Free")

# EVOLUTION MODE WITH SURGICAL INTERVAL.
fit_IDHwt_interval <- survfit(Surv(surgical_interval, rep(1, length(silver_neutrality_IDHwt$surgical_interval))) ~ binary_mode,
                              data = silver_neutrality_IDHwt)
ggsurvplot(fit_IDHwt_interval, data = silver_neutrality_IDHwt, risk.table = TRUE, pval= TRUE, pval.method = TRUE, ylab = "Proportion Progression Free")

fit_IDHwt_interval <- survfit(Surv(surgical_interval, rep(1, length(silver_neutrality_IDHwt$surgical_interval))) ~ s_binary_mode,
                              data = silver_neutrality_IDHwt)
ggsurvplot(fit_IDHwt_interval, data = silver_neutrality_IDHwt, risk.table = TRUE, pval= TRUE, pval.method = TRUE, ylab = "Proportion Progression Free")

################################################
# Neutrality stacked bar plots visualizations  #
################################################
stacked_driver = silver_neutrality %>% 
  group_by(idh_codel_subtype, primary_evolution) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n))
ggplot(stacked_driver, aes(x=idh_codel_subtype, y=freq, fill=primary_evolution)) + geom_bar(stat="identity") +
  labs(fill="Evolution status at primary") + xlab("") + theme_bw()

stacked_neutrality_recur = silver_neutrality %>% 
  group_by(idh_codel_subtype, recurrence_evolution) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n))
ggplot(stacked_neutrality_recur, aes(x=idh_codel_subtype, y=freq, fill=recurrence_evolution)) + geom_bar(stat="identity") +
  labs(fill="Evolution status at recurrence") + xlab("") + theme_bw() + ylab("Proportion of silver set pairs")

stacked_neutrality = silver_neutrality %>% 
  group_by(idh_codel_subtype, evolution_mode) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n))
ggplot(stacked_neutrality, aes(x=idh_codel_subtype, y=freq, fill=evolution_mode)) + geom_bar(stat="identity") +
  labs(fill="Evolution route") + xlab("") + theme_bw() + ylab("Proportion of silver set pairs")


#################################
# Compare R^2 values for IDHwt samples.
################################
# No substantial difference in paired R-squared values.
wilcox.test(silver_neutrality_IDHwt$model_rsq.x, silver_neutrality_IDHwt$model_rsq.y, paired = T)

# Quickly compare the R-squared values.
plot_neutrality_rsq = silver_neutrality_IDHwt %>% 
  gather("rsq_time", "rsq_value", c(model_rsq.x, model_rsq.y), -tumor_pair_barcode)
ggplot(plot_neutrality_rsq, aes(x = rsq_time, y = rsq_value, group = tumor_pair_barcode)) +
  geom_line(linetype="solid", size=1, alpha = 0.3) + 
  geom_point(color="black", size=2) + theme_bw() + ylab("R-squared values") + xlab("") + geom_hline(yintercept = 0.98, linetype="dotted") +
  geom_text(x = 1, y = 0.7, label="P=0.51") + ggtitle("Mutect2, primary_all, recurrence_all (n=61 selected pairs)")

