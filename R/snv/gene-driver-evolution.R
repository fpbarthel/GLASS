#######################################################
# Analyse gene drivers over time.
# Date: 2019.02.02 
# Author: Kevin J.
#######################################################

# Necessary packages:
library(tidyverse)
library(DBI)
library(gridExtra)
library(nlme)
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
clinal_tumor_pairs = dbReadTable(con,  Id(schema="clinical", table="clinical_by_tumor_pair"))  

# Load local table for "neutrality"
subclonal_neutrality = read_tsv("/Users/johnsk/Documents/Life-History/glass-analyses/data/glass-evolution-neutrality-20190129.txt")
subclonal_neutrality = subclonal_neutrality %>% 
  select(-idh_codel_subtype)

# These tables **MAY** change, especially the driver table.
silver_set = dbGetQuery(con, "SELECT * FROM analysis.silver_set")
all_drivers = dbGetQuery(con, "SELECT * FROM analysis.driver_status")

# Define silver set drivers to perform association tests with neutrality estimates.
silver_drivers = silver_set %>% 
  inner_join(all_drivers, by="tumor_pair_barcode") %>% 
  inner_join(clinal_tumor_pairs, by="tumor_pair_barcode") %>% 
  left_join(subclonal_neutrality, by="tumor_pair_barcode") %>% 
  inner_join(cases, by="case_barcode") 
silver_drivers$any_driver_stability = ifelse(silver_drivers$snv_driver_stability=="Driver unstable" | silver_drivers$cnv_driver_stability=="Driver unstable", "Driver unstable", "Driver stable") 
silver_drivers$bay_evo_mode = paste(silver_drivers$bayesian_evo_p, silver_drivers$bayesian_evo_r, sep="-") 
silver_drivers$bay_evo_mode = ifelse(silver_drivers$bay_evo_mode=="NA-NA", NA, silver_drivers$bay_evo_mode)

# Combine with pairs file to access "tumor_barcode".
titan_info = titan_param %>% 
  left_join(pairs, by="pair_barcode")

# There is a significant difference in the number of drivers in the tumor pair. This is probably driven by IDH being counted as driver (every codel/noncodel has +1)
kruskal.test(as.numeric(silver_drivers$snv_driver_count), as.factor(silver_drivers$idh_codel_subtype))
ggplot(silver_drivers, aes(x=as.factor(idh_codel_subtype), y=as.numeric(snv_driver_count))) + geom_boxplot() + theme_bw() +
  ylab("Total driver count") + xlab("")
  
# Examine whether there are any particular gains or losses in the non-shared IDHwt tumors.
no_shared_drivers = silver_drivers %>% 
  filter(is.na(snv_driver_shared)) 

# Bar plot tallying the total number of events shared. SNV drivers, CNV drivers, Chr7/10.
## SNV drivers.
silver_drivers_snv_sum = silver_drivers %>% 
  mutate_if(bit64::is.integer64, as.double) %>% 
  group_by(idh_codel_subtype) %>% 
  summarise(total_shared = sum(snv_driver_count_shared, na.rm = T),
            total_private_a = sum(snv_driver_count_private_a, na.rm = T),
            total_private_b = sum(snv_driver_count_private_b, na.rm = T)) %>% 
  gather(driver_type, frequency, c(total_shared, total_private_a, total_private_b)) %>% 
  mutate(type = factor(driver_type,
                       levels = c("total_shared", "total_private_a", "total_private_b"),
                       labels = c("Shared", "Primary", "Recurrence"))) %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, "IDHmut_codel" = "IDHmut codel", "IDHmut_noncodel" = "IDHmut noncodel", "IDHwt_noncodel" = "IDHwt")) 
ggplot(silver_drivers_snv_sum, aes(x=idh_codel_subtype, y=frequency, fill=type)) + geom_bar(position="dodge", stat="identity") + 
  scale_fill_manual(values=c("Shared"="#CA932F", "Primary" ="#CA2F66", "Recurrence" = "#2FB3CA"), name = "sSNV Driver Type") +
  theme_bw() + xlab("") + ylab("Frequency")
  
## CNV shared vs. private.
silver_drivers_cnv_sum = silver_drivers %>% 
  mutate_if(bit64::is.integer64, as.double) %>% 
  group_by(idh_codel_subtype) %>% 
  summarise(total_shared = sum(cnv_driver_count_shared, na.rm = T),
            total_private_a = sum(cnv_driver_count_private_a, na.rm = T),
            total_private_b = sum(cnv_driver_count_private_b, na.rm = T)) %>% 
  gather(driver_type, frequency, c(total_shared, total_private_a, total_private_b)) %>% 
  mutate(type = factor(driver_type,
                       levels = c("total_shared", "total_private_a", "total_private_b"),
                       labels = c("Shared", "Primary", "Recurrence"))) %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, "IDHmut_codel" = "IDHmut codel", "IDHmut_noncodel" = "IDHmut noncodel", "IDHwt_noncodel" = "IDHwt")) 
ggplot(silver_drivers_cnv_sum, aes(x=idh_codel_subtype, y=frequency, fill=type)) + geom_bar(position="dodge", stat="identity") + 
  scale_fill_manual(values=c("Shared"="#CA932F", "Primary" ="#CA2F66", "Recurrence" = "#2FB3CA"), name = "CNV Driver Type") +
  theme_bw() + xlab("") + ylab("Frequency")


# Chr7/Chr10 preservation.
c710data <- dbGetQuery(con, read_file("sql/heatmap/heatmap_c710.sql"))
c710data_revised = c710data %>% 
  mutate(c710_status = ifelse(is.na(c710_status), "NA", c710_status)) %>% 
  mutate(c710_status = factor(c710_status,
                                    levels = c("S", "P", "R", "WT", "NA"),
                                    labels = c("Shared", "Primary", "Recurrence", "WT", "NA"))) %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, "IDHmut_codel" = "IDHmut codel", "IDHmut_noncodel" = "IDHmut noncodel", "IDHwt_noncodel" = "IDHwt")) %>% 
  group_by(idh_codel_subtype, c710_status) %>% 
  summarise(frequency = n()) 

# Supplemental Fig.
ggplot(c710data_revised, aes(x=idh_codel_subtype, y=frequency, fill=c710_status)) + geom_bar(position="dodge", stat="identity") + 
  scale_fill_manual(values=c("Shared"="#CA932F", "Primary" ="#CA2F66", "Recurrence" = "#2FB3CA", "WT" = "darkgray", "NA" = "lightgray"), name = "Chr7+/Chr10- \n Stability") +
  theme_bw() + xlab("") + ylab("Frequency")

# Define subtype categories:
silver_drivers = silver_drivers %>% 
  mutate(snv_driver_context_change_revised = ifelse(is.na(snv_driver_context_change), "SNV-stable", snv_driver_context_change))
# The status indicator, normally is 0=alive, 1=dead.
silver_drivers$patient_vital = ifelse(silver_drivers$case_vital_status=="alive", 0, 1)
# Create a variable for post-recurrence survival (also presented in months).
silver_drivers$patient_post_recur_surv = silver_drivers$case_overall_survival_mo-silver_drivers$surgical_interval
silver_drivers$binary_mode = ifelse(silver_drivers$bay_evo_mode=="N-N", "Neutral-Neutral", "N-S|S-N|S-S")
silver_drivers$s_binary_mode = ifelse(silver_drivers$bay_evo_mode=="S-S", "Selection-Selection", "N-S|S-N|N-N")

#############################
# FILTERED out hypermutants
#############################
silver_drivers = silver_drivers %>% 
  filter(hypermutator_status == 0)

silver_drivers_codel = silver_drivers %>% 
  filter(idh_codel_subtype == "IDHmut_codel")
silver_drivers_noncodel = silver_drivers %>% 
  filter(idh_codel_subtype == "IDHmut_noncodel")
silver_drivers_IDHwt = silver_drivers %>% 
  filter(idh_codel_subtype == "IDHwt_noncodel")

# How to incorporate in-context/out-of-context? SNVs first
fisher.test(table(silver_drivers$snv_driver_stability, silver_drivers$grade_change))
# *** Interesting association with grade change.
fisher.test(table(silver_drivers$snv_driver_context_change, silver_drivers$grade_change))
# Weak association betwen snv driver stability and recurrence location. 
# *** Distal recurrence.
fisher.test(table(silver_drivers$snv_driver_stability, silver_drivers$recurrence_location))
fisher.test(table(silver_drivers$snv_driver_context_change, silver_drivers$recurrence_location))
fisher.test(table(silver_drivers$snv_driver_stability, silver_drivers$received_alkylating_agent))
fisher.test(table(silver_drivers$snv_driver_stability, silver_drivers$received_rt))
fisher.test(table(silver_drivers$snv_driver_stability, silver_drivers$bayesian_evo_r))
wilcox.test(silver_drivers$received_rt_sum_gy~silver_drivers$snv_driver_stability)
wilcox.test(silver_drivers$received_tmz_sum_cycles~silver_drivers$snv_driver_stability)

# CNV 
fisher.test(table(silver_drivers$cnv_driver_stability, silver_drivers$grade_change))
#  Out-of-context copy number change also associated with grade change.
fisher.test(table(silver_drivers$cnv_driver_context_change, silver_drivers$grade_change))
fisher.test(table(silver_drivers$cnv_driver_stability, silver_drivers$recurrence_location))
fisher.test(table(silver_drivers$cnv_driver_context_change, silver_drivers$recurrence_location))
fisher.test(table(silver_drivers$cnv_driver_stability, silver_drivers$received_alkylating_agent))
fisher.test(table(silver_drivers$cnv_driver_stability, silver_drivers$received_tmz))
# Radiotherapy was associated with cnv driver instability.
fisher.test(table(silver_drivers$cnv_driver_stability, silver_drivers$received_rt))
fisher.test(table(silver_drivers$cnv_driver_stability, silver_drivers$bayesian_evo_r))
wilcox.test(silver_drivers$received_rt_sum_gy~silver_drivers$cnv_driver_stability)
wilcox.test(silver_drivers$received_tmz_sum_cycles~silver_drivers$cnv_driver_stability)


# Are there any subtype effects for the two variables that were associated with therapy/grade change
fisher.test(table(silver_drivers$snv_driver_context_change_revised, silver_drivers$grade_change))
fisher.test(table(silver_drivers_codel$snv_driver_context_change_revised, silver_drivers_codel$grade_change))
fisher.test(table(silver_drivers_noncodel$snv_driver_context_change_revised, silver_drivers_noncodel$grade_change))
fisher.test(table(silver_drivers_IDHwt$snv_driver_context_change_revised, silver_drivers_IDHwt$grade_change))

# Does it associate with evolutionary mode?
fisher.test(table(silver_drivers$snv_driver_stability, silver_drivers$bayesian_evo_r))
fisher.test(table(silver_drivers$snv_driver_stability, silver_drivers$bay_evo_mode))
fisher.test(table(silver_drivers$snv_driver_stability, silver_drivers$binary_mode))
fisher.test(table(silver_drivers$snv_driver_stability, silver_drivers$s_binary_mode))

# RT
fisher.test(table(silver_drivers$cnv_driver_stability, silver_drivers$received_rt))
fisher.test(table(silver_drivers_codel$cnv_driver_stability, silver_drivers_codel$received_rt))
# Main change is associated IDHmut noncodels.
fisher.test(table(silver_drivers_noncodel$cnv_driver_stability, silver_drivers_noncodel$received_rt))
fisher.test(table(silver_drivers_IDHwt$cnv_driver_stability, silver_drivers_IDHwt$received_rt))



###
# What about associations with progression, OS, and post-recurrence survival?
###
# No difference in surgical interval.
wilcox.test(silver_drivers_IDHwt$surgical_interval~factor(silver_drivers_IDHwt$snv_driver_stability))
wilcox.test(silver_drivers_noncodel$surgical_interval~factor(silver_drivers_noncodel$snv_driver_stability))
wilcox.test(silver_drivers_codel$surgical_interval~factor(silver_drivers_codel$snv_driver_stability))

# *** Weak association between SNV driver instability and shorter survival.
res_cox_driver_snv <- coxph(Surv(case_overall_survival_mo, patient_vital)~ case_age_diagnosis_years + idh_codel_subtype + snv_driver_stability, data = silver_drivers)
summary(res_cox_driver_snv)
res_cox_driver_cnv <- coxph(Surv(case_overall_survival_mo, patient_vital)~ case_age_diagnosis_years + idh_codel_subtype + cnv_driver_stability, data = silver_drivers)
summary(res_cox_driver_cnv)
res_cox_driver_any <- coxph(Surv(case_overall_survival_mo, patient_vital)~ case_age_diagnosis_years + idh_codel_subtype + any_driver_stability, data = silver_drivers)
summary(res_cox_driver_any)

# No associations with post-recurrence surv.
res_cox_driver_snv <- coxph(Surv(patient_post_recur_surv, patient_vital)~ case_age_diagnosis_years + idh_codel_subtype + snv_driver_stability, data = silver_drivers)
summary(res_cox_driver_snv)
res_cox_driver_cnv <- coxph(Surv(patient_post_recur_surv, patient_vital)~ case_age_diagnosis_years + idh_codel_subtype + cnv_driver_stability, data = silver_drivers)
summary(res_cox_driver_cnv)
res_cox_driver_any <- coxph(Surv(patient_post_recur_surv, patient_vital)~ case_age_diagnosis_years + idh_codel_subtype + any_driver_stability, data = silver_drivers)
summary(res_cox_driver_any)


