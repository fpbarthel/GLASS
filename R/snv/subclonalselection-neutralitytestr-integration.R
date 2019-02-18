##################################################
# SubClonalSelection and neutralitytestR results integration
# Updated: 2019.01.30
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

##################################################
# Establish connection with database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

# Load additional tables from the database.
tumor_blocklist = dbReadTable(con,  Id(schema="analysis", table="blocklist"))
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

# Use this information in the context of the silver set for aneuploidy pairs. Note: that this includes some review CN data.
aneuploidy_pairs = dbGetQuery(con, "SELECT tumor_pair_barcode, tumor_barcode_a, tumor_barcode_b, a1.aneuploidy AS aneuploidy_a, a2.aneuploidy AS aneuploidy_b, t1.aneuploidy_score::integer AS aneuploidy_score_a,t2.aneuploidy_score::integer AS aneuploidy_score_b,idh_codel_subtype
                              FROM analysis.silver_set ss
                              LEFT JOIN analysis.aneuploidy a1 ON a1.aliquot_barcode = ss.tumor_barcode_a
                              LEFT JOIN analysis.aneuploidy a2 ON a2.aliquot_barcode = ss.tumor_barcode_b
                              LEFT JOIN analysis.taylor_aneuploidy t1 ON t1.aliquot_barcode = ss.tumor_barcode_a
                              LEFT JOIN analysis.taylor_aneuploidy t2 ON t2.aliquot_barcode = ss.tumor_barcode_b
                              LEFT JOIN biospecimen.aliquots al ON al.aliquot_barcode = ss.tumor_barcode_a
                              LEFT JOIN clinical.surgeries su ON su.sample_barcode = al.sample_barcode")


# Filter out any poor performing samples on CNV blocklist. 197 patients because I merged silver set with allow|review CN.
# Gold set has two tumor pairs (different combination of tumors) that are not in silver set.
aneuploidy_pairs_filtered = aneuploidy_pairs %>% 
  inner_join(tumor_blocklist, by=c("tumor_barcode_a"="aliquot_barcode")) %>% 
  inner_join(tumor_blocklist, by=c("tumor_barcode_b"="aliquot_barcode")) %>% 
  # Note db introduced some trailing whitespace.
  filter(cnv_exclusion.x %in%c("allow ", "review")) %>% 
  filter(cnv_exclusion.y %in% c("allow ", "review")) %>% 
  mutate(aneuploidy_diff = aneuploidy_b - aneuploidy_a) %>% 
  left_join(clinal_tumor_pairs, by="tumor_pair_barcode") %>%
  # Relabel subtypes for aesthetics.
  mutate(idh_codel_subtype = recode(idh_codel_subtype, "IDHmut_codel" = "IDHmut codel", "IDHmut_noncodel" = "IDHmut noncodel", "IDHwt_noncodel" = "IDHwt")) %>% 
  select(tumor_pair_barcode, aneuploidy_a:aneuploidy_score_b, aneuploidy_diff)


# Although all aliquot results are present here, BUT need to filter low subclonal mutation count [purity < 0.4 AND ploidy > 3].
# Combine the neutrality at the aliquot with specific tumor_pairs.
per_sample_neutrality = tumor_pairs %>% 
  left_join(aliquot_neutrality, by=c("tumor_barcode_a"="aliquot_barcode")) %>% 
  left_join(aliquot_neutrality, by=c("tumor_barcode_b"="aliquot_barcode")) %>% 
  filter(subclonal_mut.x >= 12, subclonal_mut.y >= 12, cellularity.x > 0.4, cellularity.y > 0.4, ploidy.x < 3, ploidy.y < 3) %>% 
  mutate(classification_p = ifelse(area_pval.x < 0.05, "S", "N"), 
         classification_r = ifelse(area_pval.y < 0.05, "S", "N"),
         evolution_mode = paste(classification_p, classification_r, sep="-")) %>% 
  left_join(all_drivers, by = "tumor_pair_barcode") %>% 
  left_join(clinal_tumor_pairs, by = "tumor_pair_barcode") 

# Create a standard. of care treatment group for TMZ (at least 6 cycles).
per_sample_neutrality$tmz_std_care <- ifelse(per_sample_neutrality$received_tmz_sum_cycles >= 6, "1", "0")

# Subset to silver set so that a patient is only represented once at a single time point. N = 98 samples.
silver_all_evolution = per_sample_neutrality %>% 
  inner_join(silver_set, by="tumor_pair_barcode") %>% 
  left_join(cases, by=c("case_barcode.x"="case_barcode")) %>% 
  mutate(treatment = ifelse(received_tmz == 1 | received_rt == 1, "YES", "NO"),
         treatment = ifelse(is.na(treatment), "Unknown", treatment))


# Georgette provided SubClonalSelection results on 2019.01.17.
# subclonalselect_probs = read_tsv("/Users/johnsk/Documents/Life-History/glass-analyses/data/subclonalselection_formatted_probabilities.txt")
subclonalselect_probs = read_tsv("/Users/johnsk/Documents/Life-History/glass-analyses/data/SubClonalSelection/2019.01.29-results/formatted_metadata3_probabilities.txt")

# Restrict analysis to the aliquot portion.
subclonal_aliquots = subclonalselect_probs %>% 
  separate(sample, c("aliquot_barcode", "mutation_type"), sep ="_") %>% 
  filter(mutation_type == "all") %>% 
  mutate(p0 = `0subclones`,
         p3 = `1subclones` + `2subclones`)

# 17 samples that were majority, but didn't reach 0.6 for selection.
sum(subclonal_aliquots$classification=="S")

silver_all_bayesian_evolution = silver_set %>% 
  left_join(subclonal_aliquots, by=c("tumor_barcode_a" = "aliquot_barcode")) %>% 
  left_join(subclonal_aliquots, by=c("tumor_barcode_b" = "aliquot_barcode")) %>% 
  filter(!is.na(`0subclones.x`), !is.na(`0subclones.y`))

# Visualize the distributions for probability of neutrality and 
ggplot(subclonal_aliquots, aes(x = p0)) + geom_histogram(binwidth=0.01) + geom_vline(xintercept = 0.6, color="red") +theme_bw() +
  ggtitle("SubClonalSelection Probabilities - Neutral") 
ggplot(subclonal_aliquots, aes(x = p3)) + geom_histogram(binwidth=0.01) + geom_vline(xintercept = 0.6, color="red") +theme_bw() +
  ggtitle("SubClonalSelection Probabilities - Selection") 

# Combine the two datasets (minimum 12 subclonal mutations). 
neutrality_overlap = subclonal_aliquots %>% 
  inner_join(aliquot_neutrality, by="aliquot_barcode") %>% 
  select(-evolution) %>% 
  mutate(evolution_area_pval = ifelse(area_pval < 0.05, "S", "N"))
  
# How do they compare with one another?
# 30 subclonal mutants:
neutrality_overlap %>% 
  filter(subclonal_mut > 30, cellularity > 0.4, ploidy < 3) %>% 
  group_by(classification, evolution_area_pval) %>% 
  summarise(evo_proportions = n())
# Relax criteria to 12 subclonal mutants:
neutrality_overlap %>% 
  filter(subclonal_mut > 12, cellularity > 0.4, ploidy < 3) %>% 
  group_by(classification, evolution_area_pval) %>% 
  summarise(evo_proportions = n())

# Redefine "neutrality" for those "uncertain" samples.
neutrality_overlap2 = neutrality_overlap %>% 
  select(-classification) %>% 
  mutate(classification = ifelse(p0 > p3, "N", "S"))

# Bring the samples together to see how estimates impact evolution.
silver_all_evolution = silver_set %>% 
  left_join(neutrality_overlap2, by=c("tumor_barcode_a" = "aliquot_barcode")) %>% 
  left_join(neutrality_overlap2, by=c("tumor_barcode_b" = "aliquot_barcode")) %>% 
  left_join(all_drivers, by="tumor_pair_barcode") %>% 
  mutate(bay_evo_mode = paste(classification.x, classification.y, sep="-"),
         freq_evo_mode = paste(evolution_area_pval.x, evolution_area_pval.y, sep="-")) %>% 
  filter(!is.na(classification.x), !is.na(classification.y)) %>% 
  left_join(clinal_tumor_pairs, by="tumor_pair_barcode") %>% 
  left_join(cases, by="case_barcode") %>% 
  mutate(treatment = ifelse(received_tmz == 1 | received_rt == 1, "YES", "NO"),
         treatment = ifelse(is.na(treatment), "Unknown", treatment),
         extensive_treatment = ifelse(received_tmz == 1 & received_rt == 1, "Concomitant", ifelse(treatment=="YES", "Single-agent", treatment))) #%>% 
#  filter(!grepl("U", bay_evo_mode))

colnames(silver_all_evolution) = gsub("\\.x", "_p", colnames(silver_all_evolution))
colnames(silver_all_evolution) = gsub("\\.y", "_r", colnames(silver_all_evolution))


# Chi-square test whether there is a significant difference in the estimates based on proportions.
evo_05_06p = matrix(c(95,89,32,23), nrow = 2)
chisq.test(evo_05_06p)
evo_05_06r = matrix(c(89,78,38,34), nrow = 2)
chisq.test(evo_05_06r)

evo_05_07p = matrix(c(95,75,32,23), nrow = 2)
chisq.test(evo_05_07p)
evo_05_07r = matrix(c(89,63,38,34), nrow = 2)
chisq.test(evo_05_07r)

# How well do the primary samples match up?
silver_all_evolution %>% 
  group_by(classification_p, evolution_area_pval_p) %>% 
  summarise(sample_breakdown =  n())
# Don't count `U` samples.
(78+15)/(78+15+20+8) # 75% concordance for primary tumors
# How well do the recurrent samples match up?
silver_all_evolution %>% 
  group_by(classification_r, evolution_area_pval_r) %>% 
  summarise(sample_breakdown =  n())
# Don't count `U` samples.
(60+20)/(60+20+26+16) # 66% concordance for recurrent tumors.

# What about evolution mode?
silver_all_evolution %>% 
  group_by(bay_evo_mode, freq_evo_mode) %>% 
  summarise(sample_breakdown =  n())

############################
# Test for associations with clinical covariates
############################
# First snakey plot:
## Color neutral-neutral or selection-selection.
# Use subtype, driver stability, and evolution pattern.
neutral_dat = silver_all_evolution %>% 
  select(tumor_pair_barcode, idh_codel_subtype, hypermutator_status, received_tmz, received_rt, treatment, classification_p, classification_r) %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, "IDHmut_codel" = "IDHmut codel", "IDHmut_noncodel" = "IDHmut noncodel", "IDHwt_noncodel" = "IDHwt"),
         classification_p = recode(classification_p, "N" = "Neutral", "S" = "Selected", "U" = "Unknown"),
         classification_r = recode(classification_r, "N" = "Neutral", "S" = "Selected", "U" = "Unknown"))  %>% 
  group_by(idh_codel_subtype, classification_p, treatment, classification_r) %>% 
  summarise(Freq = n())

# Adjust the alluvial plot so that different colors map to different selection strength.
# Color scheme for evolution trajectory.
#2FB3CA - IDHwt
#F1564F - IDHmut noncodel
#F69654 - IDHmut codel

# Default colors
pal = c("#619CFF", "#00BA38", "#F8766D", "lightgray")
neutral_dat$col_type = ifelse(neutral_dat$idh_codel_subtype == "IDHwt", pal[1], ifelse(neutral_dat$idh_codel_subtype == "IDHmut noncodel", pal[2], pal[3]))
neutral_dat$col_type = ifelse(neutral_dat$treatment == "Unknown", pal[4], neutral_dat$col_type)
neutral_dat$col_type = ifelse(neutral_dat$classification_p == "U"|neutral_dat$classification_r == "U", pal[4], neutral_dat$col_type)

# The alluvial algorithm doesn't want `NAs` as input. Revise about to use "Unknown".
pdf(file = "/Users/johnsk/Documents/sankey-plot-evolutionary-modes.pdf", height = 5, width = 10, bg = "transparent", useDingbats = FALSE)
alluvial(neutral_dat[,1:4], freq=neutral_dat$Freq, border=NA,
         #       hide = neutral_dat$Freq < quantile(neutral_dat$Freq, .50),
         col= neutral_dat$col_type, alpha = 0.5,
         cw = 0.25, axis_labels = c("Glioma subtype", "Primary tumor", "Treatment", "Recurrent tumor"))
dev.off()
# Create legend for alluvial plot.
plot.new()
legend("center", title="Glioma subtype",
       c("IDHwt","IDHmut noncodel","IDHmut codel", "Treatment unknown"), fill=pal, cex=1.5, bty = "n")


## Create revised survival variables for Kaplan-Meier curves and log-rank tests ##
# The status indicator, normally is 0=alive, 1=dead.
silver_all_evolution$patient_vital = ifelse(silver_all_evolution$case_vital_status=="alive", 0, 1)
# Create a variable for post-recurrence survival (also presented in months).
silver_all_evolution$patient_post_recur_surv = silver_all_evolution$case_overall_survival_mo-silver_all_evolution$surgical_interval
# Treat all "neutrality" and all "selection" as groups by themselves
silver_all_evolution$binary_mode = ifelse(silver_all_evolution$bay_evo_mode=="N-N", "Neutral-Neutral", "N-S|S-N|S-S")
silver_all_evolution$s_binary_mode = ifelse(silver_all_evolution$bay_evo_mode=="S-S", "Selection-Selection", "N-S|S-N|N-N")

# Include aneuploidy values
silver_all_evolution = silver_all_evolution %>% 
  left_join(aneuploidy_pairs_filtered, by = "tumor_pair_barcode")

# Subset to IDHwt, IDHmut noncodel, and IDHmut codel groups for silver set.
silver_all_evolution_IDHwt = silver_all_evolution %>% 
  filter(idh_codel_subtype %in% c("IDHwt_noncodel"))
silver_all_evolution_IDHmut_noncodel = silver_all_evolution %>% 
  filter(idh_codel_subtype %in% c("IDHmut_noncodel"))
silver_all_evolution_IDHmut_codel = silver_all_evolution %>% 
  filter(idh_codel_subtype %in% c("IDHmut_codel"))

# What is the median surgical interval for each subtype x mode of recurrence evolution? Seems to be a difference for IDHwt.
silver_all_evolution %>% 
  group_by(classification_p, idh_codel_subtype) %>% 
  summarise(med_int = median(surgical_interval, na.rm = T),
            sample_size = n())
# While not significant (P > 0.05), here is the breakdown for the aneuploidy differences by evolutionary trajectory.
silver_all_evolution %>% 
  group_by(bay_evo_mode, idh_codel_subtype) %>% 
  summarise(med_aneu = median(aneuploidy_diff, na.rm = T),
            sample_size = n())

# Frequency of selection at each timepoint.
table(silver_all_evolution$classification_p)
table(silver_all_evolution$classification_r)

# Stacked barplot for Bayesian neutral evolution.
stacked_neutrality = silver_all_evolution %>% 
  group_by(idh_codel_subtype, bay_evo_mode) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n))
ggplot(stacked_neutrality, aes(x=idh_codel_subtype, y=freq, fill=bay_evo_mode)) + geom_bar(stat="identity") +
  labs(fill="Mode of evolution") + xlab("") + theme_bw()

# Test by method of **EVOLUTION AT PRIMARY** since this could be explained by other clinical variables.
## IDH MUT CODELs.
fisher.test(table(silver_all_evolution$idh_codel_subtype, silver_all_evolution$bay_evo_mode))

wilcox.test(silver_all_evolution_IDHmut_codel$surgical_interval~silver_all_evolution_IDHmut_codel$classification_p)
fisher.test(table(silver_all_evolution_IDHmut_codel$classification_p, silver_all_evolution_IDHmut_codel$snv_driver_stability))
fisher.test(table(silver_all_evolution_IDHmut_codel$classification_p, silver_all_evolution_IDHmut_codel$cnv_driver_stability))
fisher.test(table(silver_all_evolution_IDHmut_codel$classification_p, silver_all_evolution_IDHmut_codel$recurrence_location))
fisher.test(table(silver_all_evolution_IDHmut_codel$classification_p, silver_all_evolution_IDHmut_codel$received_tmz))
fisher.test(table(silver_all_evolution_IDHmut_codel$classification_p, silver_all_evolution_IDHmut_codel$tmz_std_care))
fisher.test(table(silver_all_evolution_IDHmut_codel$classification_p, silver_all_evolution_IDHmut_codel$received_rt))
fisher.test(table(silver_all_evolution_IDHmut_codel$classification_p, silver_all_evolution_IDHmut_codel$treatment))
fisher.test(table(silver_all_evolution_IDHmut_codel$classification_p, silver_all_evolution_IDHmut_codel$extensive_treatment))
fisher.test(table(silver_all_evolution_IDHmut_codel$classification_p, silver_all_evolution_IDHmut_codel$hypermutator_status))
fisher.test(table(silver_all_evolution_IDHmut_codel$classification_p, silver_all_evolution_IDHmut_codel$grade_change))

# IDH MUT NON-CODELs.
wilcox.test(silver_all_evolution_IDHmut_noncodel$surgical_interval~silver_all_evolution_IDHmut_noncodel$classification_p)
fisher.test(table(silver_all_evolution_IDHmut_noncodel$classification_p, silver_all_evolution_IDHmut_noncodel$snv_driver_stability))
# Interesting, weak association between CNV driver stability and primary evolution.
fisher.test(table(silver_all_evolution_IDHmut_noncodel$classification_p, silver_all_evolution_IDHmut_noncodel$cnv_driver_stability))
fisher.test(table(silver_all_evolution_IDHmut_noncodel$classification_p, silver_all_evolution_IDHmut_noncodel$recurrence_location))
fisher.test(table(silver_all_evolution_IDHmut_noncodel$classification_p, silver_all_evolution_IDHmut_noncodel$received_tmz))
fisher.test(table(silver_all_evolution_IDHmut_noncodel$classification_p, silver_all_evolution_IDHmut_noncodel$tmz_std_care))
fisher.test(table(silver_all_evolution_IDHmut_noncodel$classification_p, silver_all_evolution_IDHmut_noncodel$received_rt))
fisher.test(table(silver_all_evolution_IDHmut_noncodel$classification_p, silver_all_evolution_IDHmut_noncodel$treatment))
fisher.test(table(silver_all_evolution_IDHmut_noncodel$classification_p, silver_all_evolution_IDHmut_noncodel$extensive_treatment))
fisher.test(table(silver_all_evolution_IDHmut_noncodel$classification_p, silver_all_evolution_IDHmut_noncodel$hypermutator_status))
fisher.test(table(silver_all_evolution_IDHmut_noncodel$classification_p, silver_all_evolution_IDHmut_noncodel$grade_change))


# IDH WT - neutrality is association with surgical interval.
wilcox.test(silver_all_evolution_IDHwt$surgical_interval~silver_all_evolution_IDHwt$classification_p)
fisher.test(table(silver_all_evolution_IDHwt$classification_p, silver_all_evolution_IDHwt$snv_driver_stability))
fisher.test(table(silver_all_evolution_IDHwt$classification_p, silver_all_evolution_IDHwt$cnv_driver_stability))
fisher.test(table(silver_all_evolution_IDHwt$classification_p, silver_all_evolution_IDHwt$recurrence_location))
fisher.test(table(silver_all_evolution_IDHwt$classification_p, silver_all_evolution_IDHwt$received_tmz))
fisher.test(table(silver_all_evolution_IDHwt$classification_p, silver_all_evolution_IDHwt$tmz_std_care))
fisher.test(table(silver_all_evolution_IDHwt$classification_p, silver_all_evolution_IDHwt$received_rt))
fisher.test(table(silver_all_evolution_IDHwt$classification_p, silver_all_evolution_IDHwt$treatment))
fisher.test(table(silver_all_evolution_IDHwt$classification_p, silver_all_evolution_IDHwt$hypermutator_status))
fisher.test(table(silver_all_evolution_IDHwt$classification_p, silver_all_evolution_IDHwt$grade_change))
wilcox.test(silver_all_evolution_IDHwt$snv_driver_count~silver_all_evolution_IDHwt$classification_p)
wilcox.test(silver_all_evolution_IDHwt$snv_driver_count_private_a~silver_all_evolution_IDHwt$classification_p)


# Test by method of **EVOLUTION AT RECURRENCE** since this could be explained by other clinical variables.
## IDH MUT CODELs.
wilcox.test(silver_all_evolution_IDHmut_codel$surgical_interval~silver_all_evolution_IDHmut_codel$classification_r)
wilcox.test(silver_all_evolution_IDHmut_codel$aneuploidy_diff~silver_all_evolution_IDHmut_codel$classification_r)
fisher.test(table(silver_all_evolution_IDHmut_codel$classification_r, silver_all_evolution_IDHmut_codel$snv_driver_stability))
fisher.test(table(silver_all_evolution_IDHmut_codel$classification_r, silver_all_evolution_IDHmut_codel$cnv_driver_stability))
fisher.test(table(silver_all_evolution_IDHmut_codel$classification_r, silver_all_evolution_IDHmut_codel$recurrence_location))
fisher.test(table(silver_all_evolution_IDHmut_codel$classification_r, silver_all_evolution_IDHmut_codel$received_tmz))
fisher.test(table(silver_all_evolution_IDHmut_codel$classification_r, silver_all_evolution_IDHmut_codel$tmz_std_care))
fisher.test(table(silver_all_evolution_IDHmut_codel$classification_r, silver_all_evolution_IDHmut_codel$received_rt))
fisher.test(table(silver_all_evolution_IDHmut_codel$classification_r, silver_all_evolution_IDHmut_codel$treatment))
fisher.test(table(silver_all_evolution_IDHmut_codel$classification_r, silver_all_evolution_IDHmut_codel$extensive_treatment))
fisher.test(table(silver_all_evolution_IDHmut_codel$classification_r, silver_all_evolution_IDHmut_codel$hypermutator_status))
fisher.test(table(silver_all_evolution_IDHmut_codel$classification_r, silver_all_evolution_IDHmut_codel$grade_change))

# IDH MUT NON-CODELs.
wilcox.test(silver_all_evolution_IDHmut_noncodel$surgical_interval~silver_all_evolution_IDHmut_noncodel$classification_r)
wilcox.test(silver_all_evolution_IDHmut_noncodel$aneuploidy_diff~silver_all_evolution_IDHmut_noncodel$classification_r)
fisher.test(table(silver_all_evolution_IDHmut_noncodel$classification_r, silver_all_evolution_IDHmut_noncodel$snv_driver_stability))
fisher.test(table(silver_all_evolution_IDHmut_noncodel$classification_r, silver_all_evolution_IDHmut_noncodel$cnv_driver_stability))
fisher.test(table(silver_all_evolution_IDHmut_noncodel$classification_r, silver_all_evolution_IDHmut_noncodel$recurrence_location))
fisher.test(table(silver_all_evolution_IDHmut_noncodel$classification_r, silver_all_evolution_IDHmut_noncodel$received_tmz))
fisher.test(table(silver_all_evolution_IDHmut_noncodel$classification_r, silver_all_evolution_IDHmut_noncodel$tmz_std_care))
fisher.test(table(silver_all_evolution_IDHmut_noncodel$classification_r, silver_all_evolution_IDHmut_noncodel$received_rt))
fisher.test(table(silver_all_evolution_IDHmut_noncodel$classification_r, silver_all_evolution_IDHmut_noncodel$treatment))
fisher.test(table(silver_all_evolution_IDHmut_noncodel$classification_r, silver_all_evolution_IDHmut_noncodel$extensive_treatment))
fisher.test(table(silver_all_evolution_IDHmut_noncodel$classification_r, silver_all_evolution_IDHmut_noncodel$hypermutator_status))
fisher.test(table(silver_all_evolution_IDHmut_noncodel$classification_r, silver_all_evolution_IDHmut_noncodel$grade_change))

# IDH WT - neutrality is association with surgical interval.
wilcox.test(silver_all_evolution_IDHwt$surgical_interval~silver_all_evolution_IDHwt$classification_r)
wilcox.test(silver_all_evolution_IDHwt$aneuploidy_diff~silver_all_evolution_IDHwt$classification_r)
fisher.test(table(silver_all_evolution_IDHwt$classification_r, silver_all_evolution_IDHwt$snv_driver_stability))
fisher.test(table(silver_all_evolution_IDHwt$classification_r, silver_all_evolution_IDHwt$cnv_driver_stability))
fisher.test(table(silver_all_evolution_IDHwt$classification_r, silver_all_evolution_IDHwt$recurrence_location))
fisher.test(table(silver_all_evolution_IDHwt$classification_r, silver_all_evolution_IDHwt$received_tmz))
fisher.test(table(silver_all_evolution_IDHwt$classification_r, silver_all_evolution_IDHwt$tmz_std_care))
fisher.test(table(silver_all_evolution_IDHwt$classification_r, silver_all_evolution_IDHwt$received_rt))
fisher.test(table(silver_all_evolution_IDHwt$classification_r, silver_all_evolution_IDHwt$treatment))
fisher.test(table(silver_all_evolution_IDHwt$classification_r, silver_all_evolution_IDHwt$extensive_treatment))
fisher.test(table(silver_all_evolution_IDHwt$classification_r, silver_all_evolution_IDHwt$hypermutator_status))
fisher.test(table(silver_all_evolution_IDHwt$evolution_area_pval_p, silver_all_evolution_IDHwt$hypermutator_status))
fisher.test(table(silver_all_evolution_IDHwt$evolution_area_pval_r, silver_all_evolution_IDHwt$hypermutator_status))
fisher.test(table(silver_all_evolution_IDHwt$freq_evo_mode, silver_all_evolution_IDHwt$hypermutator_status))

fisher.test(table(silver_all_evolution_IDHwt$classification_r, silver_all_evolution_IDHwt$grade_change))
wilcox.test(silver_all_evolution_IDHwt$snv_driver_count~silver_all_evolution_IDHwt$classification_r)
wilcox.test(silver_all_evolution_IDHwt$snv_driver_count_private_b~silver_all_evolution_IDHwt$classification_r)

# Test by method of **EVOLUTION MODE** since this could be explained by other clinical variables.
## IDH MUT CODELs.
wilcox.test(silver_all_evolution_IDHmut_codel$surgical_interval~silver_all_evolution_IDHmut_codel$bay_evo_mode)
kruskal.test(silver_all_evolution_IDHmut_codel$aneuploidy_diff~factor(silver_all_evolution_IDHmut_codel$bay_evo_mode))
fisher.test(table(silver_all_evolution_IDHmut_codel$bay_evo_mode, silver_all_evolution_IDHmut_codel$snv_driver_stability))
fisher.test(table(silver_all_evolution_IDHmut_codel$bay_evo_mode, silver_all_evolution_IDHmut_codel$cnv_driver_stability))
fisher.test(table(silver_all_evolution_IDHmut_codel$bay_evo_mode, silver_all_evolution_IDHmut_codel$recurrence_location))
fisher.test(table(silver_all_evolution_IDHmut_codel$bay_evo_mode, silver_all_evolution_IDHmut_codel$received_tmz))
fisher.test(table(silver_all_evolution_IDHmut_codel$bay_evo_mode, silver_all_evolution_IDHmut_codel$tmz_std_care))
fisher.test(table(silver_all_evolution_IDHmut_codel$bay_evo_mode, silver_all_evolution_IDHmut_codel$received_rt))
fisher.test(table(silver_all_evolution_IDHmut_codel$bay_evo_mode, silver_all_evolution_IDHmut_codel$treatment))
fisher.test(table(silver_all_evolution_IDHmut_codel$bay_evo_mode, silver_all_evolution_IDHmut_codel$hypermutator_status))
fisher.test(table(silver_all_evolution_IDHmut_codel$bay_evo_mode, silver_all_evolution_IDHmut_codel$grade_change))

# IDH MUT NON-CODELs.
wilcox.test(silver_all_evolution_IDHmut_noncodel$surgical_interval~silver_all_evolution_IDHmut_noncodel$bay_evo_mode)
kruskal.test(silver_all_evolution_IDHmut_noncodel$aneuploidy_diff~factor(silver_all_evolution_IDHmut_noncodel$bay_evo_mode))
fisher.test(table(silver_all_evolution_IDHmut_noncodel$bay_evo_mode, silver_all_evolution_IDHmut_noncodel$cnv_driver_stability))
fisher.test(table(silver_all_evolution_IDHmut_noncodel$bay_evo_mode, silver_all_evolution_IDHmut_noncodel$recurrence_location))
fisher.test(table(silver_all_evolution_IDHmut_noncodel$bay_evo_mode, silver_all_evolution_IDHmut_noncodel$received_tmz))
fisher.test(table(silver_all_evolution_IDHmut_noncodel$bay_evo_mode, silver_all_evolution_IDHmut_noncodel$tmz_std_care))
fisher.test(table(silver_all_evolution_IDHmut_noncodel$bay_evo_mode, silver_all_evolution_IDHmut_noncodel$received_rt))
fisher.test(table(silver_all_evolution_IDHmut_noncodel$bay_evo_mode, silver_all_evolution_IDHmut_noncodel$treatment))
fisher.test(table(silver_all_evolution_IDHmut_noncodel$bay_evo_mode, silver_all_evolution_IDHmut_noncodel$hypermutator_status))
fisher.test(table(silver_all_evolution_IDHmut_noncodel$bay_evo_mode, silver_all_evolution_IDHmut_noncodel$grade_change))


# IDH WT - neutrality is association with surgical interval.
wilcox.test(silver_all_evolution_IDHwt$surgical_interval~silver_all_evolution_IDHwt$bay_evo_mode)
kruskal.test(silver_all_evolution_IDHwt$aneuploidy_diff~factor(silver_all_evolution_IDHwt$bay_evo_mode))
fisher.test(table(silver_all_evolution_IDHwt$bay_evo_mode, silver_all_evolution_IDHwt$snv_driver_stability))
fisher.test(table(silver_all_evolution_IDHwt$bay_evo_mode, silver_all_evolution_IDHwt$cnv_driver_stability))
fisher.test(table(silver_all_evolution_IDHwt$bay_evo_mode, silver_all_evolution_IDHwt$recurrence_location))
fisher.test(table(silver_all_evolution_IDHwt$bay_evo_mode, silver_all_evolution_IDHwt$received_tmz))
fisher.test(table(silver_all_evolution_IDHwt$bay_evo_mode, silver_all_evolution_IDHwt$tmz_std_care))
fisher.test(table(silver_all_evolution_IDHwt$bay_evo_mode, silver_all_evolution_IDHwt$received_rt))
fisher.test(table(silver_all_evolution_IDHwt$bay_evo_mode, silver_all_evolution_IDHwt$treatment))
fisher.test(table(silver_all_evolution_IDHwt$bay_evo_mode, silver_all_evolution_IDHwt$hypermutator_status))
fisher.test(table(silver_all_evolution_IDHwt$bay_evo_mode, silver_all_evolution_IDHwt$hypermutator_status))
fisher.test(table(silver_all_evolution_IDHwt$bay_evo_mode, silver_all_evolution_IDHwt$grade_change))
wilcox.test(silver_all_evolution_IDHwt$snv_driver_count~silver_all_evolution_IDHwt$bay_evo_mode)
wilcox.test(silver_all_evolution_IDHwt$snv_driver_count_private_b~silver_all_evolution_IDHwt$bay_evo_mode)

#### IDH WT.
# *** Neutrality at primary is associated with OS survival.
fit_wt_primary <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ classification_p,
                          data = silver_all_evolution_IDHwt)
ggsurvplot(fit_wt_primary, data = silver_all_evolution_IDHwt, risk.table = FALSE, pval= TRUE, pval.method = TRUE, surv.median.line = "v", palette = c("royalblue4", "tomato3"))
fit_wt_recur <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ classification_r,
                        data = silver_all_evolution_IDHwt)
ggsurvplot(fit_wt_recur, data = silver_all_evolution_IDHwt, risk.table = TRUE, pval= TRUE, surv.median.line = "v")
fit_wt_post_recur <- survfit(Surv(patient_post_recur_surv, patient_vital) ~ classification_r,
                             data = silver_all_evolution_IDHwt)
ggsurvplot(fit_wt_post_recur, data = silver_all_evolution_IDHwt, risk.table = TRUE, pval= TRUE, surv.median.line = "v")

## Create a publication quality plot for the neutral-evolution analyses.
fit_wt_mode <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ bay_evo_mode,
                       data = silver_all_evolution_IDHwt)
pdf(file = "/Users/johnsk/Documents/survival-analysis-neutral.pdf", height = 6, width = 8, bg = "transparent", useDingbats = FALSE)
ggsurvplot(fit_wt_mode, data = silver_all_evolution_IDHwt, risk.table = FALSE, pval= TRUE, pval.method = FALSE, pval.coord = c(100, 0.75),
           surv.median.line = "v", palette = c("#67A3BD", "#8167BD", "#A3BD67", "#BD8167"), legend = "none")
dev.off()
pdf(file = "/Users/johnsk/Documents/neutral-legend.pdf", height = 6, width = 8, bg = "transparent", useDingbats = FALSE)
ggsurvplot(fit_wt_mode, data = silver_all_evolution_IDHwt, risk.table = FALSE, pval= TRUE, pval.method = FALSE,
           surv.median.line = "v", palette = c("#67A3BD", "#8167BD", "#A3BD67", "#BD8167"), legend = "right", legend.labs =c("N-N", "N-S", "S-N", "S-S"))
dev.off()


fit_wt_mode <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ binary_mode,
                        data = silver_all_evolution_IDHwt)
ggsurvplot(fit_wt_mode, data = silver_all_evolution_IDHwt, risk.table = FALSE, pval= TRUE, pval.method = TRUE, surv.median.line = "v")


#### IDHmut-noncodel
fit_noncodel_primary <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ classification_p,
                                data = silver_all_evolution_IDHmut_noncodel)
ggsurvplot(fit_noncodel_primary, data = silver_all_evolution_IDHmut_noncodel, risk.table = TRUE, pval= TRUE, surv.median.line = "v")
fit_noncodel_recur <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ classification_r,
                              data = silver_all_evolution_IDHmut_noncodel)
ggsurvplot(fit_noncodel_recur, data = silver_all_evolution_IDHmut_noncodel, risk.table = TRUE, pval= TRUE, surv.median.line = "v")
fit_noncodel_recur <- survfit(Surv(patient_post_recur_surv, patient_vital) ~ classification_r,
                          data = silver_all_evolution_IDHmut_noncodel)
ggsurvplot(fit_noncodel_recur, data = silver_all_evolution_IDHmut_noncodel, risk.table = TRUE, pval= TRUE, surv.median.line = "v")
fit_noncodel_mode <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ binary_mode,
                                data = silver_all_evolution_IDHmut_noncodel)
ggsurvplot(fit_noncodel_mode, data = silver_all_evolution_IDHmut_noncodel, risk.table = FALSE, pval= TRUE, pval.method = TRUE, surv.median.line = "v")

##### IDHwt
####### Cox proportional hazards model #####
res_cox_idh_wt <- coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + classification_p, data = silver_all_evolution_IDHwt)
summary(res_cox_idh_wt)
res_cox_idh_wt <- coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + classification_r, data = silver_all_evolution_IDHwt)
summary(res_cox_idh_wt)
res_cox_idh_wt <- coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + bay_evo_mode, data = silver_all_evolution_IDHwt)
summary(res_cox_idh_wt)
res_cox_idh_wt <- coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + binary_mode, data = silver_all_evolution_IDHwt)
summary(res_cox_idh_wt)


# PRIMARY evolution is not significant when included primary evolution as a predictor in the cox model.
res_cox_drivers <- coxph(Surv(case_overall_survival_mo, patient_vital)~ case_age_diagnosis_years + idh_codel_subtype + binary_mode + hypermutator_status, data = silver_all_evolution)
summary(res_cox_drivers)


##################################################
# For those samples where selection was indicated
# import the other selection information.
subclonal_model1 = read_tsv("/Users/johnsk/Documents/Life-History/glass-analyses/data/SubClonalSelection/2019.01.29-results/formatted_metadata3_parameters_model1.txt")
# Filter again to only use "all" samples.
subclonal_model1_results = subclonal_model1 %>% 
separate(sample, c("aliquot_barcode", "mutation_type"), sep ="_") %>% 
  filter(mutation_type == "all") 

# Restrict to any selected samples.
any_selected_tumors = silver_all_evolution %>% 
  left_join(subclonal_model1_results, by=c("tumor_barcode_a_p"="aliquot_barcode")) %>% 
  left_join(subclonal_model1_results, by=c("tumor_barcode_b_r"="aliquot_barcode")) %>% 
  filter(classification_p == "S" | classification_r == "S")

# Separate out primary from recurrence.
selected_primary_tumor = any_selected_tumors %>% 
  filter(classification_p == "S")
selected_primary_tumor %>% 
  group_by(idh_codel_subtype) %>% 
  summarise(avg_s = median(s.x),
            avg_t = median(t.x), 
            sample_size = n())
ggplot(selected_primary_tumor, aes(x= idh_codel_subtype, y = s.x)) + geom_boxplot() + ylim(0.5, 2)
ggplot(selected_primary_tumor, aes(x= idh_codel_subtype, y = t.x)) + geom_boxplot() + ylim(7, 14)
# Again, examine the average number of tumor doublings and selection.
selected_recurrence_tumor = any_selected_tumors %>% 
  filter(classification_r == "S")
selected_recurrence_tumor %>% 
  group_by(idh_codel_subtype) %>% 
  summarise(avg_s = median(s.y),
            avg_t = median(t.y), 
            sample_size = n())
ggplot(selected_recurrence_tumor, aes(x= idh_codel_subtype, y = s.y)) + geom_boxplot() + ylim(0.5, 2)
ggplot(selected_recurrence_tumor, aes(x= idh_codel_subtype, y = t.y)) + geom_boxplot() + ylim(7, 14)
 
wilcox.test(selected_recurrence_tumor$t.y, selected_primary_tumor$t.x)
wilcox.test(selected_recurrence_tumor$s.y, selected_primary_tumor$s.x)
#####################################################
### Extra
#####################################################
# Write out select columns fro Fred's analysis:
colnames(silver_all_evolution)
evolution_out = silver_all_evolution %>% 
  select(tumor_pair_barcode, frequentist_evo_p = evolution_area_pval_p, frequentist_evo_r = evolution_area_pval_r,
         bayesian_evo_p = classification_p, bayesian_evo_r = classification_r, idh_codel_subtype)
write.table(evolution_out, file = "/Users/johnsk/Documents/glass-evolution-neutrality-20190129.txt", sep="\t", row.names = F, col.names = T, quote = F)

#####################################################
#### Sample breakdown
#####################################################
# Median overall survival and group size:
silver_all_evolution %>% 
  group_by(bay_evo_mode) %>% 
  summarise(median_surv = median(case_overall_survival_mo, na.rm = T),
            sample_size = n())

silver_all_evolution %>% 
  group_by(bay_evo_mode, idh_codel_subtype) %>% 
  summarise(sample_size = n())
# Test whether there is a difference in evo mode by subtype.
fisher.test(table(silver_all_evolution$bay_evo_mode, silver_all_evolution$idh_codel_subtype))
# Collapse information so that it's IDHwt vs. IDHmut.
IDHwt_IDHmut_evo = matrix(c(36, 12, 4, 7, 40, 7, 9, 12), nrow = 4)
fisher.test(IDHwt_IDHmut_evo) # No significant differences.

silver_all_evolution %>% 
  group_by(hypermutator_status, idh_codel_subtype, bay_evo_mode) %>% 
  summarise(sample_size = n())

silver_all_evolution %>% 
  group_by(bay_evo_mode, idh_codel_subtype, hypermutator_status) %>% 
  summarise(sample_size = n())

