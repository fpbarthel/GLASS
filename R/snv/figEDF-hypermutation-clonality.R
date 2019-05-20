##################################################
# Investigate the clonality of hypermutant recurrent specific mutations
# Updated: 2019.05.20
# Author: Kevin J.
##################################################

# Working directory for this analysis in the GLASS-analysis project. 
mybasedir = "/Volumes/verhaak-lab/GLASS-analysis/"
setwd(mybasedir)

##################################################
# Necessary packages:
library(tidyverse)
library(DBI)
library(colortools)
library(survminer)
library(survival)
library(RColorBrewer)

##################################################
# Connect to version 2 of the database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB2")

# Load necessary tables.
subtypes = dbReadTable(con,  Id(schema="clinical", table="subtypes"))
tumor_pairs = dbReadTable(con,  Id(schema="analysis", table="tumor_pairs"))
mutation_freq = dbGetQuery(con, "SELECT * FROM analysis.mut_freq")
passanno = dbReadTable(con,  Id(schema="variants", table="passanno"))
mut_freq = dbReadTable(con,  Id(schema="analysis", table="mut_freq"))
variant_classifications = dbReadTable(con,  Id(schema="variants", table="variant_classifications"))
gold_set = dbReadTable(con,  Id(schema="analysis", table="gold_set"))
tumor_pairs = dbReadTable(con,  Id(schema="analysis", table="tumor_pairs"))
cases = dbReadTable(con,  Id(schema="clinical", table="cases"))

# Get the new clinical_tumor_pairs table.
clinical_tumor_pairs_query = read_file("sql/clinical/clinical-tumor-pairs-db2.sql")
clinical_tumor_pairs <- dbGetQuery(con, clinical_tumor_pairs_query)


# Floris ran PyClone on all 30X mutations from the GLASS samples. Data not available for all samples as some failed.
# Classify variants as "clonal" if CP > 0.5, variants <= 0.5 and > 0.1 are subclonal, and non-clonal (no classification) below 0.1.
pyclone_variants = dbGetQuery(con, "WITH t1 AS (SELECT 
	pg.tumor_pair_barcode,
                              pg.case_barcode,
                              pg.tumor_barcode_a,
                              pg.tumor_barcode_b,
                              pg.chrom, 
                              pg.pos,
                              pg.variant_id,
                              pg.mutect2_call_a,
                              pg.mutect2_call_b, 
                              pl1.cellular_prevalence AS cellular_prevalence_a, 
                              pl1.variant_allele_frequency AS variant_allele_frequency_a, 
                              (CASE WHEN pl1.cellular_prevalence >= 0.5 THEN 'C' WHEN pl1.cellular_prevalence >= 0.1 AND pl1.cellular_prevalence < 0.5 AND pl1.cellular_prevalence < 0.5 THEN 'S' END) AS clonality_a,
                              pl2.cellular_prevalence AS cellular_prevalence_b, 
                              pl2.variant_allele_frequency AS variant_allele_frequency_b, 
                              (CASE WHEN pl2.cellular_prevalence >= 0.5 THEN 'C' WHEN pl2.cellular_prevalence >= 0.1 AND pl2.cellular_prevalence < 0.5 AND pl2.cellular_prevalence < 0.5 THEN 'S' END) AS clonality_b,
                              (CASE WHEN mutect2_call_a AND mutect2_call_b THEN 'S'
                              WHEN mutect2_call_a AND NOT mutect2_call_b THEN 'P'
                              WHEN mutect2_call_b AND NOT mutect2_call_a THEN 'R' END) AS fraction
                              FROM variants.pgeno pg
                              LEFT JOIN variants.pyclone_loci pl1 ON pl1.variant_id = pg.variant_id AND pl1.aliquot_barcode = pg.tumor_barcode_a 
                              LEFT JOIN variants.pyclone_loci pl2 ON pl2.variant_id = pg.variant_id AND pl2.aliquot_barcode= pg.tumor_barcode_b
                              INNER JOIN analysis.gold_set gs ON pg.tumor_pair_barcode = gs.tumor_pair_barcode
                              WHERE pl1.cellular_prevalence IS NOT NULL) 
                              
                              SELECT * FROM t1 WHERE fraction = 'R'")

# How many of the silver set pairs are included in this analysis?
n_distinct(pyclone_variants$tumor_pair_barcode) # 212

# Merge the passanno with the pyclone_variants.
pyclone_annot = pyclone_variants %>% 
  left_join(passanno, by="variant_id")

# Merge with the other annotation information.
recurrent_pyclone_variants = pyclone_annot %>%
  left_join(clinical_tumor_pairs, by=c("tumor_pair_barcode", "case_barcode", "tumor_barcode_a", "tumor_barcode_b")) %>% 
  left_join(subtypes, by="case_barcode")

# Subset to IDHmut-noncodel, hypermutators.
idh_mut_hyper_clonality = recurrent_pyclone_variants %>% 
  filter(idh_codel_subtype == "IDHmut-noncodel", alk_assoc_hypermutator_status == 1)

# Stacked clonality:
stacked_pyclone_idh_hyper  = idh_mut_hyper_clonality %>% 
  group_by(clonality_b, case_barcode) %>% 
  summarise (n = n()) %>%
  mutate(freq = n / sum(n))

# Visually inspect the distribution of clonality for the IDHmut-noncodels.
ggplot(stacked_pyclone_idh_hyper, aes(x=case_barcode, y=n, fill=clonality_b)) + geom_bar(stat="identity") +
  labs(fill="Mutation Status") + xlab("") + theme_bw() + ylab("Recurrent s of Mutations") + ggtitle("IDHmut-noncodel hypermutators") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Determine for each sample which of the recurrent-only mutation clonality classes predominates.
hyper_classes_noncodels = stacked_pyclone_idh_hyper %>% 
  group_by(case_barcode) %>% 
  slice(which.max(n)) %>% 
  mutate(clonality_b = ifelse(is.na(clonality_b),"S", clonality_b)) 

# Extract the relevant clinical information and combine with the mut_freq table.
clinical_tumor_gold = clinical_tumor_pairs %>%  
  filter(tumor_pair_barcode %in% gold_set$tumor_pair_barcode) %>% 
  left_join(subtypes, by="case_barcode") %>% 
  left_join(cases, by="case_barcode") %>% 
  left_join(mut_freq, by=c("tumor_barcode_b" = "aliquot_barcode")) %>% 
  mutate(patient_vital =  ifelse(case_vital_status=="alive", 0, 1))

# Now, restrict to the IDHmut-noncodels that have received any form of alyklating agent therapy.
clinical_tumor_gold_IDHmut_noncodel = clinical_tumor_gold %>% 
  filter(idh_codel_subtype == "IDHmut-noncodel") %>% 
  filter(received_alk == 1)

# Create a df that classifies hypermutants based on whether they are clonal or subclonal hypermutation events.
idh_mut_hyper_clone = clinical_tumor_gold_IDHmut_noncodel %>% 
  left_join(hyper_classes_noncodels, by="case_barcode") %>% 
  mutate(hyper_clone = ifelse(is.na(clonality_b) & alk_assoc_hypermutator_status==0, "non-hyper", clonality_b)) 
idh_mut_hyper_clone$hyper_clone = factor(idh_mut_hyper_clone$hyper_clone, levels=c("C", "S", "non-hyper"))

# Assess whether there is any difference in survival between these three groups.
fit_IDHmut_surv <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ hyper_clone,
                          data = idh_mut_hyper_clone)
pdf(file = "/Users/johnsk/Documents/idh-noncodel-hypermutator-clonality-surv.pdf", height = 5, width = 7, bg = "transparent", useDingbats = FALSE)
ggsurvplot(fit_IDHmut_surv, data = tmp, risk.table = FALSE, pval= TRUE, pval.coord = c(100, 0.50),
           surv.median.line = "v", ylab = "Overall survival \n probability", xlab = "Time (months)", palette = c("#e31a1c", "#33a02c", "#1f78b4")) 
dev.off()

# Does this change when age is included as a variable? Answer: No.
hyper_cox_model = coxph(Surv(case_overall_survival_mo, patient_vital) ~ case_age_diagnosis_years + hyper_clone, data = tmp)
summary(hyper_cox_model) # P = 0.93

# Subset to IDHmut-codel, hypermutators.
idh_mut_codel_hyper_clonality = recurrent_pyclone_variants %>% 
  filter(idh_codel_subtype == "IDHmut-codel", alk_assoc_hypermutator_status == 1)
# Create stack.
stacked_pyclone_idh_codel_hyper  = idh_mut_codel_hyper_clonality %>% 
  group_by(clonality_b, case_barcode) %>% 
  summarise (n = n()) %>%
  mutate(freq = n / sum(n))
# Plot restricted to IDHmut-codels.
ggplot(stacked_pyclone_idh_codel_hyper, aes(x=case_barcode, y=n, fill=clonality_b)) + geom_bar(stat="identity") +
  labs(fill="Mutation Status") + xlab("") + theme_bw() + ylab("Recurrent # of Mutations") + ggtitle("IDHmut-codel hypermutators") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# Subset to IDHwt.
idh_wt_hyper_clonality = recurrent_pyclone_variants %>% 
  filter(idh_codel_subtype == "IDHwt", alk_assoc_hypermutator_status == 1)
# Create IDHwt stack.
stacked_pyclone_idh_wt_hyper  = idh_wt_hyper_clonality %>% 
  group_by(clonality_b, case_barcode) %>% 
  summarise (n = n()) %>%
  mutate(freq = n / sum(n))
# Plot the samples with PyClone and Hypermutation data.
ggplot(stacked_pyclone_idh_wt_hyper, aes(x=case_barcode, y=n, fill=clonality_b)) + geom_bar(stat="identity") +
  labs(fill="Mutation Status") + xlab("") + theme_bw() + ylab("Recurrent # of Mutations") + ggtitle("IDHwt hypermutators") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


# Combine all hypermutators in the same analysis.
all_hypers = recurrent_pyclone_variants %>% 
  filter(alk_assoc_hypermutator_status == 1) %>% 
  mutate(tmp1 = ifelse(is.na(clonality_b),"Non-clonal", clonality_b)) %>% 
  mutate(clonality_recurrence = recode(tmp1, "C" = "Clonal", "S" = "Subclonal"))
all_hypers$clonality_recurrence = factor(all_hypers$clonality_recurrence, levels = c("Clonal", "Subclonal", "Non-clonal"))
# Prepare the data for the same plot faceted across the subtypes.
stacked_all_hypers  = all_hypers %>% 
  group_by(clonality_recurrence, case_barcode) %>% 
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  left_join(subtypes, by="case_barcode")

pdf(file = "/Users/johnsk/Documents/hypermutants-recurrence-clonality.pdf", height = 5, width = 15, bg = "transparent", useDingbats = FALSE)
p = ggplot(stacked_all_hypers, aes(x=case_barcode, y=n, fill=clonality_recurrence)) + geom_bar(stat="identity") +
  labs(fill="Mutation Status") + xlab("") + theme_bw() + ylab("Privaate to R Mutations \n Hypermutators") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(~idh_codel_subtype, scales = "free", space="free_x") + scale_fill_manual(values = c("Clonal" = "#e31a1c", "Subclonal" = "#33a02c", "Non-clonal" = "#1f78b4"))
dev.off()


# Tally the total number of mutations. Some hypermutators may not share the exact same phenotype (GLSS-CU-R014 vs. GLSS-DF-3081).
all_hypers %>% 
  group_by(case_barcode) %>% 
  summarise (n = n())

# Insert whether the clonal mutations are enriched for particular variant effects.
stacked_all_hypers_class  = all_hypers %>% 
  left_join(variant_classifications, by=c("variant_classification"="variant_classification")) %>% 
  mutate(mut_class = ifelse(variant_classification_impact%in%c("HIGH", "MODERATE"), "nonsynonymous", ifelse(variant_classification_impact%in%c("LOW"), "synonymous", NA))) %>% 
  group_by(variant_effect, case_barcode) %>% 
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  left_join(subtypes, by="case_barcode")

# A coding vs. non-coding plot for the hypermutators.
pdf(file = "/Users/johnsk/Documents/hypermutants-variant-classification.pdf", height = 5, width = 15, bg = "transparent", useDingbats = FALSE)
ggplot(stacked_all_hypers_class, aes(x=case_barcode, y=n, fill=variant_effect)) + geom_bar(stat="identity") +
  labs(fill="Mutation Status") + xlab("") + theme_bw() + ylab("Private to R Mutations \n Hypermutators") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(~idh_codel_subtype, scales = "free", space="free_x") + scale_fill_brewer(palette="Dark2")
dev.off()

# Instead, offer a visualization about the breakdown clonally and functionally for these hypermutation events.
stacked_hypers_clonality_class  = all_hypers %>% 
  left_join(variant_classifications, by=c("variant_classification"="variant_classification")) %>% 
  mutate(clonal_class = paste(clonality_recurrence, variant_effect, sep=":")) %>% 
  group_by(clonal_class, case_barcode) %>% 
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  left_join(subtypes, by="case_barcode")
stacked_hypers_clonality_class$clonal_class = factor(stacked_hypers_clonality_class$clonal_class, levels = c("Clonal:Coding", "Clonal:Non-coding", "Subclonal:Coding", "Subclonal:Non-coding", "Non-clonal:Coding", "Non-clonal:Non-coding"))


pdf(file = "/Users/johnsk/Documents/hypermutants-clonal-classification.pdf", height = 7, width = 15, bg = "transparent", useDingbats = FALSE)
ggplot(stacked_hypers_clonality_class, aes(x=case_barcode, y=n, fill=clonal_class)) + geom_bar(stat="identity") +
  labs(fill="Mutation Status") + xlab("") + theme_bw() + ylab("Privaate to R Mutations \n Hypermutators") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(~idh_codel_subtype, scales = "free", space="free_x") + scale_fill_manual(values = c("Clonal:Coding" = "#e31a1c", "Clonal:Non-coding" = "#fb9a99",
                                                                                                 "Subclonal:Coding" = "#33a02c", "Subclonal:Non-coding" = "#b2df8a",
                                                                                                 "Non-clonal:Coding" = "#1f78b4", "Non-clonal:Non-coding" = "#a6cee3"), name = "Mutation Class")
dev.off()


