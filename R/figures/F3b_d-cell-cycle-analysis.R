#######################
# Cell cycle analysis #
#######################

########################
# Do necessary settings
########################
setwd("/Users/c-kocake/Box Sync/Scripts/repos/GLASS/")
# Load necessary packages
library(tidyverse)
library(DBI)
library(ggpubr)
library(survival)
library(survminer)

# Create connection to database and load necessary tables
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

# Load in variables from the database.
subtypes = dbReadTable(con,  Id(schema="clinical", table="subtypes"))
tumor_blocklist = dbReadTable(con,  Id(schema="analysis", table="blocklist"))
tumor_pairs = dbReadTable(con,  Id(schema="analysis", table="tumor_pairs"))
gold_set = dbReadTable(con,  Id(schema="analysis", table="gold_set"))
os = dbReadTable(con,  Id(schema="clinical", table="cases"))
cnv_by_gene = dbReadTable(con,  Id(schema="analysis", table="gatk_cnv_by_gene"))
mutation_freq = dbGetQuery(con, "SELECT * FROM analysis.mut_freq")

# load in mutations per aliquot from GLASS
q <- "SELECT aliquot_barcode, case_barcode, pa.chrom, lower(pa.pos) as start, upper(pa.pos) - 1 as end, ad_ref, ad_alt, af, ssm2_pass_call, gene_symbol, pa.variant_classification, variant_type FROM variants.passgeno pg
INNER JOIN variants.passanno pa ON pa.variant_id = pg.variant_id
INNER JOIN variants.variant_classifications vs ON vs.variant_classification = pa.variant_classification
WHERE vs.variant_classification_priority IS NOT NULL AND ssm2_pass_call IS TRUE"
mut_df <- DBI::dbGetQuery(con, q) 

# load in RB1 homozygous deletions from GLASS
r <- "WITH selected_genes AS (
  SELECT ge.gene_symbol,
ge.chrom,
ge.pos
FROM ref.genes ge
WHERE gene_symbol = 'RB1'
), gene_seg_intersect AS (
SELECT gs.aliquot_barcode,
t0.gene_symbol,
gs.chrom,
upper(t0.pos * gs.pos) - lower(t0.pos * gs.pos) - 1 AS w,
2::numeric ^ gs.log2_copy_ratio AS cr
FROM variants.gatk_seg gs
JOIN selected_genes t0 ON t0.chrom = gs.chrom AND t0.pos && gs.pos
), gene_sample_call AS (
SELECT gene_seg_intersect.aliquot_barcode,
gene_seg_intersect.gene_symbol,
sum(gene_seg_intersect.w::numeric * gene_seg_intersect.cr) / sum(gene_seg_intersect.w)::numeric AS wcr
FROM gene_seg_intersect
GROUP BY gene_seg_intersect.aliquot_barcode, gene_seg_intersect.gene_symbol
), seg_stats_optimized AS (
SELECT gs.aliquot_barcode,
LEAST(0.9, gs.neu_fwmean - 2::numeric * gs.neu_fwsd) AS del_thres,
GREATEST(1.1, gs.neu_fwmean + 2::numeric * gs.neu_fwsd) AS amp_thres,
CASE
WHEN gsa.max_loss_arm_wmean < 0.9 AND gsa.max_loss_arm_n >= 3 THEN GREATEST(0::numeric, gsa.max_loss_arm_wmean - 2::numeric * gsa.max_loss_arm_wsd)
WHEN gs.del_fwmean < 0.9 AND gs.del_n >= 3 THEN GREATEST(0::numeric, gs.del_fwmean - 2::numeric * gs.del_fwsd)
ELSE NULL::numeric
END AS hldel_thres,
CASE
WHEN gsa.max_gain_arm_wmean > 1.1 AND gsa.max_gain_arm_n >= 3 THEN gsa.max_gain_arm_wmean + 2::numeric * gsa.max_gain_arm_wsd
WHEN gs.amp_fwmean > 1.1 AND gs.amp_n >= 3 THEN gs.amp_fwmean + 2::numeric * gs.amp_fwsd
ELSE NULL::numeric
END AS hlamp_thres
FROM analysis.gatk_seg_stats gs
LEFT JOIN analysis.gatk_aneuploidy gsa ON gsa.aliquot_barcode = gs.aliquot_barcode
)
SELECT gc.aliquot_barcode,
gc.gene_symbol,
CASE
WHEN gc.wcr >= ss.del_thres AND gc.wcr <= ss.amp_thres THEN 0
WHEN gc.wcr < ss.hldel_thres THEN '-2'::integer
WHEN gc.wcr < ss.del_thres THEN '-1'::integer
WHEN gc.wcr > ss.hlamp_thres THEN 2
WHEN gc.wcr > ss.amp_thres THEN 1
ELSE NULL::integer
END AS hlvl_call
FROM gene_sample_call gc
LEFT JOIN seg_stats_optimized ss ON ss.aliquot_barcode = gc.aliquot_barcode
ORDER BY (
CASE
WHEN gc.wcr >= ss.del_thres AND gc.wcr <= ss.amp_thres THEN 0
WHEN gc.wcr < ss.hldel_thres THEN '-2'::integer
WHEN gc.wcr < ss.del_thres THEN '-1'::integer
WHEN gc.wcr > ss.hlamp_thres THEN 2
WHEN gc.wcr > ss.amp_thres THEN 1
ELSE NULL::integer
END);"
RB1 <- DBI::dbGetQuery(con, r) 

####################################################
# Investigate RB1 results
####################################################
#RB1 mutation
mut_RB1 <- mut_df %>% 
  filter(gene_symbol == "RB1") %>% 
  select(aliquot_barcode, gene_symbol, RB1_mutation = ssm2_pass_call) %>% 
  filter(aliquot_barcode %in% c(gold_set$tumor_barcode_a, gold_set$tumor_barcode_b))

#RB1 homozygous deletion
del_RB1 <- RB1 %>% 
  filter(aliquot_barcode %in% c(gold_set$tumor_barcode_a, gold_set$tumor_barcode_b)) %>% 
  mutate(final_cnv = ifelse(hlvl_call == -2, 1, 0)) %>% 
  select(-hlvl_call)

####################################################
# Investigate CDKN2A, CDK4, CDK6, CCND2 results
####################################################
# create dataframe with cnv tendency with expected hlvl call (-1 when deletion expected, 1 when amplification expected)
expected_cnv <- tribble(~ gene_symbol , ~ expected_tendency, 
                        "CCND2", 1, 
                        "CDK4", 1, 
                        "CDK6", 1, 
                        "CDKN2A", -1)

#make a cnv dataframe in which 0 indicates no alteration and 1 alteration (hlvl call = -2 or =2), blocked is labeled as 3
cnv <- cnv_by_gene %>% 
  left_join(tumor_blocklist, by = "aliquot_barcode") %>% # join with blocklist
  mutate(exclusion = trimws(cnv_exclusion, which = "both")) %>% # delete spaces before and after characters for cnv_exclusion column
  inner_join(expected_cnv, by = "gene_symbol") %>% # join with expected cnv list 
  mutate(final_cnv = case_when(exclusion == "block" ~  3,
                               gene_symbol == "CDKN2A" & expected_tendency == -1 & hlvl_call == -2 ~ 1, 
                               gene_symbol != "CDKN2A" & expected_tendency == 1 & hlvl_call == 2 ~ 1, 
                               TRUE ~ 0)) %>% 
  select(aliquot_barcode, gene_symbol, final_cnv) %>% # select columns of interest
  filter(aliquot_barcode %in% c(gold_set$tumor_barcode_a, gold_set$tumor_barcode_b)) %>% 
  merge(del_RB1, all = T) #add RB1 hom_del to this

###############################################################
################ Add subtype to gold-set ####################
#get annotation of Primary /Recurrence for Gold set
gold_select <- gold_set %>% 
  select(tumor_barcode_a, tumor_barcode_b)
list_all = c(gold_select$tumor_barcode_a, gold_select$tumor_barcode_b) %>% tibble() 
colnames(list_all) <- "aliquot_barcode"

#add tumor_class to tumor_barcode
list1 <- rep("P",length(1:222))
list2 <- rep("R", length(1:222))
list3 <- c(list1, list2)
gold_select_prim_rec <- add_column(list_all, sample_class = list3) %>% 
  mutate(case_barcode = substr(aliquot_barcode, 1,12)) %>% 
  left_join(subtypes, by = "case_barcode") %>% 
  select(aliquot_barcode, sample_class, idh_codel_subtype)

#########################################################
################ Combine CNV and SNV ####################
combined <- full_join(cnv, mut_RB1, by = c("aliquot_barcode", "gene_symbol")) %>%  #join snv and cnv dataframes
  distinct()
###################################################################
################ Add CNV and SNV to Gold Set ####################
gold_combined <- gold_select_prim_rec  %>% 
  left_join(combined, by = "aliquot_barcode") %>% 
  mutate(alteration = case_when(final_cnv == 1 | RB1_mutation == 1 ~ 1, 
                                           TRUE ~ 0)) %>% 
  select(-final_cnv, - RB1_mutation) %>% 
  as.data.frame() 

# further narrow down on (any) cell cycle alteration
gold_combined_cell_cycle <- gold_combined %>% 
  group_by(aliquot_barcode) %>% 
  mutate(cell_cycle_alteration = sum(alteration)) %>% 
  select(-gene_symbol, -alteration) %>% 
  distinct() %>% 
  mutate(case_barcode = substr(aliquot_barcode, 1, 12)) %>% 
  mutate(cell_cycle_subtype = ifelse(cell_cycle_alteration >= 1, "Cell_cycle_alt", "No_alt")) %>% 
  unite(idh_codel_subtype, cell_cycle_subtype, col = "final_subtype", sep = "_", remove = F)

## stats
sum(gold_combined_cell_cycle$sample_class == "P" & gold_combined_cell_cycle$cell_cycle_alteration >= 1) #primary: n=99 
sum(gold_combined_cell_cycle$sample_class == "R" & gold_combined_cell_cycle$cell_cycle_alteration >= 1) #recurrence: n=112

#########################################################
#Kaplan-Meier survival estimates for Cell_cycle status at P
#########################################################
#load in survival data and only choose primary tumors
gold_os_cell_cycle_P <- os %>% 
  select(case_barcode, case_overall_survival_mo, case_vital_status) %>% 
  inner_join(gold_combined_cell_cycle, by = "case_barcode") %>% 
  rename(time = case_overall_survival_mo, status = case_vital_status) %>%
  filter(sample_class == "P") %>% 
  mutate(end_status = ifelse(status == "dead", 2, 1))
  
#survival plot
ggsurvplot(survfit(Surv(time, end_status)
                   ~ final_subtype, data = gold_os_cell_cycle_P))

pairwise_survdiff(Surv(time, end_status) ~ final_subtype,
                  data = gold_os_cell_cycle_P)

###################################################################################
#Kaplan-Meier survival estimates for CDKN2A status at P for IDH mut noncodels
###################################################################################
gold_os_cell_cycle_P_idh_noncodel <- gold_os_cell_cycle_P %>% filter(idh_codel_subtype == "IDHmut-noncodel")

p2 <- ggsurvplot(survfit(Surv(time, end_status)
                         ~ final_subtype, data = gold_os_cell_cycle_P_idh_noncodel),
                 pval = TRUE, 
                 pval.coord = c(150, 0.5),
                 risk.table = T,
                 risk.table.col = "strata", # Change risk table color by groups
                 linetype = 1,
                 surv.median.line = "v", # Specify median survival
                 palette = c("#CD4F39","#27408B"),
                 legend.title = "Cell cycle status\nInitial tumor",
                 legend.labs = c("Alteration", "No alteration"),
                 legend = c(0.8,0.8), 
                 ylab = "Overall survival \n probability",
                 xlab = "Time (months)",
                 title = "IDHmut-noncodel"
)
summary(coxph(Surv(time, end_status)
              ~ final_subtype, data = gold_os_cell_cycle_P_idh_noncodel))

#########################################################
#Kaplan-Meier survival estimates for Cell_cycle status at R
#########################################################
gold_os_cell_cycle_R <- os %>% 
  select(case_barcode, case_overall_survival_mo, case_vital_status) %>% 
  inner_join(gold_combined_cell_cycle, by = "case_barcode") %>% 
  rename(time = case_overall_survival_mo, status = case_vital_status) %>%
  filter(sample_class == "R") %>% 
  mutate(end_status = ifelse(status == "dead", 2, 1))

#survival plot
ggsurvplot(survfit(Surv(time, end_status)
                   ~ final_subtype, data = gold_os_cell_cycle_R))

pairwise_survdiff(Surv(time, end_status) ~ final_subtype,
                  data = gold_os_cell_cycle_R)

# define list of samples with CDKN2A homozygous deletion
CDKN2A <- gold_combined %>% filter(gene_symbol == "CDKN2A" & alteration == 1) %>% select(aliquot_barcode)

####################################################
# Investigate CDKN2A results
####################################################
########################
#compare only CDKN2A (=p16), particular emphasis on hlvl_call = -2 (homozygous deletion)
p16_df <- cnv_by_gene %>%
  filter(gene_symbol == "CDKN2A") 

#join with gold set
p16_gold_df <- inner_join(gold_select_prim_rec, p16_df, by = "aliquot_barcode") %>% 
  mutate(case_barcode = substr(aliquot_barcode, 1, 12)) %>% 
  mutate(CDKN2A_subtype = ifelse(hlvl_call == -2, "CDKN2A_homdel", "CDKN2A_wt")) %>% 
  unite(idh_codel_subtype, CDKN2A_subtype, col = "final_subtype", sep = "_", remove = F)

## stats
sum(p16_gold_df$sample_class == "P" & p16_gold_df$hlvl_call == -2, na.rm = TRUE) #primary: n=53 
sum(p16_gold_df$sample_class == "R" & p16_gold_df$hlvl_call == -2, na.rm = TRUE) #recurrence: n=71


#########################################################
#Kaplan-Meier survival estimates for CDKN2A status at P
#########################################################
#load in survival data and only choose primary tumors
p16_gold_os_P <- os %>% 
  select(case_barcode, case_overall_survival_mo, case_vital_status) %>% 
  inner_join(p16_gold_df, by = "case_barcode") %>% 
  rename(time = case_overall_survival_mo, status = case_vital_status) %>%
  filter(sample_class == "P") %>% 
  mutate(end_status = ifelse(status == "dead", 2, 1))

#survival plot
ggsurvplot(survfit(Surv(time, end_status)
           ~ final_subtype, data = p16_gold_os_P))

pairwise_survdiff(Surv(time, end_status) ~ final_subtype,
                  data = p16_gold_os_P)

###################################################################################
#Kaplan-Meier survival estimates for CDKN2A status at P for IDH mut noncodels
###################################################################################
p16_gold_os_P_idh_noncodel <- p16_gold_os_P %>% filter(idh_codel_subtype == "IDHmut-noncodel")

#survival plot
p1 <- ggsurvplot(survfit(Surv(time, end_status)
                         ~ final_subtype, data = p16_gold_os_P_idh_noncodel),
           pval = TRUE, 
           pval.coord = c(150, 0.5),
           risk.table = T,
           risk.table.col = "strata", # Change risk table color by groups
           linetype = 1,
           surv.median.line = "v", # Specify median survival
           palette = c("#CD4F39","#27408B"),
           legend.title = "CDKN2A status\nInitial tumor",
           legend.labs = c("CDKN2A homdel", "CDKN2A WT"),
           legend = c(0.8,0.8), 
           ylab = "Overall survival \n probability",
           xlab = "Time (months)",
           title = "IDHmut-noncodel"
)

summary(coxph(Surv(time, end_status)
      ~ final_subtype, data = p16_gold_os_P_idh_noncodel))
#########################################################
#Kaplan-Meier survival estimates for CDKN2A status at R
#########################################################
#load in survival data and only choose primary tumors
p16_gold_os_R <- os %>% 
  select(case_barcode, case_overall_survival_mo, case_vital_status) %>% 
  inner_join(p16_gold_df, by = "case_barcode") %>% 
  rename(time = case_overall_survival_mo, status = case_vital_status) %>%
  filter(sample_class == "R") %>% 
  mutate(end_status = ifelse(status == "dead", 2, 1))

#survival plot

ggsurvplot(survfit(Surv(time, end_status)
                   ~ final_subtype, data = p16_gold_os_R))

pairwise_survdiff(Surv(time, end_status) ~ final_subtype,
                  data = p16_gold_os_R)

#########################################################
# Build 4 groups based on CDKN2A status 
# 1. No alterations at P and R (none)
# 2. Alteration at P and R (shared)
# 3. No alteration at P, but gained alteration at R (recurrence-only)
# 4. Alteration at P, but no alteration at R (primary-only)
#########################################################
# load CDKN2A data into gold set 
p16_gold_comp <- gold_set %>% 
  rename(aliquot_barcode = tumor_barcode_a) %>% 
  inner_join(p16_df, by = "aliquot_barcode") %>% 
  select(case_barcode, Primary = aliquot_barcode, P = hlvl_call, aliquot_barcode  = tumor_barcode_b) %>% 
  inner_join(p16_df, by = "aliquot_barcode") %>% 
  select(case_barcode, P, R = hlvl_call) %>% 
  mutate(CDKN2A_homdel = case_when(P != -2 & R != -2 ~ "None", 
                                   P == -2 & R == -2 ~ "Shared", 
                                   P != -2 & R == -2 ~ "Recurrence-only", 
                                   P == -2 & R != -2 ~ "Primary-only", 
                                   TRUE ~ "NA")) %>% 
  inner_join(subtypes, by = "case_barcode") %>% 
  unite(idh_codel_subtype, CDKN2A_homdel, sep = " ", col = "final_subtype", remove = F) %>% 
  inner_join(os, by = "case_barcode") %>% 
  select(case_barcode, final_subtype, status = case_vital_status, time = case_overall_survival_mo, idh_codel_subtype) %>% 
  mutate(end_status = ifelse(status == "dead", 2, 1))


#### test for each subtype seperately #####
# IDHwt
p16_gold_comp_IDHwt <- p16_gold_comp %>% filter(idh_codel_subtype == "IDHwt")
ggsurvplot(survfit(Surv(time, end_status)
                   ~ final_subtype, data = p16_gold_comp_IDHwt))
#pairwise test
pairwise_survdiff(Surv(time, end_status) ~ final_subtype,
                  data = p16_gold_comp_IDHwt)

# IDHmut-codel
p16_gold_comp_IDHmut_codel <- p16_gold_comp %>% filter(idh_codel_subtype == "IDHmut-codel")

ggsurvplot(survfit(Surv(time, end_status)
                   ~ final_subtype, data = p16_gold_comp_IDHmut_codel))
#pairwise test
pairwise_survdiff(Surv(time, end_status) ~ final_subtype,
                  data = p16_gold_comp_IDHmut_codel)

# IDHmut-noncodel
p16_gold_comp_IDHmut_noncodel <- p16_gold_comp %>% filter(idh_codel_subtype == "IDHmut-noncodel")

ggsurvplot(survfit(Surv(time, end_status)
                   ~ final_subtype, data = p16_gold_comp_IDHmut_noncodel))
#pairwise test
pairwise_survdiff(Surv(time, end_status) ~ final_subtype,
                  data = p16_gold_comp_IDHmut_noncodel) 

#######################################################################
# Correlation of Aneuploidy-change with CDKN2A status
#######################################################################
# Construct a table that provides subject-level information about clinical variables between two timepoints.
clinical_tumor_pairs_query = read_file("sql/clinical/clinical-tumor-pairs-db2.sql")
clinical_tumor_pairs <- dbGetQuery(con, clinical_tumor_pairs_query)

# Define aneuploidy pairs from the gold set.
aneuploidy_pairs <- dbGetQuery(con, "SELECT tumor_pair_barcode, case_barcode, tumor_barcode_a, tumor_barcode_b, a1.prop_aneuploidy AS aneuploidy_a, a1.aneuploidy_score::integer AS aneuploidy_score_a, 
                               a1.aneuploidy_amp_score::integer AS aneuploidy_amp_score_a, a1.aneuploidy_del_score::integer AS aneuploidy_del_score_a, 
                               a2.prop_aneuploidy AS aneuploidy_b, a2.aneuploidy_score::integer AS aneuploidy_score_b, a2.aneuploidy_amp_score::integer AS aneuploidy_amp_score_b, a2.aneuploidy_del_score::integer AS aneuploidy_del_score_b
                               FROM analysis.gold_set gs
                               LEFT JOIN analysis.gatk_aneuploidy a1 ON a1.aliquot_barcode = gs.tumor_barcode_a
                               LEFT JOIN analysis.gatk_aneuploidy a2 ON a2.aliquot_barcode = gs.tumor_barcode_b")

# Incorporate tables containing age, survival, and subtype.
aneuploidy_pairs_clin = aneuploidy_pairs %>% 
  left_join(clinical_tumor_pairs, by=c("tumor_pair_barcode", "case_barcode", "tumor_barcode_a", "tumor_barcode_b")) %>% 
  filter(tumor_pair_barcode %in% gold_set$tumor_pair_barcode) %>% 
  left_join(subtypes, by="case_barcode") %>% 
  left_join(os, by = "case_barcode")  

# Define aneuploidy pairs by IDH subtype.
aneuploidy_pairs_IDHwt = aneuploidy_pairs_clin %>% 
  filter(idh_codel_subtype == "IDHwt") %>% 
  mutate(aneuploidy_diff = aneuploidy_b-aneuploidy_a)
aneuploidy_pairs_IDH_codel = aneuploidy_pairs_clin %>% 
  filter(idh_codel_subtype == "IDHmut-codel") %>% 
  mutate(aneuploidy_diff = aneuploidy_b-aneuploidy_a)
aneuploidy_pairs_IDH_noncodel = aneuploidy_pairs_clin %>% 
  filter(idh_codel_subtype == "IDHmut-noncodel") %>% 
  mutate(aneuploidy_diff = aneuploidy_b-aneuploidy_a)

# include CDKN2A status with further information (hypermutation, surgical intervals)
aneuploidy_idh_noncodel = aneuploidy_pairs_IDH_noncodel %>% 
  mutate(post_recurrence_surv = case_overall_survival_mo-surgical_interval,
         acquired_aneuploidy = ifelse(aneuploidy_diff >= 0.12285, 1, 0),
         aneuploidy_levels = ntile(aneuploidy_diff, 3),
         patient_vital = ifelse(case_vital_status=="alive", 0, 1)) %>% 
  left_join(mutation_freq, by=c("tumor_barcode_b" = "aliquot_barcode")) %>% 
  mutate(hypermutation = ifelse(coverage_adj_mut_freq > 10, 1, 0)) %>% 
  inner_join(p16_gold_comp) %>% 
  dplyr::select(tumor_pair_barcode, tumor_barcode_a, tumor_barcode_b, aneuploidy_a, aneuploidy_b, aneuploidy_diff, acquired_aneuploidy, hypermutation, mutation_count, patient_vital, case_overall_survival_mo, post_recurrence_surv, final_subtype)   

### Association studies CDKN2A hom del vs. Acquired Aneuploidy
#1
# Are CDKN2A hom del tumors also the acquired_aneuploidy tumors? (high-level)
table(aneuploidy_idh_noncodel$final_subtype, aneuploidy_idh_noncodel$acquired_aneuploidy)
# Acquired aneupploidy and CDKN2A status seem to be associated (high-level)
fisher.test(table(aneuploidy_idh_noncodel$final_subtype, aneuploidy_idh_noncodel$acquired_aneuploidy))

#2
# Are CDKN2A hom del tumors also the acquired_aneuploidy tumors? (None vs any CDKN2A hom del)
aneuploidy_idh_noncodel_anyp16 <- aneuploidy_idh_noncodel %>% mutate(final_subtype = recode(final_subtype, 
                                                                                            "IDHmut-noncodel Primary-only" = "IDHmut-noncodel CDKN2A homdel", 
                                                                                            "IDHmut-noncodel Recurrence-only" = "IDHmut-noncodel CDKN2A homdel", 
                                                                                            "IDHmut-noncodel Shared" = "IDHmut-noncodel CDKN2A homdel"))

table(aneuploidy_idh_noncodel_anyp16$final_subtype, aneuploidy_idh_noncodel_anyp16$acquired_aneuploidy)
# Acquired aneupploidy and CDKN2A hom_del seem to be highly associated (None vs any CDKN2A hom del)
fisher.test(table(aneuploidy_idh_noncodel_anyp16$final_subtype, aneuploidy_idh_noncodel_anyp16$acquired_aneuploidy))

#3
# Are aqcuired CDKN2A hom del tumors also the acquired_aneuploidy tumors? (None vs CDKN2A hom del at recurrence only)
aneuploidy_idh_noncodel_recp16 <- aneuploidy_idh_noncodel %>% filter(final_subtype %in% c("IDHmut-noncodel Recurrence-only", "IDHmut-noncodel None"))

table(aneuploidy_idh_noncodel_recp16$final_subtype, aneuploidy_idh_noncodel_recp16$acquired_aneuploidy)
# Acquired aneupploidy and gained CDKN2A hom_del at recurrence seem to be highly associated (None vs CDKN2A hom del at recurrence)
fisher.test(table(aneuploidy_idh_noncodel_recp16$final_subtype, aneuploidy_idh_noncodel_recp16$acquired_aneuploidy))

# ---> Note that CDKN2A homozygous deletion is highly associated with acquired aneuploidy


### Association studies CDKN2A hom del vs. Hypermutation
#1
# Are aqcuired CDKN2A hom del tumors also the hypermutator tumors? (high-level)
table(aneuploidy_idh_noncodel$final_subtype, aneuploidy_idh_noncodel$hypermutation)
# Hypermutation and CDKN2A status seem not to be associated (high-level)
fisher.test(table(aneuploidy_idh_noncodel$final_subtype, aneuploidy_idh_noncodel$hypermutation))

#2
# Are aqcuired CDKN2A hom del tumors also the hypermutator tumors? (None vs any CDKN2A hom del)
table(aneuploidy_idh_noncodel_anyp16$final_subtype, aneuploidy_idh_noncodel_anyp16$hypermutation)
# Hypermutation and CDKN2A status seem not to be associated (None vs any CDKN2A hom del)
fisher.test(table(aneuploidy_idh_noncodel_anyp16$final_subtype, aneuploidy_idh_noncodel_anyp16$hypermutation))

#3
#  Are aqcuired CDKN2A hom del tumors also the hypermutator tumors (None vs CDKN2A hom del at recurrence only)
table(aneuploidy_idh_noncodel_recp16$final_subtype, aneuploidy_idh_noncodel_recp16$hypermutation)
# Hypermutation and CDKN2A status seem not to be associated  (None vs CDKN2A hom del at recurrence only)
fisher.test(table(aneuploidy_idh_noncodel_recp16$final_subtype, aneuploidy_idh_noncodel_recp16$hypermutation))

# ---> Note that CDKN2A homozygous deletion is not associated with hypermutation

###########################################################################################
# Overall survival analysis of tumors with/without CDKN2A hom_del at Recurrence
###########################################################################################
rec_p16_idh_noncodel <- aneuploidy_idh_noncodel %>% mutate(final_subtype = recode(final_subtype, 
                                                                          "IDHmut-noncodel None" = "IDHmut-noncodel CDKN2A wt", 
                                                                          "IDHmut-noncodel Primary-only" = "IDHmut-noncodel CDKN2A wt", 
                                                                          "IDHmut-noncodel Recurrence-only" = "IDHmut-noncodel CDKN2A homdel", 
                                                                          "IDHmut-noncodel Shared" = "IDHmut-noncodel CDKN2A homdel"))

p3 <- ggsurvplot(survfit(Surv(case_overall_survival_mo, patient_vital)
                         ~ final_subtype, data = rec_p16_idh_noncodel),
           pval = TRUE, 
           pval.coord = c(150, 0.5),
           risk.table = T,
           risk.table.col = "strata", # Change risk table color by groups
           linetype = 1,
           surv.median.line = "v", # Specify median survival
           palette = c("#CD4F39", "#27408B"),
           legend.title = "CDKN2A status\nRecurrence",
           legend.labs = c("CDKN2A homdel", "CDKN2A WT"),
           legend = c(0.8,0.8), 
           ylab = "Overall survival \n probability",
           xlab = "Time (months)",
           title = "IDHmut-noncodel"
)
summary(coxph(Surv(case_overall_survival_mo, patient_vital)
        ~ final_subtype, data = rec_p16_idh_noncodel))

#########################################################
# Build 4 groups based on cell cycle status 
# 1. No alterations at P and R (none)
# 2. Alteration at P and R (shared)
# 3. No alteration at P, but gained alteration at R (recurrence-only)
# 4. Alteration at P, but no alteration at R (primary-only)
#########################################################
# load cell_cycle data into gold set 
cell_cycle_gold_comp <- gold_set %>% 
  rename(aliquot_barcode = tumor_barcode_a) %>% 
  inner_join(gold_combined_cell_cycle, by = c("aliquot_barcode", "case_barcode")) %>% 
  select(case_barcode, Primary = aliquot_barcode, P = cell_cycle_alteration, aliquot_barcode  = tumor_barcode_b) %>% 
  inner_join(gold_combined_cell_cycle, by = c("aliquot_barcode", "case_barcode")) %>% 
  select(case_barcode, P, R = cell_cycle_alteration) %>% 
  mutate(CC_alteration = case_when(P == 0 & R == 0 ~ "None", 
                                   P != 0 & R != 0 ~ "Shared", 
                                   P == 0 & R != 0 ~ "Recurrence-only", 
                                   P != 0 & R == 0 ~ "Primary-only", 
                                   TRUE ~ "NA")) %>% 
  inner_join(subtypes, by = "case_barcode") %>% 
  unite(idh_codel_subtype, CC_alteration, sep = " ", col = "final_subtype", remove = F) %>% 
  inner_join(os, by = "case_barcode") %>% 
  select(case_barcode, final_subtype, status = case_vital_status, time = case_overall_survival_mo, idh_codel_subtype) %>% 
  mutate(end_status = ifelse(status == "dead", 2, 1))


#### test for each subtype seperately #####
# IDHwt
cell_cycle_gold_comp_IDHwt <- cell_cycle_gold_comp %>% filter(idh_codel_subtype == "IDHwt")

ggsurvplot(survfit(Surv(time, end_status)
                   ~ final_subtype, data = cell_cycle_gold_comp_IDHwt))
#pairwise test
pairwise_survdiff(Surv(time, end_status) ~ final_subtype,
                  data = cell_cycle_gold_comp_IDHwt) # no significant differences

# IDHmut-codel
cell_cycle_gold_comp_IDHmut_codel <- cell_cycle_gold_comp %>% filter(idh_codel_subtype == "IDHmut-codel")

ggsurvplot(survfit(Surv(time, end_status)
                   ~ final_subtype, data = cell_cycle_gold_comp_IDHmut_codel))
#pairwise test
pairwise_survdiff(Surv(time, end_status) ~ final_subtype,
                  data = cell_cycle_gold_comp_IDHmut_codel) # no significant differences

# IDHmut-noncodel
cell_cycle_gold_comp_IDHmut_noncodel <- cell_cycle_gold_comp %>% filter(idh_codel_subtype == "IDHmut-noncodel")

ggsurvplot(survfit(Surv(time, end_status)
                   ~ final_subtype, data = cell_cycle_gold_comp_IDHmut_noncodel))
#pairwise test
pairwise_survdiff(Surv(time, end_status) ~ final_subtype,
                  data = cell_cycle_gold_comp_IDHmut_noncodel) # no significant difference between groups

#further subdivision of cell cycle alteration in different groups has no effect
#######################################################################
# Correlation of Aneuploidy-change with Cell-cycle-alteration status (at Recurrence)
#######################################################################
aneuploidy_idh_noncodel_cell_cycle = aneuploidy_pairs_IDH_noncodel %>% 
  mutate(post_recurrence_surv = case_overall_survival_mo-surgical_interval,
         acquired_aneuploidy = ifelse(aneuploidy_diff >= 0.12285, 1, 0),
         aneuploidy_levels = ntile(aneuploidy_diff, 3),
         patient_vital = ifelse(case_vital_status=="alive", 0, 1)) %>% 
  left_join(mutation_freq, by=c("tumor_barcode_b" = "aliquot_barcode")) %>% 
  mutate(hypermutation = ifelse(coverage_adj_mut_freq > 10, 1, 0)) %>% 
  inner_join(cell_cycle_gold_comp) %>% 
  dplyr::select(tumor_pair_barcode, tumor_barcode_a, tumor_barcode_b, aneuploidy_a, aneuploidy_b, aneuploidy_diff, acquired_aneuploidy, hypermutation, mutation_count, patient_vital, case_overall_survival_mo, post_recurrence_surv, final_subtype)   

### Association studies cell cycle alteration vs. Acquired Aneuploidy
#1
# Are cell cycle altered tumors also the acquired_aneuploidy tumors? (high-level)
table(aneuploidy_idh_noncodel_cell_cycle$final_subtype, aneuploidy_idh_noncodel_cell_cycle$acquired_aneuploidy)
# Acquired aneupploidy and cell cycle status seem to be associated (high-level)
fisher.test(table(aneuploidy_idh_noncodel_cell_cycle$final_subtype, aneuploidy_idh_noncodel_cell_cycle$acquired_aneuploidy))

#2
# Are cell cycle altered tumors also the acquired_aneuploidy tumors? (None vs any cell cycle alt)
aneuploidy_idh_noncodel_any_cellcycle <- aneuploidy_idh_noncodel_cell_cycle %>% mutate(final_subtype = recode(final_subtype, 
                                                                                            "IDHmut-noncodel Primary-only" = "IDHmut-noncodel Cell_cycle_alteration", 
                                                                                            "IDHmut-noncodel Recurrence-only" = "IDHmut-noncodel Cell_cycle_alteration", 
                                                                                            "IDHmut-noncodel Shared" = "IDHmut-noncodel Cell_cycle_alteration"))

table(aneuploidy_idh_noncodel_any_cellcycle$final_subtype, aneuploidy_idh_noncodel_any_cellcycle$acquired_aneuploidy)
# Acquired aneupploidy and Cell cycle alteration seem to be highly associated (None vs any CDKN2A hom del)
fisher.test(table(aneuploidy_idh_noncodel_any_cellcycle$final_subtype, aneuploidy_idh_noncodel_any_cellcycle$acquired_aneuploidy))

#3
# Are aqcuired cell cycle alteration tumors also the acquired_aneuploidy tumors? (None vs cell cycle alt at recurrence only)
aneuploidy_idh_noncodel_rec_cellcycle <- aneuploidy_idh_noncodel_cell_cycle %>% filter(final_subtype %in% c("IDHmut-noncodel Recurrence-only", "IDHmut-noncodel None"))

table(aneuploidy_idh_noncodel_rec_cellcycle$final_subtype, aneuploidy_idh_noncodel_rec_cellcycle$acquired_aneuploidy)
# Acquired aneupploidy and gained cell cycle alteration at recurrence seem to be highly associated (None vs cell cycle alt at recurrence only)
fisher.test(table(aneuploidy_idh_noncodel_rec_cellcycle$final_subtype, aneuploidy_idh_noncodel_rec_cellcycle$acquired_aneuploidy))


# ---> Note that Cell cycle alteration is highly associated with acquired aneuploidy

### Association studies cell cycle alteration vs. Hypermutation
#1
# Are aqcuired cell cycle alt. tumors also the hypermutator tumors? (high-level)
table(aneuploidy_idh_noncodel_cell_cycle$final_subtype, aneuploidy_idh_noncodel_cell_cycle$hypermutation)
# Hypermutation and cell cycle alteration seem not to be associated (high-level)
fisher.test(table(aneuploidy_idh_noncodel_cell_cycle$final_subtype, aneuploidy_idh_noncodel_cell_cycle$hypermutation))

#2
# Are aqcuired cell cycle alt.  tumors also the hypermutator tumors? (None vs any cell cycle alteration)
table(aneuploidy_idh_noncodel_any_cellcycle$final_subtype, aneuploidy_idh_noncodel_any_cellcycle$hypermutation)
# Hypermutation and cell cycle alteration seem not to be associated (None vs any cell cycle alteration)
fisher.test(table(aneuploidy_idh_noncodel_any_cellcycle$final_subtype, aneuploidy_idh_noncodel_any_cellcycle$hypermutation))

#3
# Are aqcuired cell cycle alt. tumors also the hypermutator tumors? (None vs cell cycle alteration at recurrence only)
table(aneuploidy_idh_noncodel_rec_cellcycle$final_subtype, aneuploidy_idh_noncodel_rec_cellcycle$hypermutation)
# Hypermutation and cell cycle alteration seem not to be associated  (None vs cell cycle alteration at recurrence only)
fisher.test(table(aneuploidy_idh_noncodel_rec_cellcycle$final_subtype, aneuploidy_idh_noncodel_rec_cellcycle$hypermutation))

# ---> Note that cell cycle alteration in general is not associated with hypermutation

###########################################################################################
# Post-recurrence survival analysis of tumors with/without cell cycle alteration at Recurrence
###########################################################################################
rec_cell_cycle_idh_noncodel <- aneuploidy_idh_noncodel_cell_cycle %>% mutate(final_subtype = recode(final_subtype, 
                                                                                       "IDHmut-noncodel None" = "IDHmut-noncodel No alteration", 
                                                                                       "IDHmut-noncodel Primary-only" = "IDHmut-noncodel No alteration", 
                                                                                       "IDHmut-noncodel Recurrence-only" = "IDHmut-noncodel Cell cycle alteration", 
                                                                                       "IDHmut-noncodel Shared" = "IDHmut-noncodel Cell cycle alteration"))

p4 <- ggsurvplot(survfit(Surv(case_overall_survival_mo, patient_vital)
                         ~ final_subtype, data = rec_cell_cycle_idh_noncodel),
           pval = TRUE, 
           pval.coord = c(150, 0.5),
           risk.table = T,
           risk.table.col = "strata", # Change risk table color by groups
           linetype = 1,
           surv.median.line = "v", # Specify median survival
           palette = c("#CD4F39","#27408B"),
           legend.title = "Cell cycle status\nRecurrence",
           legend.labs = c("Alteration", "No alteration"),
           legend = c(0.8,0.8), 
           ylab = "Overall survival \n  probability",
           xlab = "Time (months)",
           title = "IDHmut-noncodel"
)
summary(coxph(Surv(case_overall_survival_mo, patient_vital)
              ~ final_subtype, data = rec_cell_cycle_idh_noncodel))

#################################################################################
#################################################################################
#Aneuploidy value comparison and highlighting samples with cell cycle alterations
#################################################################################
#################################################################################
# combine aneuploidy values with cell cycle information
cell_cycle <- gold_combined %>%  unite("class_gene_symbol", sample_class, gene_symbol) %>% mutate(case_barcode = substr(aliquot_barcode, 1, 12)) %>%  select(-aliquot_barcode) %>%
  spread(key = class_gene_symbol, value = alteration) %>% mutate(P_cc = P_CDKN2A + P_CCND2 + P_CDK4 + P_CDK6 + P_RB1, 
                                                                 R_cc = R_CDKN2A + R_CCND2 + R_CDK4 + R_CDK6 + R_RB1, 
                                                                 P = case_when(P_CDKN2A == 1 & P_cc >= 1 ~ "cdkn2a", P_CCND2 == 1 | P_CDK4 == 1 | P_CDK6 == 1 | P_RB1 == 1 & P_CDKN2A == 0 ~ "other_cc", TRUE ~ "none"), 
                                                                 R = case_when(R_CDKN2A == 1 & R_cc >= 1 ~ "cdkn2a", R_CCND2 == 1 | R_CDK4 == 1 | R_CDK6 == 1 | R_RB1 == 1 & R_CDKN2A == 0 ~ "other_cc", TRUE ~ "none"), 
                                                                 acquired = case_when(P == "none" & R == "cdkn2a" ~ "cc", P == "none" & R == "other_cc" ~ "cc", P %in% c("cdkn2a", "other_cc") & R %in% c("cdkn2a", "other_cc") ~ "shared", TRUE ~ "none"))

aneuploidy_cell_cycle <- aneuploidy_pairs %>% left_join(cell_cycle, by = "case_barcode")

# significance testing
aneuploidy_cell_cycle_noncodel <- aneuploidy_cell_cycle %>% filter(idh_codel_subtype == "IDHmut-noncodel")
wilcox_none <- aneuploidy_cell_cycle_noncodel %>% filter(acquired == "none")
wilcox_cc <- aneuploidy_cell_cycle_noncodel %>% filter(acquired == "cc")
wilcox_shared <- aneuploidy_cell_cycle_noncodel %>% filter(acquired == "shared")

wilcox.test(aneuploidy_cell_cycle_noncodel$aneuploidy_a, aneuploidy_cell_cycle_noncodel$aneuploidy_b, paired = T)
wilcox.test(wilcox_none$aneuploidy_a, wilcox_none$aneuploidy_b, paired = T)
wilcox.test(wilcox_cc$aneuploidy_a, wilcox_cc$aneuploidy_b, paired = T)
wilcox.test(wilcox_shared$aneuploidy_a, wilcox_shared$aneuploidy_b, paired = T)

aneuploidy_cell_cycle_plot = aneuploidy_cell_cycle %>% 
  gather(sample, aneuploidy, c(aneuploidy_a, aneuploidy_b)) %>% 
  mutate(sample = recode(sample,  "aneuploidy_a" = "Initial",  "aneuploidy_b" = "Recurrence"))


ladder <- aneuploidy_cell_cycle_plot %>% filter(idh_codel_subtype == "IDHmut-noncodel") %>% ggplot(aes(x = sample, y = aneuploidy, group = case_barcode, color = acquired)) +
  geom_line(linetype="solid", size=1, alpha = 0.3) + ylab("Aneuploidy value") + xlab("") + ggtitle("IDHmut-noncodel") + theme_bw() +
  geom_point(size=2) +
  annotate("text", x=1.5, y=1.03, label= "P = 1.41e-06", fontface= "italic", size = 5) + 
  scale_color_manual(values=c("cc" = "#FF7100", "none" = "black", "shared" = "#008EFF"), 
                     name="Cell cycle alteration",
                     labels=c("cc" = "Newly acquired", "shared" = "Shared", "none" = "None"))
ladder

###########################
###########################
#Exporting plots
###########################
###########################
### Enumarate CDKN2A status at primary vs recurrence
p16gold_df_P <- p16_gold_df %>% filter(sample_class == "P")
p16gold_df_R <- p16_gold_df %>% filter(sample_class == "R")

table(p16gold_df_P$idh_codel_subtype, p16gold_df_P$CDKN2A_subtype)
table(p16gold_df_R$idh_codel_subtype, p16gold_df_R$CDKN2A_subtype)


# IDH mut noncodel
splots <- list()
splots[[1]] <- p1
splots[[2]] <- p2

pdf(file = "/Users/c-kocake/Box Sync/Projects/GLASS/Cell_cycle/idhmut-noncodel_initial.pdf", height = 5, width = 12, bg = "transparent", useDingbats = FALSE)
arrange_ggsurvplots(splots, print = TRUE,
                    ncol = 2, nrow = 1, risk.table.height = 0.3)
dev.off()

splots_2 <- list()
splots_2[[1]] <- p3
splots_2[[2]] <- p4

pdf(file = "/Users/c-kocake/Box Sync/Projects/GLASS/Cell_cycle/idhmut-noncodel_recurrence.pdf", height = 5, width = 12, bg = "transparent", useDingbats = FALSE)
arrange_ggsurvplots(splots_2, print = TRUE,
                    ncol = 2, nrow = 1, risk.table.height = 0.3)
dev.off()

pdf(file = "/Users/c-kocake/Box Sync/Projects/GLASS/Cell_cycle/cell_cycle_aneuploidy.pdf", height = 5, width = 5, bg = "transparent", useDingbats = FALSE)
ladder
dev.off()

### END ###