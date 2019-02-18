##################################################
# Determine the mutational frequencies at each time point.
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
library(survival)
library(ggExtra)
library(EnvStats)
library(cowplot)
##################################################
# Establish connection with database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

# Essential tables to load:
clinal_tumor_pairs = dbReadTable(con,  Id(schema="clinical", table="clinical_by_tumor_pair"))  
tumor_mut_compare = dbReadTable(con,  Id(schema="analysis",table="tumor_mut_comparison_anno"))
cases = dbReadTable(con,  Id(schema="clinical", table="cases"))
silver_set = dbGetQuery(con, "SELECT * FROM analysis.silver_set")
titan_param = dbReadTable(con,  Id(schema="analysis", table="titan_params"))

# Analyze only those tumors included in the silver set.
mut_freq_silver = silver_set %>% 
  left_join(cases, by="case_barcode") %>% 
  inner_join(clinal_tumor_pairs, by="tumor_pair_barcode") %>% 
  inner_join(tumor_mut_compare, by="tumor_pair_barcode") 

# Combine with pairs file to access "tumor_barcode".
titan_info = titan_param %>% 
  left_join(pairs, by="pair_barcode")
titan_purity = titan_info %>% 
  select(tumor_barcode, purity)

mut_freq_silver_purity =  %>% 
  
  
# Stacked barplot for mutational frequency.
mut_freq_bar = mut_freq_silver %>% 
  filter(hypermutator_status==0) %>% 
  mutate(primary_only = (count_a-intersection_ab)/union_ab,
         recurrence_only =  (count_b-intersection_ab)/union_ab,
         shared = intersection_ab/union_ab) %>% 
  select(primary_only, recurrence_only, shared, tumor_pair_barcode, union_ab, idh_codel_subtype) %>% 
  gather(mutation_type, mutation_percent, c(primary_only, recurrence_only, shared), -tumor_pair_barcode, -union_ab, -idh_codel_subtype)
mut_freq_bar$mutation_type = factor(mut_freq_bar$mutation_type, levels = c("recurrence_only", "primary_only", "shared"))

# Total mutational frequency. 
top_plot = ggplot(mut_freq_bar, aes(x =reorder(tumor_pair_barcode, as.numeric(union_ab)), y=as.numeric(log10(union_ab)))) + geom_bar(stat="identity") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_text(size=rel(0.6))) +
  xlab("") + ylab("Log10(mutation union)") + facet_grid(~idh_codel_subtype, scales="free")
# Stacked barplot.
bot_plot = ggplot(mut_freq_bar, aes(x =reorder(tumor_pair_barcode, as.numeric(union_ab)), y=as.numeric(mutation_percent), fill=mutation_type)) + geom_bar(stat="identity") + 
  scale_fill_manual(values=c("#2FB3CA", "#CA2F66", "#CA932F")) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_text(size=rel(0.9))) +
  xlab("") + ylab("% SNVs and Indels") + guides(fill=guide_legend("Mutation")) + facet_grid(~idh_codel_subtype, scales="free") 

# Combine the plot to incorporate both mutation burden and proportions.
plot_grid(top_plot, bot_plot, align = "v", axis= "rl" , nrow = 2, rel_heights = c(1/4, 3/4), label_size = 12)
# Some pretty large differences across subtypes.
mut_freq_bar %>% 
  group_by(mutation_type, idh_codel_subtype) %>% 
  summarise(proportion_shared = median(mutation_percent))
# Are there any statistical differences.
prop_shared = mut_freq_bar %>% 
  filter(mutation_type=="shared")
# Test between subtypes for proportions of shared mutations.
kruskal.test(prop_shared$mutation_percent~factor(prop_shared$idh_codel_subtype))

# Restrict analyses initially to the non-hypermutators:
mut_freq_silver_non_hyper = mut_freq_silver %>% 
  filter(hypermutator_status==0) %>% 
  gather(sample_type, mut_freq, c(mf_a, mf_b)) %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, "IDHmut_codel" = "IDHmut codel", "IDHmut_noncodel" = "IDHmut noncodel", "IDHwt_noncodel" = "IDHwt"),
         sample_type = recode(sample_type, "mf_a" = "Primary", "mf_b" = "Recurrence"))



# Create specific colors for primary/recurrence designation.
ggplot(mut_freq_silver_non_hyper, aes(x=sample_type, y=mut_freq, fill=sample_type)) +
  geom_boxplot() + theme_bw() + xlab("") + ylab("") + labs("Tumor") + scale_fill_manual(values=c("#CA6C18", "#005496")) +
  facet_wrap(~idh_codel_subtype)

# Analysis only on the IDHwt tumors:
mut_freq_silver_non_hyper_IDHwt = mut_freq_silver %>% 
  filter(hypermutator_status==0)  %>% 
  filter(idh_codel_subtype == "IDHwt_noncodel")
mut_freq_silver_non_hyper_IDHwt %>% 
  summarise(med_mf_p = median(mf_a),
            med_mf_r = median(mf_b))
wilcox.test(mut_freq_silver_non_hyper_IDHwt$mf_a, mut_freq_silver_non_hyper_IDHwt$mf_b, paired = T)

# Only the IDHmut noncodel tumors:
mut_freq_silver_non_hyper_IDHmut_non = mut_freq_silver %>% 
  filter(hypermutator_status==0) %>% 
  filter(idh_codel_subtype == "IDHmut_noncodel")
mut_freq_silver_non_hyper_IDHmut_non %>% 
  summarise(med_mf_p = median(mf_a),
            med_mf_r = median(mf_b))
wilcox.test(mut_freq_silver_non_hyper_IDHmut_non$mf_a, mut_freq_silver_non_hyper_IDHmut_non$mf_b, paired = T)

# Only the IDHmut codel tumors:
mut_freq_silver_non_hyper_IDHmut_codel = mut_freq_silver %>% 
  filter(hypermutator_status==0) %>% 
  filter(idh_codel_subtype == "IDHmut_codel")
mut_freq_silver_non_hyper_IDHmut_codel %>% 
  summarise(med_mf_p = median(mf_a),
            med_mf_r = median(mf_b))
wilcox.test(mut_freq_silver_non_hyper_IDHmut_codel$mf_a, mut_freq_silver_non_hyper_IDHmut_codel$mf_b, paired = T)


# How many samples had a higher mutational frequency after recurrence?
mut_freq_silver_non_hyper = mut_freq_silver %>% 
  filter(hypermutator_status==0) 
sum(mut_freq_silver_non_hyper$delta_mf > 0)/dim(mut_freq_silver_non_hyper)[1]

#### Same approach, but present results for hypermutators:
# Restrict analyses initially to the non-hypermutators:
mut_freq_silver_hyper = mut_freq_silver %>% 
  filter(hypermutator_status==1) %>% 
  gather(sample_type, mut_freq, c(mf_a, mf_b)) %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, "IDHmut_codel" = "IDHmut codel", "IDHmut_noncodel" = "IDHmut noncodel", "IDHwt_noncodel" = "IDHwt"),
         sample_type = recode(sample_type, "mf_a" = "Primary", "mf_b" = "Recurrence"))
# Create specific colors for primary/recurrence designation.
ggplot(mut_freq_silver_hyper, aes(x=sample_type, y=mut_freq, fill=sample_type)) +
  geom_boxplot() + theme_bw() + xlab("") + ylab("") + labs("Tumor") + scale_fill_manual(values=c("#CA6C18", "#005496")) +
  facet_wrap(~idh_codel_subtype)

# Analysis only on the IDHwt tumors:
mut_freq_silver_hyper_IDHwt = mut_freq_silver %>% 
  filter(hypermutator_status==1)  %>% 
  filter(idh_codel_subtype == "IDHwt_noncodel")
mut_freq_silver_hyper_IDHwt %>% 
  summarise(med_mf_p = median(mf_a),
            med_mf_r = median(mf_b))

# Only the IDHmut noncodel tumors:
mut_freq_silver_hyper_IDHmut_non = mut_freq_silver %>% 
  filter(hypermutator_status==1) %>% 
  filter(idh_codel_subtype == "IDHmut_noncodel")
mut_freq_silver_hyper_IDHmut_non %>% 
  summarise(med_mf_p = median(mf_a),
            med_mf_r = median(mf_b))

# Only the IDHmut codel tumors:
mut_freq_silver_hyper_IDHmut_codel = mut_freq_silver %>% 
  filter(hypermutator_status==1) %>% 
  filter(idh_codel_subtype == "IDHmut_codel")
mut_freq_silver_hyper_IDHmut_codel %>% 
  summarise(med_mf_p = median(mf_a),
            med_mf_r = median(mf_b))

# Possible to include mutiple-sector data?
tumor_mut_multi = tumor_mut_compare %>% 
  filter(comparison_type!="longitudinal") %>% 
  gather(sample_type, mut_freq, c(mf_a, mf_b)) 

# Quick visualization of the multipl-sector data mutational frequencies.
ggplot(tumor_mut_multi, aes(x=case_barcode, y=mut_freq)) +
  geom_point() + theme_bw() + xlab("") + ylab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 10) +
  facet_wrap(~idh_codel_subtype, scales = "free") 
