##################################################
# Determine whether purity is associated with mutational frequency.
# Updated: 2019.05.17
# Author: Kevin J.
##################################################

# Working directory for this analysis in the GLASS-analysis project. 
mybasedir = "/Users/johnsk/Documents/Life-History/GLASS-WG/" 
setwd(mybasedir)

##################################################
# Necessary packages:
library(tidyverse)
library(DBI)
library(openxlsx)
#library(ggpubr)

##################################################
# Establish connection with database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB2")

# Load in the necessary information regarding purity
mutation_freq = dbReadTable(con,  Id(schema="analysis", table="mut_freq"))
silver_set = dbReadTable(con,  Id(schema="analysis", table="silver_set"))
gold_set = dbReadTable(con,  Id(schema="analysis", table="gold_set"))
blocklist = dbReadTable(con,  Id(schema="analysis", table="blocklist"))
ssm2_count = dbGetQuery(con, "SELECT * FROM variants.ssm2_count")
pairs = dbReadTable(con,  Id(schema="analysis", table="pairs"))
subtypes = dbReadTable(con,  Id(schema="clinical", table="subtypes"))
seqz_params = dbReadTable(con,  Id(schema="variants", table="seqz_params"))
titan_params = dbReadTable(con,  Id(schema="variants", table="titan_params"))

# Get the new clinical_tumor_pairs table.
clinical_tumor_pairs_query = read_file("sql/clinical/clinical-tumor-pairs-db2.sql")
clinical_tumor_pairs <- dbGetQuery(con, clinical_tumor_pairs_query)
clinical_tumor_gold = clinical_tumor_pairs %>%  
  filter(tumor_pair_barcode %in% gold_set$tumor_pair_barcode)

# Results from anti-join are dependent upon the order of the objects. Notice not all samples have purity estimates for each method.
anti_join(seqz_params, titan_params, by="pair_barcode")
anti_join(titan_params, seqz_params, by="pair_barcode")

# Currently, we use the gold_set for all major analyses.
gold_set_aliquots = gold_set %>% 
  gather(sample_type, aliquot_barcode, c(tumor_barcode_a, tumor_barcode_b)) %>% 
  arrange(case_barcode)

# Restrict to gold_set purity estimates.
seqz_purity = seqz_params %>% 
  inner_join(pairs, by= "pair_barcode") %>% 
  inner_join(gold_set_aliquots, by=c("tumor_barcode"="aliquot_barcode")) %>% 
  inner_join(mutation_freq, by=c("tumor_barcode"="aliquot_barcode")) %>% 
  mutate(sample_barcode = substr(pair_barcode, 1, 15),
         seq_type = substr(pair_barcode, 27, 29)) %>% 
  left_join(subtypes, by="case_barcode") 

# Correlation between tumor purity and mutational frequency in general. A weak positive association.
cor.test(seqz_purity$cellularity, seqz_purity$coverage_adj_mut_freq, method="spearman")

# Visualize the relationship between Sequenza tumor purity and mutational frequency.
pdf(file = "/Users/johnsk/Documents/sequenza-purity-mf-gold.pdf", height = 5, width = 7, bg = "transparent", useDingbats = FALSE)
ggplot(seqz_purity, aes(x=cellularity, y=log10(coverage_adj_mut_freq), color=idh_codel_subtype)) + geom_point() + theme_bw() + 
  ylab("log10(mutation freq.)") + xlab("Sequenza - cellularity") +
  annotate("text", x = 0.5, y=-1, label = "n = 442, rho=0.18, p=2.1E-04") +  geom_smooth(method='lm') + scale_color_discrete(name = "Glioma subtype")
dev.off()

# Determine the relationship between tumor purity and mutational freq. adjusted in the primary tumors.
seqz_purity_primary = seqz_purity %>% 
  filter(grepl("-TP-", pair_barcode))
seqz_purity_primary_codel = seqz_purity_primary %>% 
  filter(idh_codel_subtype == "IDHmut-codel")
seqz_purity_primary_noncodel = seqz_purity_primary %>% 
  filter(idh_codel_subtype == "IDHmut-noncodel")
seqz_purity_primary_wt = seqz_purity_primary %>% 
  filter(idh_codel_subtype == "IDHwt")

# Independently apply each statistical test.
cor.test(seqz_purity_primary_codel$cellularity, seqz_purity_primary_codel$coverage_adj_mut_freq, method="spearman")
cor.test(seqz_purity_primary_noncodel$cellularity, seqz_purity_primary_noncodel$coverage_adj_mut_freq, method="spearman")
# It's only strongly related with purity in the primary IDHwt.
cor.test(seqz_purity_primary_wt$cellularity, seqz_purity_primary_wt$coverage_adj_mut_freq, method="spearman")

# Breakdown the subtypes at recurrence.
seqz_purity_recurrences = seqz_purity %>% 
  filter(!grepl("-TP-", pair_barcode))
seqz_purity_recurrences_codel = seqz_purity_recurrences %>% 
  filter(idh_codel_subtype == "IDHmut-codel")
seqz_purity_recurrences_noncodel = seqz_purity_recurrences %>% 
  filter(idh_codel_subtype == "IDHmut-noncodel")
seqz_purity_recurrences_wt = seqz_purity_recurrences %>% 
  filter(idh_codel_subtype == "IDHwt")

# Apply to each recurrent set.
cor.test(seqz_purity_recurrences_codel$cellularity, seqz_purity_recurrences_codel$coverage_adj_mut_freq, method="spearman")
# Purity is more strongly associated with mut. freq. in the non-codels and IDHwt.
cor.test(seqz_purity_recurrences_noncodel$cellularity, seqz_purity_recurrences_noncodel$coverage_adj_mut_freq, method="spearman")
cor.test(seqz_purity_recurrences_wt$cellularity, seqz_purity_recurrences_wt$coverage_adj_mut_freq, method="spearman")

# Sequenza compare:
seqz_compare = seqz_purity %>% 
  mutate(sample_type = recode(sample_type,  "tumor_barcode_a" = "Primary", "tumor_barcode_b" = "Recurrence")) %>% 
  select(tumor_pair_barcode, case_barcode, sample_type, cellularity) %>% 
  spread(sample_type, cellularity) %>% 
  left_join(subtypes, by="case_barcode")

# Tests for significant differences in Sequenza tumor purity estimates.
seqz_compare_codel = seqz_compare %>% filter(idh_codel_subtype =="IDHmut-codel")
seqz_compare_noncodel = seqz_compare %>% filter(idh_codel_subtype =="IDHmut-noncodel")
seqz_compare_IDHwt  = seqz_compare %>% filter(idh_codel_subtype =="IDHwt")
wilcox.test(seqz_compare_codel$Primary, seqz_compare_codel$Recurrence, paired = TRUE)
# A significant difference in the primary vs. recurrence for IDHmut-noncodel.
wilcox.test(seqz_compare_noncodel$Primary, seqz_compare_noncodel$Recurrence, paired = TRUE)
wilcox.test(seqz_compare_IDHwt$Primary, seqz_compare_IDHwt$Recurrence, paired = TRUE)

# Sequenza has broader estimates of purity.
seqz_plot_purity = seqz_compare %>% 
  gather(sample, cellularity, c(Primary, Recurrence))

# Visualize the primary/recurrent differences for the GLASS Sequenza estimates.
ggplot(seqz_plot_purity, aes(x=idh_codel_subtype, y=cellularity, fill=sample)) + 
  geom_boxplot()  + ylab("Sequenza Purity") + xlab("") + theme_bw() + scale_fill_manual(values = c("Primary" = "#a6611a", "Recurrence" = "#018571")) + labs(fill="Sample Type")


########################
# TITAN purity
########################
titan_purity = titan_params %>% 
  inner_join(pairs, by= "pair_barcode") %>% 
  inner_join(gold_set_aliquots, by=c("tumor_barcode"="aliquot_barcode")) %>% 
  inner_join(mutation_freq, by=c("tumor_barcode"="aliquot_barcode")) %>% 
  mutate(sample_barcode = substr(pair_barcode, 1, 15),
         seq_type = substr(pair_barcode, 27, 29)) %>% 
  left_join(subtypes, by="case_barcode") 

# What's the spread of the purity estimates in GLASS.
summary(titan_purity$purity)

# Correlation between tumor purity and mutational frequency in general, not significant.
cor.test(titan_purity$purity, titan_purity$coverage_adj_mut_freq, method="spearman")

# Visualize the relationship between tumor purity and mutational frequency.
pdf(file = "/Users/johnsk/Documents/titan-purity-mf-goldset.pdf", height = 5, width = 7, bg = "transparent", useDingbats = FALSE)
ggplot(titan_purity, aes(x=purity, y=log10(coverage_adj_mut_freq), color=idh_codel_subtype)) + geom_point() + theme_bw() + 
  ylab("log10(mutation freq.)") + xlab("TITAN - purity") + xlim(0, 1) +
  annotate("text", x = 0.5, y=-1, label = "n = 442, rho = -0.01, p = 0.85", size =5) +  geom_smooth(method='lm') + scale_color_discrete(name = "Glioma subtype")
dev.off()

# Test whether this relationship persists in only the non-hypermutator group.
titan_purity_hypermut_status = titan_purity %>% 
  left_join(clinical_tumor_gold, by=c("tumor_barcode"="tumor_barcode_b")) 

# Restrict to non-hypermutator samples.
titan_purity_nonhyper = titan_purity_hypermut_status %>% 
  filter(hypermutator_status %in% c(0, NA))
# Still not significant.
cor.test(titan_purity_nonhyper$purity, titan_purity_nonhyper$coverage_adj_mut_freq, method="spearman")

# Determine the relationship between tumor purity and mutational freq. adjusted 
titan_purity_primary = titan_purity %>% 
  filter(grepl("-TP-", pair_barcode))
titan_purity_primary_codel = titan_purity_primary %>% 
  filter(idh_codel_subtype == "IDHmut-codel")
titan_purity_primary_noncodel = titan_purity_primary %>% 
  filter(idh_codel_subtype == "IDHmut-noncodel")
titan_purity_primary_wt = titan_purity_primary %>% 
  filter(idh_codel_subtype == "IDHwt")

# For the primary tumors apply correlation tests.
cor.test(titan_purity_primary_codel$purity, titan_purity_primary_codel$coverage_adj_mut_freq, method="spearman")
cor.test(titan_purity_primary_noncodel$purity, titan_purity_primary_noncodel$coverage_adj_mut_freq, method="spearman")
cor.test(titan_purity_primary_wt$purity, titan_purity_primary_wt$coverage_adj_mut_freq, method="spearman")

# Now just the recurrences.
titan_purity_recurrences = titan_purity %>% 
  filter(!grepl("-TP-", pair_barcode))
titan_purity_recurrences_codel = titan_purity_recurrences %>% 
  filter(idh_codel_subtype == "IDHmut-codel")
titan_purity_recurrences_noncodel = titan_purity_recurrences %>% 
  filter(idh_codel_subtype == "IDHmut-noncodel")
titan_purity_recurrences_wt = titan_purity_recurrences %>% 
  filter(idh_codel_subtype == "IDHwt")

# Is there a difference for either initial or recurrent tumors across time points.
kruskal.test(titan_purity_primary$purity, as.factor(titan_purity_primary$idh_codel_subtype))
kruskal.test(titan_purity_recurrences$purity, as.factor(titan_purity_recurrences$idh_codel_subtype))

# Recurrence-specific purity and mutational frequency. No TITAN-specific association.
cor.test(titan_purity_recurrences_codel$purity, titan_purity_recurrences_codel$coverage_adj_mut_freq, method="spearman")
cor.test(titan_purity_recurrences_noncodel$purity, titan_purity_recurrences_noncodel$coverage_adj_mut_freq, method="spearman")
cor.test(titan_purity_recurrences_wt$purity, titan_purity_recurrences_wt$coverage_adj_mut_freq, method="spearman")

# How many gold_set samples are above 50% purity?
sum(titan_purity$purity>0.5)/442

# Create tumor_pairs.
titan_purity_pairs = titan_purity %>% 
  select(case_barcode, sample_type, purity, idh_codel_subtype) %>% 
  # group_by(case_barcode) %>% 
  spread(sample_type, purity)

# Define TITAN gold_set by IDH subtype.
titan_pairs_IDHwt = titan_purity_pairs %>% 
  filter(idh_codel_subtype == "IDHwt")
titan_pairs_IDH_codel = titan_purity_pairs %>% 
  filter(idh_codel_subtype == "IDHmut-codel")
titan_pairs_IDH_noncodel = titan_purity_pairs %>% 
  filter(idh_codel_subtype == "IDHmut-noncodel")

# Test purity for statistical difference.
wilcox.test(titan_pairs_IDH_codel$tumor_barcode_a, titan_pairs_IDH_codel$tumor_barcode_b, paired = TRUE)
wilcox.test(titan_pairs_IDH_noncodel$tumor_barcode_a, titan_pairs_IDH_noncodel$tumor_barcode_b, paired = TRUE)
wilcox.test(titan_pairs_IDHwt$tumor_barcode_a, titan_pairs_IDHwt$tumor_barcode_b, paired = TRUE)

# Visualize the primary/recurrent differences for the GLASS TITAN estimates.
pdf(file = "/Users/johnsk/Documents/titan-purity-goldset-pairs.pdf", height = 5, width = 7, bg = "transparent", useDingbats = FALSE)
titan_purity = titan_purity %>% mutate(sample = recode(sample_type, "tumor_barcode_a" = "Initial", "tumor_barcode_b" = "Recurrence"))
ggplot(titan_purity, aes(x=idh_codel_subtype, y=purity, fill=sample)) + 
  geom_boxplot()  + ylab("TITAN Purity") + xlab("") + theme_bw() + scale_fill_manual(values = c("Recurrence" = "#a6611a", "Initial" = "#018571")) + labs(fill="Sample Type") +
  annotate("text", x = 1, y= 1.1, label = "n = 25 pairs\np = 0.50", size = 3) + annotate("text", x = 2, y= 1.1, label = "n = 63 pairs\np = 0.37", size = 3) +
  annotate("text", x = 3, y= 1.1, label = "n = 132 pairs\np = 0.67", size = 3) 
dev.off()



#############################################
# Compare available TCGA samples with Taylor et al
# ABSOLUTE array-based purity classification
#############################################
# Cancer Cell data:
taylor_data = readWorkbook("data/published/taylor-cancer-cell-2018.xlsx", sheet = 1, rowNames = F, colNames = TRUE)

# Extract only the iniformation for GBM and LGG. Recode some of their metadata to fit that of GLASS.
glioma_dat = taylor_data %>% 
  filter(Type%in%c("GBM", "LGG")) %>% 
  mutate(sample_type_num = substr(Sample, 14, 15), 
         sample_type = recode(sample_type_num, "01" = "TP", "02" = "R1"),
         sample_barcode = paste(substr(Sample, 1, 12), sample_type, sep ="-")) %>% 
  select(Sample, sample_type, sample_barcode, ABSOLUTE_purity = Purity)

# Examine how the purity data determined by ABSOLUTE match with GLASS estimates of purity for shared samples.
glass_taylor_seqz_purity = seqz_purity %>%
  inner_join(glioma_dat, by="sample_barcode")
cor.test(glass_taylor_seqz_purity$cellularity, glass_taylor_seqz_purity$ABSOLUTE_purity, method="s")

# Plot Sequenza vs. ABSOLUTE for shared TCGA samples.
ggplot(glass_taylor_seqz_purity, aes(ABSOLUTE_purity, cellularity)) + geom_point() + theme_bw() + xlab("ABSOLUTE purity") +
  ylab("Sequenza purity") +  ggtitle("TCGA glioma samples (n=40 )") + geom_abline(slope = 1) + ylim(0,1) + xlim(0, 1)

# Purity for TITAN vs. ABSOLUTE.
glass_taylor_titan_purity = titan_purity %>%
  inner_join(glioma_dat, by="sample_barcode")

# A weaker correlation than the Sequenza samples.
cor.test(glass_taylor_titan_purity$purity, glass_taylor_titan_purity$ABSOLUTE_purity, method="s")

ggplot(glass_taylor_titan_purity, aes(ABSOLUTE_purity, purity)) + geom_point() + theme_bw() + xlab("ABSOLUTE purity") +
  ylab("TITAN purity") + ggtitle("TCGA glioma samples (n=44)") + geom_abline(slope = 1) + ylim(0,1) + xlim(0, 1)

# Plot how the the purity estimates change between primary and recurrent tumors.
titan_compare = titan_purity %>% 
  mutate(sample_type = recode(sample_type,  "tumor_barcode_a" = "Primary", "tumor_barcode_b" = "Recurrence")) %>% 
  select(tumor_pair_barcode, case_barcode, sample_type, purity) %>% 
  spread(sample_type, purity) %>% 
  left_join(subtypes, by="case_barcode")
# How many samples have all the data? n = 250.
sum(complete.cases(titan_compare))

# Test for subtypes across timepoints.
kruskal.test(titan_compare$Primary ~ as.factor(titan_compare$idh_codel_subtype))
kruskal.test(titan_compare$Recurrence ~ as.factor(titan_compare$idh_codel_subtype))

# Perform statistical tests for purity differences between two samples.
titan_compare_codel = titan_compare %>% filter(idh_codel_subtype =="IDHmut-codel")
titan_compare_noncodel = titan_compare %>% filter(idh_codel_subtype =="IDHmut-noncodel")
titan_compare_IDHwt  = titan_compare %>% filter(idh_codel_subtype =="IDHwt")
wilcox.test(titan_compare_codel$Primary, titan_compare_codel$Recurrence, paired = TRUE)
wilcox.test(titan_compare_noncodel$Primary, titan_compare_noncodel$Recurrence, paired = TRUE)
wilcox.test(titan_compare_IDHwt$Primary, titan_compare_IDHwt$Recurrence, paired = TRUE)

# Prepare the data for a side-by-side boxplot.
titan_plot_purity = titan_compare %>% 
  gather(sample, purity, c(Primary, Recurrence))
pdf(file = "/Users/johnsk/Documents/TITAN_purity_silverset_n250.pdf", height = 5, width = 7, bg = "transparent", useDingbats = FALSE)
ggplot(titan_plot_purity, aes(x=idh_codel_subtype, y=purity, fill=sample)) + 
  geom_boxplot()  + ylab("TITAN Purity") + xlab("") + theme_bw() + scale_fill_manual(values = c("Primary" = "#a6611a", "Recurrence" = "#018571")) + labs(fill="Sample Type") +
  ylim(0, 1.1) + annotate("text", label = c("P = 0.68", "P = 0.33", "P = 0.55"),  x = c(1, 2, 3), y = c(1.1, 1.1, 1.1), size = 4) + ggtitle("GLASS silver set (n = 250 subjects)")
dev.off()



#########################
# Comparisons of purity in gold_set versus non-gold-set samples
#########################
titan_purity = titan_params %>% 
  inner_join(pairs, by= "pair_barcode") %>% 
  inner_join(gold_set_aliquots, by=c("tumor_barcode"="aliquot_barcode")) %>% 
  inner_join(mutation_freq, by=c("tumor_barcode"="aliquot_barcode")) %>% 
  mutate(sample_barcode = substr(pair_barcode, 1, 15),
         seq_type = substr(pair_barcode, 27, 29)) %>% 
  left_join(subtypes, by="case_barcode") 

# Remove any gold set samples and only keep the review|allow samples
titan_anti_purity = titan_params %>% 
  inner_join(pairs, by= "pair_barcode") %>% 
  anti_join(gold_set_aliquots, by=c("tumor_barcode"="aliquot_barcode")) %>% 
  left_join(blocklist, by=c("tumor_barcode"="aliquot_barcode")) %>% 
  filter(cnv_exclusion!="allow ")

# Is there a difference?
wilcox.test(titan_purity$purity, titan_anti_purity$purity)

# TITAN seems to provide higher estimates for these difficult samples.
ggplot(titan_purity, aes(y=purity)) + geom_boxplot()
ggplot(titan_anti_purity, aes(y=purity)) + geom_boxplot()

## Try Sequenza ##
# Restrict to gold_set purity estimates.
  seqz_purity = seqz_params %>% 
    inner_join(pairs, by= "pair_barcode") %>% 
    inner_join(gold_set_aliquots, by=c("tumor_barcode"="aliquot_barcode")) %>% 
    inner_join(mutation_freq, by=c("tumor_barcode"="aliquot_barcode")) %>% 
    mutate(sample_barcode = substr(pair_barcode, 1, 15),
           seq_type = substr(pair_barcode, 27, 29)) %>% 
    left_join(subtypes, by="case_barcode") 
  
# Same as with TITAN.
  seqz_anti_purity = seqz_params %>% 
    inner_join(pairs, by= "pair_barcode") %>% 
    anti_join(gold_set_aliquots, by=c("tumor_barcode"="aliquot_barcode")) %>% 
    left_join(blocklist, by=c("tumor_barcode"="aliquot_barcode")) %>% 
    filter(cnv_exclusion!="allow ")
  
# Strong statistical association in the expected direction.
wilcox.test(seqz_purity$cellularity, seqz_anti_purity$cellularity)
ggplot(seqz_purity, aes(y=cellularity)) + geom_boxplot()
ggplot(seqz_anti_purity, aes(y=cellularity)) + geom_boxplot()
  
  