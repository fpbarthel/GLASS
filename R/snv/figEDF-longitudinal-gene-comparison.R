#########################################
# Longitudinal gene comparison analysis #
#########################################

###########################################
################ Settings #################
setwd("/Users/c-kocake/Box Sync/Scripts/repos/GLASS/")
library(tidyverse)
library(DBI)
library(ggpubr)

# Create connection to database and load necessary tables
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

##################################################
################ Loading in data #################
# Load in variables from the database.
gold_set = dbReadTable(con,  Id(schema="analysis", table="gold_set"))
tumor_mut_comp_anno = dbReadTable(con, Id(schema="analysis", table = "tumor_mut_comparison_anno"))
hypermut_anno = tumor_mut_comp_anno %>% select(case_barcode, hypermutator_status)

# Load in mutation data from the database for all genes
mutation_query <- read_file("sql/snv/longitudinal_gene_comparison_snv_all_genes.sql")
mutation <- DBI::dbGetQuery(con, mutation_query)

# load in SNV and CNV data from the database for driver genes 
heatmap_snv_query <- read_file("sql/snv/longitudinal_gene_comparison_snv_smg.sql")
heatmap_snv_df <- DBI::dbGetQuery(con, heatmap_snv_query)
heatmap_snv_df_unique <- heatmap_snv_df %>% unique()

heatmap_cnv_query <- read_file("sql/heatmap/heatmap_cnv.sql")
heatmap_cnv_df <- DBI::dbGetQuery(con, heatmap_cnv_query)

#define genes of interest 
genes_of_interest_snv <- c("EGFR", "ATRX", "PTEN", "TP53", "TERT","PIK3CA", "PIK3R1", "CIC", "RB1", "NF1", "FUBP1", "IDH1")

genes_of_interest_cnv <- c("EGFR", "ATRX", "PDGFRA", "CCND2", "CDK4", "CDK6", "CDKN2A", "PTEN", "RB1", "MDM2", "MDM4", "MET", "MYCN")

mutated_genes <- mutation %>% filter(!is.na(variant_call)) %>% select(gene_symbol) %>% unique()
################### Significance basic dataframe ###################
sig <- mutation %>% filter(gene_symbol %in% mutated_genes$gene_symbol) %>% 
  mutate(variant_call = case_when(covered == "Coverage < 5x" & is.na(variant_call) ~ "low_coverage", 
                                  covered != "Coverage < 5x" & is.na(variant_call) ~ "none",
                                  TRUE ~ variant_call)) %>% 
  filter(variant_call != "low_coverage") %>% left_join(hypermut_anno, by = "case_barcode")

sig_nonhyper <- sig %>% filter(hypermutator_status == 0)
sig_hyper <- sig %>% filter(hypermutator_status == 1)
sig_nonhyper_smg <- sig_nonhyper %>% filter(gene_symbol %in% heatmap_snv_df_unique$gene_symbol)
sig_hyper_smg <- sig_hyper %>% filter(gene_symbol %in% heatmap_snv_df_unique$gene_symbol)
#########################################################

################### fisher function ###################
#create function for fisher.test (p-value)
function_fisher_p_value <- function(y){
  table <-  matrix(as.numeric(c(y[7], y[8], y[9], y[10])), ncol = 2, byrow = TRUE)
  if(any(is.na(table))) p <- "error" else p <- fisher.test(table,alternative = "less")$p.value
  p}  
function_fisher_p_value_subtype <- function(y){
  table <-  matrix(as.numeric(c(y[8], y[9], y[10], y[11])), ncol = 2, byrow = TRUE)
  if(any(is.na(table))) p <- "error" else p <- fisher.test(table,alternative = "two.sided" )$p.value
  p} 

#create new functions for applying fishers exact in combined table
fisher_hyper <- function(y){
  table <-  matrix(as.numeric(c(y[3], y[4], y[5], y[6])), ncol = 2, byrow = TRUE)
  if(any(is.na(table))) p <- "error" else p <- fisher.test(table,alternative = "less")$p.value
  p}  
fisher_nonhyper <- function(y){
  table <-  matrix(as.numeric(c(y[9], y[10], y[11], y[12])), ncol = 2, byrow = TRUE)
  if(any(is.na(table))) p <- "error" else p <- fisher.test(table,alternative = "less")$p.value
  p}  
#########################################################

################### Significance testing with all genes ###################
#significance test snv nonhyper all genes
s_nonhyper <- as.data.frame(table(sig_nonhyper$gene_symbol, sig_nonhyper$variant_call)) %>% spread(key = Var2, value = Freq) %>% mutate(patients = none + P + R + S, 
                                                                                                                                                P_alt = P, 
                                                                                                                                                P_no_alt = patients - P - S , 
                                                                                                                                                R_alt = R,
                                                                                                                                                R_no_alt = patients - R - S) %>% 
  mutate(fishers_exact = apply(., 1, function_fisher_p_value), adjusted_FDR = p.adjust(fishers_exact, method = "fdr"))

#significance test snv hyper all genes
s_hyper <- as.data.frame(table(sig_hyper$gene_symbol, sig_hyper$variant_call)) %>% spread(key = Var2, value = Freq) %>% mutate(patients = none + P + R + S, 
                                                                                                                                       P_alt = P, 
                                                                                                                                       P_no_alt = patients - P - S , 
                                                                                                                                       R_alt = R,
                                                                                                                                       R_no_alt = patients - R - S) %>% 
  mutate(fishers_exact = apply(., 1, function_fisher_p_value), adjusted_FDR = p.adjust(fishers_exact, method = "fdr"))
#########################################################

################### Significance testing with smgs ###################
#significance test snv nonhyper smgs
s_nonhyper_smg <- as.data.frame(table(sig_nonhyper_smg$gene_symbol, sig_nonhyper_smg$variant_call)) %>% spread(key = Var2, value = Freq) %>% mutate(patients = none + P + R + S, 
                                                                                                                                           P_alt = P, 
                                                                                                                                           P_no_alt = patients - P - S , 
                                                                                                                                           R_alt = R,
                                                                                                                                           R_no_alt = patients - R - S) %>% 
  mutate(fishers_exact = apply(., 1, function_fisher_p_value), adjusted_FDR = p.adjust(fishers_exact, method = "fdr"))

#significance test snv hyper smgs
s_hyper_smg <- as.data.frame(table(sig_hyper_smg$gene_symbol, sig_hyper_smg$variant_call)) %>% spread(key = Var2, value = Freq) %>% mutate(patients = none + P + R + S, 
                                                                                                                                  P_alt = P, 
                                                                                                                                  P_no_alt = patients - P - S , 
                                                                                                                                  R_alt = R,
                                                                                                                                  R_no_alt = patients - R - S) %>% 
  mutate(fishers_exact = apply(., 1, function_fisher_p_value), adjusted_FDR = p.adjust(fishers_exact, method = "fdr"))
#########################################################

################### Significance testing with smgs Hypermutant vs Non-Hypermutant ###################
#combine s_hyper and s_nonhyper for simultaneous visualization
s_both_smg <- s_hyper_smg %>% left_join(s_nonhyper_smg, by = c("Var1")) %>% 
  select(gene_symbol = Var1, 
         patients_hyper = patients.x, P_alt_hyper = P_alt.x, P_no_alt_hyper = P_no_alt.x, R_alt_hyper = R_alt.x, R_no_alt_hyper = R_no_alt.x, S_hyper = S.x,
         patients_nonhyper = patients.y, P_alt_nonhyper = P_alt.y, P_no_alt_nonhyper = P_no_alt.y, R_alt_nonhyper = R_alt.y, R_no_alt_nonhyper = R_no_alt.y, S_nonhyper = S.y) %>% 
  mutate(fishers_exact_hyper = apply(., 1, fisher_hyper), fishers_exact_nonhyper = apply(., 1, fisher_nonhyper), n_hyper = patients_hyper - S_hyper - P_alt_hyper, n_nonhyper = patients_nonhyper - S_nonhyper - P_alt_nonhyper) %>% filter(R_alt_hyper > 3 | R_alt_nonhyper > 3) %>% 
  mutate(adjusted_FDR_hyper = p.adjust(fishers_exact_hyper, method = "fdr"), adjusted_FDR_nonhyper = p.adjust(fishers_exact_nonhyper, method = "fdr")) %>% 
  gather(key = hypermutation, value = value, `adjusted_FDR_hyper`, `adjusted_FDR_nonhyper`) %>% group_by(hypermutation) %>% mutate(gene_symbol = reorder(gene_symbol, -value))

#plot significances
sig_plot_both <- ggplot(s_both_smg,aes(x =-log10(value), y=gene_symbol,  color = hypermutation)) + 
  geom_point() +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(y = "Gene", x = "-log10(adjusted p-value)", subtitle = "Gold-Set (Hypermutators n=35, Non-Hypermutators n=187)", title = "Acquired Mutations at Recurrence") +
  geom_text(data = subset(s_both_smg, hypermutation == "adjusted_FDR_hyper"),aes(x = -log10(value) + 0.2,label = sprintf("n = %s/%s", R_alt_hyper, n_hyper), hjust = 0.5, vjust = 0.3, fontface = "italic"), show.legend = F, check_overlap = T) +
  geom_text(data = subset(s_both_smg, hypermutation == "adjusted_FDR_nonhyper"),aes(x = -log10(value) + 0.2,label = sprintf("n = %s/%s", R_alt_nonhyper, n_nonhyper), hjust = 0.5, vjust = 0.3,fontface = "italic"), show.legend = F, check_overlap = T) + theme_bw() +
  scale_colour_manual(values=c("adjusted_FDR_hyper" = "#FF5733", "adjusted_FDR_nonhyper" = "#000000"), 
                      name = "Hypermutation", 
                      labels=c("adjusted_FDR_hyper" = "Yes", "adjusted_FDR_nonhyper" = "No"))

#PDF for sig_plot_both
pdf(file = "/Users/c-kocake/Box Sync/Projects/GLASS/Pathway/significance_plot.pdf", height = 10, width = 10, bg = "transparent", useDingbats = FALSE)
sig_plot_both
dev.off()
#########################################################


################### SNV Plot by subtypes ###################
#snv Non-Hypermutants
snv <- heatmap_snv_df_unique %>% filter(gene_symbol %in% genes_of_interest_snv) %>% 
  mutate(variant_call = case_when(covered == "Coverage < 5x" & is.na(variant_call) ~ "low_coverage", 
                                  covered != "Coverage < 5x" & is.na(variant_call) ~ "none",
                                  TRUE ~ variant_call)) %>% filter(variant_call != "low_coverage") %>% 
  left_join(hypermut_anno, by = "case_barcode") %>% filter(hypermutator_status == 0) %>% 
  add_count(idh_codel_subtype, gene_symbol, variant_call, name = "n_PRS")  %>% select(idh_codel_subtype, gene_symbol, variant_call, n_PRS) %>% unique() %>% 
  spread(key = variant_call, value = n_PRS) %>% mutate_all(~replace(., is.na(.), 0)) %>% 
  mutate( patients = none + R + S + P, P_alt = P, P_no_alt = patients - S - P, R_alt = R, R_no_alt = patients - S - R) %>% 
  mutate(fishers_exact = apply(., 1, function_fisher_p_value_subtype), adjusted_FDR = p.adjust(fishers_exact, method = "fdr")) %>% 
  gather(key = variant_call, value = n_PRS, `none`, `P`, `R`, `S`) %>% 
  mutate(variant_call = factor(variant_call, levels = c("none", "R", "P", "S")), prop = n_PRS/patients)

snv_plot <- ggplot(snv, aes(x=gene_symbol, y = prop, fill = variant_call)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "Proportion", x = "Gene Symbol", title = "SNV", subtitle = "Gold-Set (Non-Hypermutant)") +
  scale_fill_manual(values = c("none" = "grey", "R" = "#2FB3CA", "P" ="#CA2F66", "S"="#CA932F"), 
                    name = " ",
                    labels = c("none" = "None", "R" = "Recurrence", "P" = "Initial", "S" = "Shared")) +
  coord_flip() + theme_classic2() + geom_vline(xintercept = 1.5:11.5) + facet_grid(~idh_codel_subtype) + 
  geom_text(data = subset(snv, variant_call == "P"), aes(y = 0.35,label = sprintf("P = %s", signif(fishers_exact, 2)), hjust = 0.0, vjust = 0.3, fontface = "italic"), show.legend = F, check_overlap = T) +
  geom_errorbar(aes(ymax=0.325, ymin=0.325), width = 0.3, size = 0.1)

segment_snv_1 <- c(1.15:12.15)
for (i in segment_snv_1) {
  snv_plot <- snv_plot + geom_segment(x=i, y=0.3, xend=i, yend=0.325, size = 0.1)
}
segment_snv_2 <- c(0.85:11.85)
for (i in segment_snv_2) {
  snv_plot <- snv_plot + geom_segment(x=i, y=0.3, xend=i, yend=0.325, size = 0.1)
}
snv_plot
#########################################################

################### CNV Plot by subtypes ###################
#cnv
cnv <- heatmap_cnv_df %>% filter(gene_symbol %in% genes_of_interest_cnv) %>% mutate(cnv_change = ifelse(is.na(cnv_change), "none", cnv_change)) %>% 
  add_count(idh_codel_subtype, gene_symbol, cnv_change, name = "n_PRS") %>% select(idh_codel_subtype, gene_symbol, cnv_change, n_PRS) %>% unique() %>% 
  spread(key = cnv_change, value = n_PRS) %>% mutate_all(~replace(., is.na(.), 0)) %>% 
  mutate(patients = none + R + S + P, P_alt = P, P_no_alt = patients - S - P, R_alt = R, R_no_alt = patients - S - R) %>% 
  mutate(fishers_exact = apply(., 1, function_fisher_p_value_subtype), adjusted_FDR = p.adjust(fishers_exact, method = "fdr")) %>% 
  gather(key = cnv_change, value = n_PRS, `none`, `P`, `R`, `S`) %>% 
  mutate(cnv_change = factor(cnv_change, levels = c("none", "R", "P", "S")), prop = n_PRS/patients)

cnv_plot <- ggplot(cnv, aes(x=gene_symbol, y = prop, fill = cnv_change)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "Proportion", x = "Gene Symbol", title = "CNV", subtitle = "Gold-Set") +
  scale_fill_manual(values = c("none" = "grey", "R" = "#2FB3CA", "P" ="#CA2F66", "S"="#CA932F"), 
                    name = " ",
                    labels = c("none" = "None", "R" = "Recurrence", "P" = "Initial", "S" = "Shared")) +
  coord_flip() + theme_classic2() + geom_vline(xintercept = 1.5:11.5) + facet_grid(~idh_codel_subtype) + 
  geom_text(data = subset(cnv, cnv_change == "P"), aes(y = 0.35,label = sprintf("P = %s", signif(fishers_exact, 2)), hjust = 0.0, vjust = 0.3, fontface = "italic"), show.legend = F, check_overlap = T) +
  geom_errorbar(aes(ymax=0.325, ymin=0.325), width = 0.3, size = 0.1)

segment_cnv_1 <- c(1.15:12.15)
for (i in segment_cnv_1) {
  cnv_plot <- cnv_plot + geom_segment(x=i, y=0.3, xend=i, yend=0.325, size = 0.1)
}
segment_cnv_2 <- c(0.85:11.85)
for (i in segment_cnv_2) {
  cnv_plot <- cnv_plot + geom_segment(x=i, y=0.3, xend=i, yend=0.325, size = 0.1)
}
cnv_plot
#########################################################

pdf(file = "/Users/c-kocake/Box Sync/Projects/GLASS/Pathway/gene_comparison_wide.pdf", height = 10, width = 18, bg = "transparent", useDingbats = FALSE)
ggarrange(snv_plot, cnv_plot, ncol = 2, nrow = 1, common.legend = T, align = "hv", legend = "right")
dev.off()

pdf(file = "/Users/c-kocake/Box Sync/Projects/GLASS/Pathway/gene_comparison_long.pdf", height = 15, width = 10, bg = "transparent", useDingbats = FALSE)
ggarrange(snv_plot, cnv_plot, ncol = 1, nrow = 2, common.legend = F, align = "hv", legend = "right")
dev.off()

##### END #######