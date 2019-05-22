##### Enrichment of genes in Recurrence analysis ######
# Which genes are enriched in the recurrence?
# Specific focus on Recurrence only (R) vs Primary only (P)

###########################################
################ Settings #################
setwd("/Users/c-kocake/Box Sync/Projects/GLASS/Pathway/")
library(tidyverse)
library(DBI)

# Create connection to database and load necessary tables
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

##################################################
################ Loading in data #################
# Load in variables from the database.
gold_set = dbReadTable(con,  Id(schema="analysis", table="gold_set"))
tumor_mut_comp_anno = dbReadTable(con, Id(schema="analysis", table = "tumor_mut_comparison_anno"))
hypermut_anno = tumor_mut_comp_anno %>% select(case_barcode, hypermutator_status)

# read in heatmap query (Fig 3 mutational landscape)
heatmap_snv_query <- "/*
                                                                                       */
                                                                                       WITH
                                                                                       selected_tumor_pairs AS
                                                                                       (
                                                                                       SELECT * FROM analysis.gold_set
                                                                                       ),
                                                                                       selected_aliquots AS
                                                                                       (
                                                                                       SELECT tumor_barcode_a AS aliquot_barcode, case_barcode FROM selected_tumor_pairs
                                                                                       UNION
                                                                                       SELECT tumor_barcode_b AS aliquot_barcode, case_barcode FROM selected_tumor_pairs
                                                                                       ),
                                                                                       selected_genes AS
                                                                                       (
                                                                                       SELECT DISTINCT sn.gene_symbol, ensembl_gene_id, variant_id, chrom, pos, alt, sn.variant_classification, variant_classification_priority, protein_change
                                                                                       FROM variants.anno sn
                                                                                       INNER JOIN ref.driver_genes ds ON ds.gene_symbol = sn.gene_symbol OR sn.gene_symbol IN ('NOTCH1','TCF12','BCOR','STAG2','PLCG1','ZBTB20','NUP210L','GABRA6','ABCD1','SLC6A3','SETD2','PTPN11','TRPA1','NRAS','ST3GAL6','KRAS','ARID2','RPL5','DNMT3A','COL6A3','TYRP1','LZTR1','EEF1A1','ZNF292','KRT15','FAM126B','TLR6','MAX','BRAF','NIPBL','TP63','RNF168','AOX1','DDX5','QKI','MUC17','F5','SMARCA4','EIF1AX','IL1RL1','ZDHHC4','ACAN','ARID1A','C10orf76','C10orf48','CLCN7','CLIP1','CREBZF','EDAR','EPHA3','FHOD1','GLT8D2','HMCN1','IL4R','KEL','KRT13','MYT1','PIK3C2G','PPM1J','RBBP6','SEMG1','TGFA','TPTE2') 
                                                                                       INNER JOIN ref.ensembl_gene_mapping gm ON gm.gene_symbol = sn.gene_symbol
                                                                                       LEFT JOIN variants.variant_classifications vc ON sn.variant_classification = vc.variant_classification
                                                                                       WHERE
                                                                                       has_mut IS TRUE AND
                                                                                       ((sn.gene_symbol NOT IN ('TERT','IDH1') AND variant_classification_priority IS NOT NULL) OR
                                                                                       (sn.gene_symbol = 'TERT' AND sn.variant_classification = 'FIVE_PRIME_FLANK' AND lower(sn.pos) IN (1295228,1295250)) OR
                                                                                       (sn.gene_symbol = 'IDH1' AND sn.protein_change IN ('p.R132C','p.R132G','p.R132H','p.R132S')))
                                                                                       ),
                                                                                       selected_genes_uq AS
                                                                                       (
                                                                                       SELECT DISTINCT sg.gene_symbol, ensembl_gene_id
                                                                                       FROM selected_genes sg
                                                                                       ),
                                                                                       selected_genes_samples AS
                                                                                       (
                                                                                       SELECT aliquot_barcode, case_barcode, gene_symbol, ensembl_gene_id, chrom, lower(pos) AS start_pos, upper(pos)-1 AS end_pos, alt
                                                                                       FROM selected_aliquots, selected_genes
                                                                                       ),
                                                                                       selected_genes_samples_uq AS
                                                                                       (
                                                                                       SELECT DISTINCT aliquot_barcode, case_barcode, gene_symbol, ensembl_gene_id
                                                                                       FROM selected_genes_samples
                                                                                       ),
                                                                                       selected_variants_samples AS
                                                                                       (
                                                                                       SELECT variant_id, tumor_pair_barcode, variant_classification_priority, protein_change
                                                                                       FROM selected_genes
                                                                                       CROSS JOIN selected_tumor_pairs
                                                                                       ),
                                                                                       hotspot_coverage AS -- For IDH1 and TERT we don't want genic coverage but coverage at the hotspot sites
                                                                                       (
                                                                                         SELECT aliquot_barcode, case_barcode, gene_symbol, sum(ad_alt + ad_ref)::decimal / COUNT(*) AS gene_cov
                                                                                         FROM variants.geno pg
                                                                                         INNER JOIN variants.passanno pa ON pa.variant_id = pg.variant_id
                                                                                         LEFT JOIN variants.variant_classifications vc ON vc.variant_classification = pa.variant_classification
                                                                                         WHERE
                                                                                         (pa.gene_symbol = 'TERT' AND pa.variant_classification = 'FIVE_PRIME_FLANK' AND lower(pa.pos) IN (1295228,1295250)) OR
                                                                                         (pa.gene_symbol = 'IDH1' AND pa.protein_change IN ('p.R132C','p.R132G','p.R132H','p.R132S'))
                                                                                         GROUP BY 1,2,3
                                                                                         ORDER BY 1,2,3
                                                                                       ),
                                                                                       gene_sample_coverage AS -- For all other genes get genic coverage
                                                                                       (
                                                                                         SELECT sgu.aliquot_barcode, case_barcode, sgu.gene_symbol, gene_coverage::double precision / gene_size AS gene_cov
                                                                                         FROM selected_genes_samples_uq sgu
                                                                                         INNER JOIN analysis.gencode_coverage gc ON gc.ensembl_gene_id = sgu.ensembl_gene_id AND gc.aliquot_barcode = sgu.aliquot_barcode
                                                                                         INNER JOIN ref.ensembl_genes eg ON eg.ensembl_gene_id = sgu.ensembl_gene_id
                                                                                         WHERE sgu.gene_symbol NOT IN ('IDH1','TERT')
                                                                                         
                                                                                         UNION
                                                                                         
                                                                                         SELECT * FROM hotspot_coverage
                                                                                       ),
                                                                                       gene_pair_coverage AS
                                                                                       (
                                                                                         SELECT
                                                                                         stp.tumor_pair_barcode,
                                                                                         stp.case_barcode,
                                                                                         sg.gene_symbol,
                                                                                         c1.gene_cov AS cov_a,
                                                                                         c2.gene_cov AS cov_b
                                                                                         FROM selected_tumor_pairs stp
                                                                                         CROSS JOIN selected_genes_uq sg
                                                                                         LEFT JOIN gene_sample_coverage c1 ON c1.aliquot_barcode = stp.tumor_barcode_a AND c1.gene_symbol = sg.gene_symbol
                                                                                         LEFT JOIN gene_sample_coverage c2 ON c2.aliquot_barcode = stp.tumor_barcode_b AND c2.gene_symbol = sg.gene_symbol
                                                                                       ),
                                                                                       variants_by_case_and_gene AS
                                                                                       (
                                                                                         SELECT
                                                                                         gtc.gene_symbol,
                                                                                         gtc.case_barcode,
                                                                                         gtc.tumor_pair_barcode,
                                                                                         gtc.chrom,
                                                                                         gtc.pos,
                                                                                         gtc.alt,
                                                                                         vc.variant_classification_vep AS variant_classification,
                                                                                         protein_change,
                                                                                         alt_count_a > 0 AS selected_call_a, --(alt_count_a::decimal / (alt_count_a+ref_count_a) > 0.05)
                                                                                         alt_count_b > 0 AS selected_call_b, --(alt_count_b::decimal / (alt_count_b+ref_count_b) > 0.05)
                                                                                         row_number() OVER (PARTITION BY gtc.gene_symbol, gtc.case_barcode ORDER BY mutect2_call_a::integer + mutect2_call_b::integer = 2, vc.variant_classification_priority, mutect2_call_a::integer + mutect2_call_b::integer DESC, ref_count_a + ref_count_b + alt_count_a + alt_count_b DESC) AS priority
                                                                                         FROM variants.pgeno gtc
                                                                                         INNER JOIN selected_variants_samples svs ON svs.tumor_pair_barcode = gtc.tumor_pair_barcode AND svs.variant_id = gtc.variant_id
                                                                                         INNER JOIN variants.variant_classifications vc ON vc.variant_classification = gtc.variant_classification
                                                                                         WHERE
                                                                                         (mutect2_call_a OR mutect2_call_b) AND (alt_count_a+ref_count_a) >= 5 AND (alt_count_b+ref_count_b) >= 5
                                                                                       ),
                                                                                       squared_variants AS
                                                                                       (
                                                                                         SELECT
                                                                                         gc.gene_symbol,
                                                                                         gc.case_barcode,
                                                                                         vcg.chrom,
                                                                                         vcg.pos,
                                                                                         vcg.alt,
                                                                                         vcg.variant_classification,
                                                                                         vcg.protein_change,
                                                                                         vcg.selected_call_a,
                                                                                         vcg.selected_call_b,
                                                                                         gc.cov_a,
                                                                                         gc.cov_b
                                                                                         FROM (SELECT * FROM variants_by_case_and_gene WHERE priority = 1) vcg
                                                                                         RIGHT JOIN gene_pair_coverage gc ON gc.tumor_pair_barcode = vcg.tumor_pair_barcode AND gc.gene_symbol = vcg.gene_symbol
                                                                                       )
                                                                                       --SELECT * --sa.aliquot_barcode, case_barcode, gm.gene_symbol, gene_coverage::double precision / gene_size AS gene_cov
                                                                                       --FROM selected_genes_samples_uq sgu
                                                                                       --INNER JOIN analysis.gencode_coverage gc ON gc.ensembl_gene_id = sgu.ensembl_gene_id AND gc.aliquot_barcode = sgu.aliquot_barcode-- 
                                                                                         --SELECT * from gene_sample_coverage --gene_pair_coverage gc
                                                                                       --SELECT * FROM variants_by_case_and_gene vcg
                                                                                       --RIGHT JOIN gene_pair_coverage gc ON gc.tumor_pair_barcode = vcg.tumor_pair_barcode AND gc.gene_symbol = vcg.gene_symbol
                                                                                       SELECT
                                                                                       gene_symbol,
                                                                                       var.case_barcode,
                                                                                       idh_codel_subtype,
                                                                                       variant_classification,
                                                                                       protein_change,
                                                                                       (CASE
                                                                                        WHEN cov_a >= 5 AND cov_b >= 5 AND selected_call_a AND selected_call_b 				THEN 'S'
                                                                                        WHEN cov_a >= 5 AND cov_b >= 5 AND selected_call_a AND NOT selected_call_b 			THEN 'P'
                                                                                        WHEN cov_a >= 5 AND cov_b >= 5 AND selected_call_b AND NOT selected_call_a 			THEN 'R'
                                                                                        ELSE NULL END) AS variant_call,
                                                                                       (CASE
                                                                                        WHEN cov_a >= 30 AND cov_b >= 30 THEN 'Coverage >= 30x' 
                                                                                        --WHEN cov_a >= 20 AND cov_b >= 20 THEN 'Coverage >= 20x' 
                                                                                        --WHEN cov_a >= 10 AND cov_b >= 5 THEN 'Coverage >= 10x' 
                                                                                        WHEN cov_a >= 15 AND cov_b >= 15 THEN 'Coverage >= 15x' 
                                                                                        WHEN cov_a >= 5 AND cov_b >= 5 THEN 'Coverage >= 5x' 
                                                                                        ELSE 'Coverage < 5x' END) AS covered
                                                                                       FROM squared_variants var
                                                                                       LEFT JOIN clinical.subtypes st ON st.case_barcode = var.case_barcode
                                                                                       --LEFT JOIN selected_tumor_pairs stp ON stp.case_barcode = var.case_barcode AND priority = 1
                                                                                       --LEFT JOIN clinical.surgeries su ON su.sample_barcode = substring(tumor_barcode_a from 1 for 15)"
heatmap_snv_df <- DBI::dbGetQuery(con, heatmap_snv_query)
heatmap_snv_df_unique <- heatmap_snv_df %>% unique()

heatmap_cnv_query <- read_file('/Users/c-kocake/Box Sync/Scripts/repos/GLASS/sql/heatmap/heatmap_cnv.sql')
heatmap_cnv_df <- DBI::dbGetQuery(con, heatmap_cnv_query)

#define genes of interest 
genes_of_interest_snv <- c("EGFR", "ATRX", "PTEN", "TP53", "TERT","PIK3CA", "PIK3R1", "CIC", "RB1", "NF1", "FUBP1", "IDH1")

genes_of_interest_cnv <- c("EGFR", "ATRX", "PDGFRA", "CCND2", "CDK4", "CDK6", "CDKN2A", "PTEN", "RB1", "MDM2", "MDM4", "MET", "MYCN")

###########################
#snv
heatmap_snv <- heatmap_snv_df_unique %>% 
  filter(case_barcode %in% gold_set$case_barcode, gene_symbol %in% genes_of_interest_snv) %>% 
  mutate(variant_call = case_when(covered == "Coverage < 5x" & is.na(variant_call) ~ "low_coverage", 
                                  covered != "Coverage < 5x" & is.na(variant_call) ~ "none",
                                  TRUE ~ variant_call)) %>% 
  filter(variant_call != "low_coverage") %>% 
  add_count(idh_codel_subtype, gene_symbol, name = "n_genes") %>% 
  add_count(idh_codel_subtype, gene_symbol, variant_call, name = "n_PRS") %>% 
  mutate(prop = n_PRS/n_genes) %>% 
  select(gene_symbol, idh_codel_subtype, variant_call, prop) %>% distinct() %>% 
  spread(key = variant_call, value = prop) %>% 
  mutate_all(~replace(., is.na(.), 0)) %>%
  gather(key = variant_call, value = prop, -idh_codel_subtype, -gene_symbol) %>% 
  mutate(variant_call = factor(variant_call,
                               levels = c("none", "R", "P", "S")))

snv_plot <- ggplot(heatmap_snv, aes(x=gene_symbol, y = prop, fill = variant_call)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "Proportion", x = "Gene Symbol", title = "SNV Gold-Set (all)") +
  scale_fill_manual(values = c("none" = "grey", "R" = "#2FB3CA", "P" ="#CA2F66", "S"="#CA932F"), 
                    name = " ",
                    labels = c("none" = "None", "R" = "Recurrence", "P" = "Initial", "S" = "Shared")) +
  coord_flip() + theme_classic2() + geom_vline(xintercept = 1.5:11.5) +
  facet_grid(~idh_codel_subtype)


#snv non_hyper
heatmap_snv_nonhypermutant <- heatmap_snv_df_unique %>% 
  filter(case_barcode %in% gold_set$case_barcode, gene_symbol %in% genes_of_interest_snv) %>% 
  mutate(variant_call = case_when(covered == "Coverage < 5x" & is.na(variant_call) ~ "low_coverage", 
                                  covered != "Coverage < 5x" & is.na(variant_call) ~ "none",
                                  TRUE ~ variant_call)) %>% 
  filter(variant_call != "low_coverage") %>% 
  left_join(hypermut_anno, by = "case_barcode") %>% filter(hypermutator_status == 0) %>% 
  add_count(idh_codel_subtype, gene_symbol, name = "n_genes") %>% 
  add_count(idh_codel_subtype, gene_symbol, variant_call, name = "n_PRS") %>% 
  mutate(prop = n_PRS/n_genes) %>% 
  select(gene_symbol, idh_codel_subtype, variant_call, prop) %>% distinct() %>% 
  spread(key = variant_call, value = prop) %>% 
  mutate_all(~replace(., is.na(.), 0)) %>%
  gather(key = variant_call, value = prop, -idh_codel_subtype, -gene_symbol) %>% 
  mutate(variant_call = factor(variant_call,
                               levels = c("none", "R", "P", "S")))

snv_plot_nonhypermutant <- ggplot(heatmap_snv_nonhypermutant, aes(x=gene_symbol, y = prop, fill = variant_call)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "Proportion", x = "Gene Symbol", title = "SNV Gold-Set (Non-Hypermutant)") +
  scale_fill_manual(values = c("none" = "grey", "R" = "#2FB3CA", "P" ="#CA2F66", "S"="#CA932F"), 
                    name = " ",
                    labels = c("none" = "None", "R" = "Recurrence", "P" = "Initial", "S" = "Shared")) +
  coord_flip() + theme_classic2() + geom_vline(xintercept = 1.5:11.5) +
  facet_grid(~idh_codel_subtype) #+
#geom_text(data = s_nonhyper, aes(label = sprintf("n=%s", fishers_exact)), x= gene_symbol,  y = 0.90)

pdf(file = "/Users/c-kocake/Box Sync/Projects/GLASS/Pathway/SNV.pdf", height = 8, width = 17, bg = "transparent", useDingbats = FALSE)
ggarrange(snv_plot, snv_plot_nonhypermutant, ncol=2, nrow=1, common.legend = T, legend = "right")
dev.off()




###################
#significance test snv nonhyper
sig_snv_nonhyper <- heatmap_snv_df_unique %>% 
  filter(case_barcode %in% gold_set$case_barcode) %>% 
  mutate(variant_call = case_when(covered == "Coverage < 5x" & is.na(variant_call) ~ "low_coverage", 
                                  covered != "Coverage < 5x" & is.na(variant_call) ~ "none",
                                  TRUE ~ variant_call)) %>% 
  filter(variant_call != "low_coverage") %>% 
  left_join(hypermut_anno, by = "case_barcode") %>% filter(hypermutator_status == 0)

s_nonhyper <- as.data.frame(table(sig_snv_nonhyper$gene_symbol, sig_snv_nonhyper$variant_call, sig_snv_nonhyper$idh_codel_subtype)) %>% spread(key = Var2, value = Freq) %>% mutate(patients = none + P + R + S, 
                                                                                                                                                                                    P_alt = P, 
                                                                                                                                                                                    P_no_alt = patients - P - S , 
                                                                                                                                                                                    R_alt = R,
                                                                                                                                                                                    R_no_alt = patients - R - S)
#create function for fisher.test (p-value)
function_fisher_p_value <- function(y){
  # include here as.numeric to be sure that your values are numeric:
  table <-  matrix(as.numeric(c(y[8], y[9], y[10], y[11])), ncol = 2, byrow = TRUE)
  if(any(is.na(table))) p <- "error" else p <- fisher.test(table, alternative="two.sided")$p.value
  p
}  
#create new column within existing dataframe and apply tests including post-hoc tests
s_nonhyper$fishers_exact <- apply(s_nonhyper, 1, function_fisher_p_value)
s_nonhyper$adjusted_FDR <- p.adjust(s_nonhyper$fishers_exact, method = "fdr")
s_nonhyper$adjusted_BH <- p.adjust(s_nonhyper$fishers_exact, method = "BH")

#significance test snv hyper
sig_snv_hyper <- heatmap_snv_df_unique %>% 
  filter(case_barcode %in% gold_set$case_barcode) %>% 
  mutate(variant_call = case_when(covered == "Coverage < 5x" & is.na(variant_call) ~ "low_coverage", 
                                  covered != "Coverage < 5x" & is.na(variant_call) ~ "none",
                                  TRUE ~ variant_call)) %>% 
  filter(variant_call != "low_coverage") %>% 
  left_join(hypermut_anno, by = "case_barcode") %>% filter(hypermutator_status == 1)

s_hyper <- as.data.frame(table(sig_snv_hyper$gene_symbol, sig_snv_hyper$variant_call, sig_snv_hyper$idh_codel_subtype)) %>% spread(key = Var2, value = Freq) %>% mutate(patients = none + P + R + S, 
                                                                                                                                                                        P_alt = P, 
                                                                                                                                                                        P_no_alt = patients - P - S , 
                                                                                                                                                                        R_alt = R,
                                                                                                                                                                        R_no_alt = patients - R - S) 




#create new column within existing dataframe and apply tests including post-hoc tests
s_hyper$fishers_exact <- apply(s_hyper, 1, function_fisher_p_value)
s_hyper$adjusted_FDR <- p.adjust(s_hyper$fishers_exact, method = "fdr")
s_hyper$adjusted_BH <- p.adjust(s_hyper$fishers_exact, method = "BH")


#### plot significances for hyper/nonhyper fishers exact and fdr
sig_plot_nonhyper <- ggplot(s_nonhyper, aes(x = -log10(fishers_exact), y = Var1)) +
  geom_point() +
  geom_vline(xintercept = -log10(0.05)) +
  facet_grid(~Var3) +labs (y = "Gene", x = "-log10(p-value)", subtitle = "SNV Non-Hypermutants (Fisher's exact)")

sig_plot_hyper_fisher <- ggplot(s_hyper, aes(x = -log10(fishers_exact), y = Var1)) +
  geom_point() +
  geom_vline(xintercept = -log10(0.05)) +
  facet_grid(~Var3) +labs (y = "Gene", x = "-log10(p-value)", subtitle = "SNV Hypermutants (Fisher's exact)")

sig_plot_hyper_fdr <- ggplot(s_hyper, aes(x = -log10(adjusted_FDR), y = Var1)) +
  geom_point() +
  geom_vline(xintercept = -log10(0.05)) +
  facet_grid(~Var3) +labs (y = "Gene", x = "-log10(adjusted p-value)", subtitle = "SNV Hypermutants (Adjusted FDR")

pdf(file = "/Users/c-kocake/Box Sync/Projects/GLASS/Pathway/SNV_significance_plots.pdf", height =12, width = 17, bg = "transparent", useDingbats = FALSE)
ggarrange(sig_plot_nonhyper, sig_plot_hyper_fisher, sig_plot_hyper_fdr, ncol=3, nrow=1)
dev.off()

######################################################
#cnv
heatmap_cnv <- heatmap_cnv_df %>% 
  filter(gene_symbol %in% genes_of_interest_cnv) %>% 
  add_count(idh_codel_subtype, gene_symbol, name = "n_genes") %>% 
  add_count(idh_codel_subtype, gene_symbol, cnv_change, name = "n_PRS") %>% 
  mutate(prop = n_PRS/n_genes) %>% 
  select(gene_symbol, idh_codel_subtype, cnv_change, prop) %>% distinct() %>%
  spread(key = cnv_change, value = prop) %>% 
  mutate_all(~replace(., is.na(.), 0)) %>% 
  rename("No_alteration" = "<NA>") %>% 
  gather(key=type, value=prop, -idh_codel_subtype, -gene_symbol) %>% 
  mutate(type = factor(type,
                       levels = c("No_alteration", "R", "P", "S")))

cnv_plot <- ggplot(heatmap_cnv,aes(x=gene_symbol, y = prop, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "Proportion", x = "Gene Symbol", title = "CNV") +
  scale_fill_manual(values=c("No_alteration" = "grey", "R" = "#2FB3CA", "P" ="#CA2F66", "S"="#CA932F")) +
  coord_flip() + theme_classic2() + geom_vline(xintercept = 1.5:12.5) +
  facet_grid(~idh_codel_subtype)

#significance test cnv
sig_cnv <- heatmap_cnv_df %>% filter(gene_symbol %in% genes_of_interest_cnv) %>% mutate(variant_call = ifelse(is.na(cnv_change), "No_alteration", cnv_change))
c <- as.data.frame(table(sig_cnv$gene_symbol, sig_cnv$variant_call, sig_cnv$idh_codel_subtype)) %>% spread(key = Var2, value = Freq) %>% mutate(patients = No_alteration + P + R + S, 
                                                                                                                                                P_alt = P, 
                                                                                                                                                P_no_alt = patients - P - S, 
                                                                                                                                                R_alt = R,
                                                                                                                                                R_no_alt = patients - R - S)
#create new column within existing dataframe and apply tests including post-hoc tests
c$fishers_exact <- apply(c, 1, function_fisher_p_value)
c$adjusted_FDR <- p.adjust(c$fishers_exact, method = "fdr")
c$adjusted_BH <- p.adjust(c$fishers_exact, method = "BH")

# ----> only significant is CDKN2A(noncodels), EGFR(wt) and this disappears after multiple correction

pdf(file = "/Users/c-kocake/Box Sync/Projects/GLASS/Pathway/gene_comparison.pdf", height = 11, width = 18, bg = "transparent", useDingbats = FALSE)
ggarrange(snv_plot_nonhypermutant, cnv_plot, ncol = 2, nrow = 1, common.legend = T, align = "hv", legend = "right")
dev.off()
### save environment
save.image("GLASS_gene_comparison_clean.RData")

###########################################################################
################################### END ###################################
###########################################################################
