library(DBI)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

con <- DBI::dbConnect(odbc::odbc(), "GLASSv2")  
res <- dbGetQuery(con, read_file("sql/timing/ccf_shared.sql"))

############################################################################################################################################
## Plot ladder plot comparing initial and recurrence CCF
## For all genes w/ colored driver genes
## Faceted by driver gene
############################################################################################################################################

## Munge data
tmp <- res %>% filter(rank==1) %>% 
  group_by(case_barcode)  %>%
  mutate(w = n()) %>%
  ungroup() %>%
  select(case_barcode, w, idh_codel_subtype, gene_symbol, variant_classification_vep, cellular_prevalence_a, cellular_prevalence_b) %>%
  gather(c(cellular_prevalence_a, cellular_prevalence_b), key = "sample_type", value = "ccf") %>%
  mutate(sample_type = factor(sample_type, levels = c("cellular_prevalence_a", "cellular_prevalence_b"), labels = c("P","R")),
         plot_id = paste(case_barcode,gene_symbol),
         gene_label = factor(case_when(gene_symbol == "TP53" ~ "TP53",
                                gene_symbol == "IDH1" ~ "IDH1",
                                gene_symbol == "ATRX" ~ "ATRX",
                                gene_symbol == "PTEN" ~ "PTEN",
                                gene_symbol == "PIK3CA" ~ "PIK3CA",
                                gene_symbol == "PIK3R1" ~ "PIK3R1",
                                gene_symbol == "NF1" ~ "NF1",
                                gene_symbol == "EGFR" ~ "EGFR",
                                gene_symbol == "CIC" ~ "CIC",
                                gene_symbol == "FUBP1" ~ "FUBP1",
                                gene_symbol == "TERT" ~ "TERT",
                                gene_symbol == "RB1" ~ "RB1",
                                TRUE ~ NA_character_)))

## Perform statistical testing
testWilcoxGroup <- function(df) {
  wtest = wilcox.test(df$ccf ~ df$sample_type, paired = TRUE, conf.int = TRUE)
  data.frame(n = nrow(df)/2,
             median_a = median(df$ccf[df$sample_type=="P"]),
             median_b = median(df$ccf[df$sample_type=="R"]),
             statistic = wtest$statistic,
             estimate = wtest$estimate,
             lcl = wtest$conf.int[1],
             ucl = wtest$conf.int[2],
             wilcox_p = wtest$p.value,
             test_str = sprintf("n=%s\n%s",
                                nrow(df),
                                case_when(wtest$p.value < 0.0001 ~ "P<0.0001",
                                          wtest$p.value > 0.05 ~ sprintf("P=%s", format(round(wtest$p.value, 2),scientific=F)),
                                          TRUE ~ sprintf("P=%s", format(round(wtest$p.value, 4),scientific=F)))),
             stringsAsFactors = FALSE)
}

test_case <- tmp %>% 
  group_by(case_barcode,idh_codel_subtype) %>%
  do(testWilcoxGroup(.)) %>%
  ungroup()

test_subtype <- tmp %>% 
  group_by(idh_codel_subtype) %>%
  do(testWilcoxGroup(.)) %>%
  ungroup()

test_gene <- tmp %>% 
  filter(gene_symbol %in% tmp$gene_label) %>%
  group_by(gene_symbol) %>% 
  do(testWilcoxGroup(.)) %>%
  ungroup() %>%
  mutate(gene_label = gene_symbol)

test_gene_subtype <- tmp %>% 
  filter(gene_symbol %in% c("TP53","IDH1")) %>%
  group_by(gene_symbol,idh_codel_subtype) %>% 
  do(testWilcoxGroup(.)) %>%
  ungroup()

## Plot for all genes
g1 <- ggplot(tmp %>% filter(complete.cases(gene_label)), aes(x=sample_type, y=ccf, group=plot_id, color = gene_label)) + 
  geom_point(data = tmp %>% filter(!complete.cases(gene_label))) + 
  geom_line(data = tmp %>% filter(!complete.cases(gene_label))) + 
  geom_point() + 
  geom_line(na.rm = TRUE) + 
  geom_text(data = test_subtype, aes(x=1.5, y=0.05, label = test_str), group = NA, color = "black") +
  facet_wrap(~idh_codel_subtype) +
  scale_color_manual(values = c("#B4464B", "#B47846", "#B4AF46", "#82B446", "#4BB446", "#46B478", "#46B4AF", "#4682B4", "#464BB4", "#7846B4", "#AF46B4", "#B44682"), na.value = "gray75") +
  theme_bw(base_size = 12) +
  labs(x = "Sample Type", y = "Cancer Cell Fraction", color = "Gene Symbol")

g1

pdf("~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/shared_ccf_paired_ladderplot.pdf", width=12, height = 8, useDingbats = FALSE)
plot(g1)
dev.off()

## Plot faceted for each gene seperately
g2 <- tmp %>% 
  filter(complete.cases(gene_label), !(gene_symbol %in% c("CIC","FUBP1","TERT","RB1"))) %>%
  ggplot() + 
  geom_point(aes(x=sample_type, y=ccf, group=plot_id, color = variant_classification_vep)) + 
  geom_line(aes(x=sample_type, y=ccf, group=plot_id, color = variant_classification_vep), na.rm = TRUE) + 
  geom_text(data=filter(test_gene, !(gene_symbol %in% c("CIC","FUBP1","TERT","RB1"))), aes(x=1.5, y=0.1, label = test_str)) +
  facet_wrap(~gene_label, ncol = 4) +
  scale_color_manual(values=c("5'Flank" = brewer.pal(9, "Paired")[9],
                              "Frame_Shift_Del" = brewer.pal(7, "Paired")[1],
                              "Frame_Shift_Ins" = brewer.pal(7, "Paired")[2],
                              "In_Frame_Del" = brewer.pal(7, "Paired")[3],
                              "In_Frame_Ins" = brewer.pal(7, "Paired")[4],
                              "Missense_Mutation" = brewer.pal(7, "Paired")[5],
                              "Nonsense_Mutation" = brewer.pal(7, "Paired")[6],
                              "Splice_Site" = brewer.pal(9, "Paired")[7],
                              "Translation_Start_Site" = brewer.pal(9, "Paired")[8])) +
  theme_bw(base_size = 12) +
  labs(x = "Sample Type", y = "Cancer Cell Fraction", color = "Variant Classification") +
  coord_cartesian(ylim = c(0,1))

g2

pdf("~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/shared_ccf_paired_ladderplot_gene.pdf", width=12, height = 6, useDingbats = FALSE)
plot(g2)
dev.off()

## Plot for TP53 between subtypes
g2b <- tmp %>% filter(gene_label=="TP53") %>%
  ggplot() + 
  geom_point(aes(x=sample_type, y=ccf, group=plot_id, color = variant_classification_vep)) + 
  geom_line(aes(x=sample_type, y=ccf, group=plot_id, color = variant_classification_vep), na.rm = TRUE) + 
  geom_text(data=filter(test_gene_subtype, gene_symbol == "TP53"), aes(x=1.5, y=0.1, label = test_str)) +
  facet_wrap(~idh_codel_subtype) +
  scale_color_manual(values=c("5'Flank" = brewer.pal(9, "Paired")[9],
                              "Frame_Shift_Del" = brewer.pal(7, "Paired")[1],
                              "Frame_Shift_Ins" = brewer.pal(7, "Paired")[2],
                              "In_Frame_Del" = brewer.pal(7, "Paired")[3],
                              "In_Frame_Ins" = brewer.pal(7, "Paired")[4],
                              "Missense_Mutation" = brewer.pal(7, "Paired")[5],
                              "Nonsense_Mutation" = brewer.pal(7, "Paired")[6],
                              "Splice_Site" = brewer.pal(9, "Paired")[7],
                              "Translation_Start_Site" = brewer.pal(9, "Paired")[8])) +
  theme_bw(base_size = 12) +
  labs(x = "Sample Type", y = "Cancer Cell Fraction", color = "Variant Classification", title = "TP53") +
  coord_cartesian(ylim = c(0,1))

g2b

pdf("~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/shared_ccf_paired_ladderplot_TP53.pdf", width=6, height = 3, useDingbats = FALSE)
plot(g2b)
dev.off()

## Plot for IDH1 between subtypes
g2c <- tmp %>% filter(gene_label=="IDH1") %>%
  ggplot() + 
  geom_point(aes(x=sample_type, y=ccf, group=plot_id, color = variant_classification_vep)) + 
  geom_line(aes(x=sample_type, y=ccf, group=plot_id, color = variant_classification_vep), na.rm = TRUE) + 
  geom_text(data=filter(test_gene_subtype, gene_symbol == "IDH1"), aes(x=1.5, y=0.1, label = test_str)) +
  facet_wrap(~idh_codel_subtype) +
  scale_color_manual(values=c("5'Flank" = brewer.pal(9, "Paired")[9],
                              "Frame_Shift_Del" = brewer.pal(7, "Paired")[1],
                              "Frame_Shift_Ins" = brewer.pal(7, "Paired")[2],
                              "In_Frame_Del" = brewer.pal(7, "Paired")[3],
                              "In_Frame_Ins" = brewer.pal(7, "Paired")[4],
                              "Missense_Mutation" = brewer.pal(7, "Paired")[5],
                              "Nonsense_Mutation" = brewer.pal(7, "Paired")[6],
                              "Splice_Site" = brewer.pal(9, "Paired")[7],
                              "Translation_Start_Site" = brewer.pal(9, "Paired")[8])) +
  theme_bw(base_size = 12) +
  labs(x = "Sample Type", y = "Cancer Cell Fraction", color = "Variant Classification", title = "IDH1") +
  coord_cartesian(ylim = c(0,1))

g2c

pdf("~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/shared_ccf_paired_ladderplot_IDH1.pdf", width=6, height = 3, useDingbats = FALSE)
plot(g2c)
dev.off()

## Plot CCF wilcoxon results for each patient
test_case <- test_case %>%
  mutate(direction = case_when(wilcox_p < 0.05 & median_b > median_a ~ "CCF increase",
                               wilcox_p < 0.05 & median_b < median_a ~ "CCF decrease",
                               wilcox_p > 0.05 ~ "CCF stable",
                               TRUE ~ "???"))

test_case_counts <- test_case %>% count(idh_codel_subtype,direction) %>%
  group_by(idh_codel_subtype) %>% 
  mutate(n_total=sum(n),
         prop = n/sum(n),
         prop_txt = sprintf("%s%%", round(100*n/sum(n),1))) %>% 
  ungroup()

g2d <- ggplot(test_case, aes(x = idh_codel_subtype, fill = direction)) +
  geom_bar(position = position_dodge(width = 1)) +
  geom_text(data = test_case_counts, aes(label = prop_txt, y = n + 2), position = position_dodge(width = 1)) +
  theme_bw(base_size = 10) +
  labs(x = "Subtype", y = "Number of Patients", fill = "CCF Change") +
  scale_fill_manual(values = c( "#00A3FF", "#5C00FF", "#FFDC00")) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  coord_cartesian(ylim=c(0,65)) +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        axis.line = element_blank()) #+
  #facet_wrap(~idh_codel_subtype)

g2d

pdf("~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/shared_ccf_paired_patient_wilcoxon.pdf", width=6, height = 4, useDingbats = FALSE)
plot(g2d)
dev.off()

############################################################################################################################################
## Plot scatter plot of initial and recurrent CCF
## For all genes w/ colored driver genes
## Faceted by driver gene
############################################################################################################################################

## Munge data
tmp2 <- res %>% filter(rank==1) %>% 
  select(case_barcode, idh_codel_subtype, variant_classification_vep, gene_symbol, cellular_prevalence_a, cellular_prevalence_b) %>%
  mutate(plot_id = paste(case_barcode,gene_symbol),
         gene_label = factor(case_when(gene_symbol == "TP53" ~ "TP53",
                                       gene_symbol == "IDH1" ~ "IDH1",
                                       gene_symbol == "ATRX" ~ "ATRX",
                                       gene_symbol == "PTEN" ~ "PTEN",
                                       gene_symbol == "PIK3CA" ~ "PIK3CA",
                                       gene_symbol == "PIK3R1" ~ "PIK3R1",
                                       gene_symbol == "NF1" ~ "NF1",
                                       gene_symbol == "EGFR" ~ "EGFR",
                                       gene_symbol == "CIC" ~ "CIC",
                                       gene_symbol == "FUBP1" ~ "FUBP1",
                                       gene_symbol == "TERT" ~ "TERT",
                                       gene_symbol == "RB1" ~ "RB1",
                                       TRUE ~ NA_character_)))

## Plot scatter plot facetted per gene
g3 <- tmp2 %>% filter(complete.cases(gene_label)) %>%
  ggplot(aes(x=cellular_prevalence_a, y=cellular_prevalence_b, color = variant_classification_vep)) + 
  geom_point() + 
  facet_wrap(~gene_label) +
  theme_bw(base_size = 12) +
  geom_abline(slope=1, alpha=0.2, linetype=2) +
  labs(x="Initial CCF", y="Recurrence CCF", color = "Variant Classification") +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) + 
  scale_color_manual(values=c("5'Flank" = brewer.pal(9, "Paired")[9],
                              "Frame_Shift_Del" = brewer.pal(7, "Paired")[1],
                              "Frame_Shift_Ins" = brewer.pal(7, "Paired")[2],
                              "In_Frame_Del" = brewer.pal(7, "Paired")[3],
                              "In_Frame_Ins" = brewer.pal(7, "Paired")[4],
                              "Missense_Mutation" = brewer.pal(7, "Paired")[5],
                              "Nonsense_Mutation" = brewer.pal(7, "Paired")[6],
                              "Splice_Site" = brewer.pal(9, "Paired")[7],
                              "Translation_Start_Site" = brewer.pal(9, "Paired")[8]))

g3

g4 <- ggplot(tmp2, aes(x=cellular_prevalence_a, y=cellular_prevalence_b, color = gene_label)) + 
  geom_point() + 
  facet_wrap(~idh_codel_subtype) + 
  facet_wrap(~idh_codel_subtype) +
  scale_color_manual(values = c("#B4464B", "#B47846", "#B4AF46", "#82B446", "#4BB446", "#46B478", "#46B4AF", "#4682B4", "#464BB4", "#7846B4", "#AF46B4", "#B44682"), na.value = "gray75") +
  theme_bw(base_size = 12) +
  labs(x = "Sample Type", y = "Cancer Cell Fraction", color = "Gene Symbol")

g4

############################################################################################################################################
## Plot bar plots counting clonality
## For all genes w/ colored driver genes
## Faceted by driver gene
############################################################################################################################################

tmp3 <- res %>% filter(rank==1) %>% 
  select(case_barcode, idh_codel_subtype, gene_symbol, clonality_a, clonality_b) %>%
  gather(c(clonality_a, clonality_b), key = "sample_type", value = "clonality") %>%
  mutate(sample_type = factor(sample_type, levels = c("clonality_a", "clonality_b"), labels = c("P","R")),
         plot_id = paste(case_barcode,gene_symbol),
         gene_label = factor(case_when(gene_symbol == "TP53" ~ "TP53",
                                       gene_symbol == "IDH1" ~ "IDH1",
                                       gene_symbol == "ATRX" ~ "ATRX",
                                       gene_symbol == "PTEN" ~ "PTEN",
                                       gene_symbol == "PIK3CA" ~ "PIK3CA",
                                       gene_symbol == "NF1" ~ "NF1",
                                       gene_symbol == "EGFR" ~ "EGFR",
                                       TRUE ~ NA_character_))) %>%
  filter(complete.cases(clonality))

g5 <- ggplot(tmp3, aes(x=clonality, fill = sample_type)) + 
  geom_bar(position = "dodge") +
  facet_wrap(~idh_codel_subtype, scales = "free_y") + 
  labs(x = "Clonality", y = "Number of Mutations", fill = "Sample Type") +
  theme_bw(base_size = 12) + 
  scale_fill_manual(values = c("#B47846", "#4682B4"))

g5

tmp4 <- res %>% 
  filter(rank==1, clonality_a != "ND", clonality_b != "ND") %>% 
  select(case_barcode, idh_codel_subtype, gene_symbol, clonality_a, clonality_b) %>%
  mutate(clonality = sprintf("%s-%s", clonality_a, clonality_b)) %>%
  filter(complete.cases(clonality_a, clonality_b))

tmp4_counts <- tmp4 %>% 
  count(idh_codel_subtype, clonality) %>% 
  group_by(idh_codel_subtype) %>% 
  mutate(x = seq(0.75,1.5,0.25), 
         n_total=sum(n),
         prop = n/sum(n),
         prop_txt = sprintf("%s%%", round(100*n/sum(n),1))) %>% 
  ungroup()

g6 <- ggplot(tmp4, aes(x=1, fill=clonality)) + 
  geom_bar(position = position_dodge(width = 1)) +
  geom_text(data = tmp4_counts, aes(label = prop_txt, y = n + 0.02 * n_total), position = position_dodge(width = 1)) +
  facet_wrap(~idh_codel_subtype, scales = "free") + 
  labs(x = "Subtype", y = "Number of Mutations", fill = "Clonality") +
  theme_bw(base_size = 12)  + 
  scale_fill_manual(values = c("#CA6720", "#2ECA20", "#2083CA", "#BC20CA")) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_x_discrete(expand = c(0,0)) +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        axis.line = element_blank()) #+

g6

pdf("~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/shared_frac_clonality_paired.pdf", width=12, height = 8)
g6
dev.off()

