library(DBI)
library(tidyverse)
library(ggplot2)

con <- DBI::dbConnect(odbc::odbc(), "GLASSv2")  
res <- dbGetQuery(con, read_file("sql/pyclone/pyclone_cluster_stats2.sql"))

tmp <- res %>% 
  mutate(known_driver = ifelse(is.na(num_drivers), "No known driver", "Known driver mutation"),
         plot_id = paste(case_barcode,cluster_id))

#wheel("#B4464B",7)
#test = tmp %>% count(case_barcode, gene_symbol)

g1 <- ggplot(tmp, aes(x=sample_type, y=mean, group=plot_id, color = known_driver)) + 
  geom_point() + 
  geom_line(na.rm = TRUE) + 
  facet_wrap(~idh_codel_subtype) + 
  scale_color_manual(values = c("#2FB3CA","#F1564F")) +
  theme_bw(base_size = 12) +
  labs(x = "Sample Type", y = "Cancer Cell Fraction", color = "Gene Symbol")

g1

pdf("figures/pyclone_cluster_ccf_paired_ladderplot.pdf", width=12, height = 8)
g1
dev.off()


tmp <- res %>%
  group_by(case_barcode, idh_codel_subtype) %>%
  summarize(n_cluster = n_distinct(cluster_id)) %>%
  ungroup()

g2 <- ggplot(tmp, aes(x=idh_codel_subtype, y = n_cluster)) + geom_boxplot()

g2

test_subgroup <- function(df) {
  wtest = wilcox.test(df$mean ~ df$sample_type, paired = TRUE)
  data.frame(n = nrow(df), median_a = median(df$mean[df$sample_type == "P"]), median_b = median(df$mean[df$sample_type == "R"]), wilcox_p = wtest$p.value, wilcox_v = wtest$statistic)
}

tmp %>% group_by(idh_codel_subtype) %>% do(test_subgroup(.))
tmp %>% group_by(idh_codel_subtype,known_driver) %>% do(test_subgroup(.))


### PLOT barplot of clonal/subclobal

tmp <- res %>% filter(rank==1) %>% 
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

g2 <- ggplot(tmp, aes(x=clonality, fill = sample_type)) + 
  geom_bar(position = "dodge") +
  facet_wrap(~idh_codel_subtype, scales = "free_y") + 
  labs(x = "Clonality", y = "Number of Mutations", fill = "Sample Type") +
  theme_bw(base_size = 12) + 
  scale_fill_manual(values = c("#B47846", "#4682B4"))

pdf("figures/shared_frac_clonality.pdf", width=12, height = 8)
g2
dev.off()

tmp <- res %>% filter(rank==1) %>% 
  select(case_barcode, idh_codel_subtype, gene_symbol, clonality_a, clonality_b) %>%
  mutate(clonality = sprintf("%s-%s", clonality_a, clonality_b)) %>%
  filter(complete.cases(clonality_a, clonality_b))

g3 <- ggplot(tmp, aes(x=1, fill=clonality)) + 
  geom_bar() +
  facet_wrap(~idh_codel_subtype, scales = "free") + 
  labs(x = "Subtype", y = "Number of Mutations", fill = "Clonality") +
  theme_bw(base_size = 12) 

g3

pdf("figures/shared_frac_clonality_paired.pdf", width=12, height = 8)
g3
dev.off()

