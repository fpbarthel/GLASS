library(DBI)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

## Load heatmap data
hmapdata <- dbGetQuery(con, read_file("sql/build_heatmap_data_mutation.sql"))
genedata <- dbGetQuery(con, read_file("sql/build_heatmap_data_mutation_by_gene.sql"))
cnv_data <- dbGetQuery(con, read_file("sql/build_heatmap_data_cnv.sql"))
## Load IDH-codel subtypes
clindata <- dbGetQuery(con, "SELECT DISTINCT case_barcode, idh_codel_subtype FROM clinical.surgeries WHERE idh_codel_subtype IS NOT NULL")

## Sort
df <- hmapdata %>% 
  complete(gene_symbol,case_barcode) %>%
  left_join(clindata) %>%
  group_by(gene_symbol) %>% 
  mutate(n_mut_gene = sum(!is.na(variant_call))) %>% 
  ungroup() %>%
  arrange(n_mut_gene) %>%
  mutate(gene_symbol = factor(gene_symbol, levels = unique(gene_symbol))) %>%
  group_by(case_barcode) %>%
  mutate(n_mut_case = sum(!is.na(variant_call))) %>% 
  ungroup() %>%
  arrange(desc(n_mut_case)) %>%
  mutate(case_barcode = factor(case_barcode, levels = unique(case_barcode))) %>%
  mutate(covered = ifelse(is.na(covered), "Coverage < 5", covered))

## CNV 
cnvdf <- cnv_data %>%
  left_join(clindata)

ggplot(df, aes(x=case_barcode, y=gene_symbol)) + 
  geom_tile(aes(fill = covered)) + 
  geom_point(aes(shape = variant_call, color = variant_classification)) +
  facet_grid(. ~ idh_codel_subtype, scales = "free_x", space = "free") +
  scale_fill_manual(values=c("gray90","white"), na.value = "gray90") + 
  scale_color_manual(values=c("5'Flank" = brewer.pal(9, "Paired")[9],
                             "Frame_Shift_Del" = brewer.pal(7, "Paired")[1],
                             "Frame_Shift_Ins" = brewer.pal(7, "Paired")[2],
                             "In_Frame_Del" = brewer.pal(7, "Paired")[3],
                             "In_Frame_Ins" = brewer.pal(7, "Paired")[4],
                             "Missense_Mutation" = brewer.pal(7, "Paired")[5],
                             "Nonsense_Mutation" = brewer.pal(7, "Paired")[6],
                             "Splice_Site" = brewer.pal(9, "Paired")[7],
                             "Translation_Start_Site" = brewer.pal(9, "Paired")[8])) +
  labs(color = "Variant Classification", fill = "Coverage", shape = "Variant Retention") +
  theme_bw()

genedf <- genedata %>%
  mutate(gene_symbol = factor(gene_symbol, levels = levels(df$gene_symbol)))

ggplot(genedf, aes(x=gene_symbol)) +
  geom_bar(aes(y=count_a), stat="identity") + 
  geom_bar(aes(y=shared), stat="identity", fill = "yellow") + 
  geom_bar(aes(y=-count_b), stat="identity") +
  geom_bar(aes(y=-shared), stat="identity", fill = "yellow") +
  geom_hline(yintercept = 0, linetype = 2) +
  coord_flip() +
  theme_bw() + 
  labs(y = "Number of Variants\n<-- recurrence -- primary -->")


## plot cnv
ggplot(cnvdf, aes(x=case_barcode, y=gene_symbol)) + 
  geom_tile(aes(fill = cnv_class)) +
  geom_point(aes(shape = cnv_call), color = "#f7f7f7") +
  facet_grid(. ~ idh_codel_subtype, scales = "free_x", space = "free") +
  scale_fill_manual(values=c("Homozygous deletion" = "#0571b0",
                              "Deletion" = "#92c5de",
                              "Heterozygous" = "#f7f7f7",
                              "Amplification" = "#f4a582",
                              "High-level amplification" = "#ca0020")) +
  labs(color = "Variant Classification", fill = "Coverage", shape = "Variant Retention") +
  theme_bw()
