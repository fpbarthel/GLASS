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

## genes
genedf <- genedata %>%
  mutate(gene_symbol = factor(gene_symbol, levels = levels(df$gene_symbol)))

## Plot mutations
p1a <- 
  ggplot(df, aes(x=case_barcode, y=gene_symbol)) + 
  geom_tile(aes(fill = covered)) + 
  geom_point(aes(shape = variant_call, color = variant_classification)) +
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

## Plot genes
p1b <-
  ggplot(genedf, aes(x=gene_symbol)) +
  geom_bar(aes(y=count_a), stat="identity") + 
  geom_bar(aes(y=shared), stat="identity", fill = "yellow") + 
  geom_bar(aes(y=-count_b), stat="identity") +
  geom_bar(aes(y=-shared), stat="identity", fill = "yellow") +
  geom_hline(yintercept = 0, linetype = 2) +
  coord_flip() +
  labs(y = "Number of Variants\n<-- recurrence -- primary -->")


## plot cnv
p2 <-
  ggplot(cnvdf, aes(x=case_barcode, y=gene_symbol)) + 
  geom_tile(fill = "#f7f7f7") +
  geom_point(aes(color = cnv_class, shape = cnv_retention), size = 2) +
  scale_color_manual(values=c("Deletion" = "#0571b0",
                              "Heterozygous" = "#f7f7f7",
                              "Amplification" = "#ca0020")) +
  labs(color = "Variant Classification", fill = "Coverage", shape = "Variant Retention")

## Plot elements
plot_grid      = facet_grid(. ~ idh_codel_subtype, scales = "free_x", space = "free")
plot_theme     <- theme_bw()
null_legend <- theme(legend.position = 'none')
null_x      <- theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank()) 
bottom_x    <- theme(axis.text.x=element_blank())
null_facet = theme(strip.background = element_blank(),
                   strip.text.x = element_blank())
top_margin    <- theme(plot.margin= unit(c(1, 1, 0.1, 1), "lines")) ## Top, Right, Bottom, Left
middle_margin <- theme(plot.margin= unit(c(0, 1, 0.1, 1), "lines"))
bottom_margin <- theme(plot.margin= unit(c(0, 1, 1, 1), "lines"))

## Bind
g1 = ggplotGrob(p1a + plot_grid + plot_theme + null_legend + null_x + null_facet + top_margin + theme(axis.title = element_text(size = 10), axis.text = element_text(size=10)))  %>% gtable_frame()
g2 = ggplotGrob(g2 + plot_grid + plot_theme + null_legend + bottom_x + null_facet + bottom_margin + theme(axis.title = element_text(size = 10), axis.text = element_text(size=10)))  %>% gtable_frame()

g = gtable_rbind(g1, g2)
