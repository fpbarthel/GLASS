library(DBI)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(egg)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

## Load heatmap data
hmapdata <- dbGetQuery(con, read_file("sql/build_heatmap_data_mutation.sql"))
genedata <- dbGetQuery(con, read_file("sql/build_heatmap_data_mutation_by_gene.sql"))
cnv_data <- dbGetQuery(con, read_file("sql/build_heatmap_data_cnv.sql"))
mutfdata <- dbGetQuery(con, read_file("sql/mutation_freq_private_shared.sql"))

## Common cases
common_cases <- intersect(hmapdata$case_barcode, cnv_data$case_barcode)

## Load IDH-codel subtypes
clindata <- dbGetQuery(con, "SELECT DISTINCT case_barcode, idh_codel_subtype FROM clinical.surgeries WHERE idh_codel_subtype IS NOT NULL")

## Plot

# Generate a df with information about the proportion of unique and shared mutations.
mut_freq_bar = mut_freq_gold %>% 
  mutate(primary_only = (count_a-intersection_ab)/union_ab,
         recurrence_only =  (count_b-intersection_ab)/union_ab,
         shared = intersection_ab/union_ab) %>% 
  select(primary_only, recurrence_only, shared, tumor_pair_barcode, union_ab, idh_codel_subtype) %>% 
  gather(mutation_type, mutation_percent, c(primary_only, recurrence_only, shared), -tumor_pair_barcode, -union_ab, -idh_codel_subtype)
mut_freq_bar$mutation_type = factor(mut_freq_bar$mutation_type, levels = c("recurrence_only", "primary_only", "shared"))

# Choose complementary colors from the colortools package. http://www.sthda.com/english/wiki/the-elements-of-choosing-colors-for-great-data-visualization-in-r
splitComp("#2fb3ca")
tetradic("#005295")

# Order by total number of mutations (union_ab) AND only present for each class.
top_plot = ggplot(mut_freq_bar, aes(x =reorder(tumor_pair_barcode, as.numeric(union_ab)), y=as.numeric(log10(union_ab)))) + geom_bar(stat="identity") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_text(size=rel(0.6))) +
  xlab("") + ylab("Log10(mutation union)") + facet_grid(~idh_codel_subtype, scales="free")

bot_plot = ggplot(mut_freq_bar, aes(x =reorder(tumor_pair_barcode, as.numeric(union_ab)), y=as.numeric(mutation_percent), fill=mutation_type)) + geom_bar(stat="identity") + 
  scale_fill_manual(values=c("#2FB3CA", "#CA2F66", "#CA932F")) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_text(size=rel(0.9))) +
  xlab("") + ylab("% SNVs and Indels") + guides(fill=guide_legend("Mutation")) + facet_grid(~idh_codel_subtype, scales="free") 

# Combine the plot to incorporate both mutation burden and proportions.
plot_grid(top_plot, bot_plot, align = "v", axis= "rl" , nrow = 2, rel_heights = c(1/4, 3/4), label_size = 12)

## Sort
mutdf <- hmapdata %>%
  filter(case_barcode %in% common_cases) %>%
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
  filter(case_barcode %in% common_cases) %>%
  left_join(clindata) %>%
  mutate(case_barcode = factor(case_barcode, levels = levels(mutdf$case_barcode)))

## Mut Freq
mutfdf <- 
  mutfdata %>%
  filter(case_barcode %in% common_cases) %>%
  left_join(clindata) %>%
  mutate(case_barcode = factor(case_barcode, levels = levels(mutdf$case_barcode)))

## genes
genedf <- genedata %>%
  mutate(gene_symbol = factor(gene_symbol, levels = levels(mutdf$gene_symbol)))

## Plot mutations
p1a <- 
  ggplot(mutdf, aes(x=case_barcode, y=gene_symbol)) + 
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
  labs(color = "Variant Classification", fill = "Coverage", shape = "Variant Retention")

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
p2a <-
  ggplot(cnvdf, aes(x=case_barcode, y=gene_symbol)) + 
  geom_tile(fill = "#f7f7f7") +
  geom_point(aes(color = cnv_class, shape = cnv_retention), size = 2) +
  scale_color_manual(values=c("Deletion" = "#0571b0",
                              "Heterozygous" = "#f7f7f7",
                              "Amplification" = "#ca0020")) +
  labs(color = "Variant Classification", fill = "Coverage", shape = "Variant Retention")

p2b <-
  ggplot(cnvdf, aes(y=gene_symbol)) +
  geom_blank()

## Plot elements
plot_grid      = facet_grid(. ~ idh_codel_subtype, scales = "free_x", space = "free")
plot_theme     <- theme_bw()
null_legend <- theme(legend.position = 'none')
null_x      <- theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank())
null_y      <- theme(axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())
bottom_x    <- theme(axis.text.x=element_blank())
null_facet = theme(strip.background = element_blank(),
                   strip.text.x = element_blank())
top_margin    <- theme(plot.margin= unit(c(1, 1, 0.1, 1), "lines")) ## Top, Right, Bottom, Left
middle_margin <- theme(plot.margin= unit(c(0, 1, 0.1, 1), "lines"))
bottom_margin <- theme(plot.margin= unit(c(0, 1, 1, 1), "lines"))

## Bind
g1a <- ggplotGrob(p1a + plot_grid + plot_theme + null_legend + null_x + null_facet + top_margin + theme(axis.title = element_text(size = 10), axis.text = element_text(size=10))) %>% gtable_frame()
g1b <- ggplotGrob(p1b + plot_theme + null_legend + null_x + null_y + null_facet + top_margin + theme(axis.title = element_text(size = 10), axis.text = element_text(size=10))) %>% gtable_frame()
g2a <- ggplotGrob(p2a + plot_grid + plot_theme + null_legend + bottom_x + null_facet + bottom_margin + theme(axis.title = element_text(size = 10), axis.text = element_text(size=10))) %>% gtable_frame()
g2b <- ggplotGrob(p2b + plot_theme + null_legend + null_x + null_y + null_facet + top_margin + theme(axis.title = element_text(size = 10), axis.text = element_text(size=10))) %>% gtable_frame()

## PLot mutations
g1 <- gtable_cbind(g1a, g1b)
panels = g1$layout$t[grep("panel", g1$layout$name)]
g1$widths[panels] <- unit(c(0.1,6), "null")

g2 <- gtable_cbind(g2a, g2b)
panels = g2$layout$t[grep("panel", g2$layout$name)]
g2$widths[panels] <- unit(c(0.1,6), "null")

g = gtable_rbind(g1, g2)
plot(g)
