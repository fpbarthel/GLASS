library(DBI)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(egg)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

########################
## Load heatmap data
########################

hmapdata <- dbGetQuery(con, read_file("sql/build_heatmap_data_mutation.sql"))
genedata <- dbGetQuery(con, read_file("sql/build_heatmap_data_mutation_by_gene.sql"))
cnv_data <- dbGetQuery(con, read_file("sql/build_heatmap_data_cnv.sql"))
mutfdata <- dbGetQuery(con, read_file("sql/mutation_freq_private_shared.sql"))
clindata <- dbGetQuery(con, "SELECT * FROM clinical.surgeries WHERE idh_codel_subtype IS NOT NULL")
subjdata <- dbGetQuery(con, "SELECT * FROM clinical.cases")
anpldata <- dbGetQuery(con, "SELECT * FROM analysis.titan_seg_paired_cn2_comparison")
c710data <- dbGetQuery(con, read_file("sql/compute_chr7_10.sql"))

## Common cases
common_cases <- intersect(hmapdata$case_barcode, cnv_data$case_barcode)

########################
# Data pre-processing: transform data into a form that can be ingested and used by ggplot plotting engine
########################

idh_codel_data <- clindata %>% 
  select(case_barcode, idh_codel_subtype) %>% 
  distinct() %>%
  filter(case_barcode %in% common_cases)

c710_data <- c710data %>%
  filter(case_barcode %in% common_cases) %>%
  left_join(idh_codel_data)

## Time data
surv_data <- subjdata %>%
  filter(case_barcode %in% common_cases) %>%
  select(case_barcode, surgery_number = case_vital_status, surgical_interval_mo = case_overall_survival_mo) %>%
  left_join(idh_codel_data)

time_data <- clindata %>%
  filter(case_barcode %in% common_cases) %>%
  mutate(surgery_number = as.character(surgery_number)) %>%
  bind_rows(surv_data) %>%
  mutate(event_type = ifelse(surgery_number %in% as.character(1:5), "surgery", surgery_number))

anpl_data <- anpldata %>% 
  mutate(case_barcode = substr(tumor_pair_barcode, 1, 12)) %>%
  filter(case_barcode %in% common_cases, !duplicated(case_barcode)) %>%
  left_join(idh_codel_data)# %>%
  #mutate(prop_het_a = ifelse(diamond_set==1,prop_het_a,NA),
  #       prop_het_b = ifelse(diamond_set==1,prop_het_b,NA))

## Generate a data frame with information about the proportion of unique and shared mutations.
mut_freq_case <- mutfdata %>%
  filter(case_barcode %in% common_cases) %>%
  left_join(idh_codel_data) %>%
  mutate(mfchange = factor(ifelse(mf_b > mf_a, "+", ifelse(mf_a > mf_b, "-", NA))))

mut_freq_prop_case = mutfdata %>% 
  filter(case_barcode %in% common_cases) %>%
  mutate(P = (count_a-intersection_ab)/union_ab,
         R =  (count_b-intersection_ab)/union_ab,
         S = intersection_ab/union_ab) %>% 
  select(P, R, S, case_barcode, union_ab, idh_codel_subtype) %>% 
  gather(mutation_type, mutation_percent, c(P, R, S), -case_barcode, -union_ab, -idh_codel_subtype) %>%
  mutate(mutation_type = factor(mutation_type, levels = c("R", "P", "S")))

mut_freq_prop_gene = genedata %>% 
  mutate(P = (count_a-shared)/total,
         R =  (count_b-shared)/total,
         S = shared/total) %>% 
  select(P, R, S, gene_symbol) %>% 
  gather(mutation_type, mutation_percent, c(P, R, S), -gene_symbol) %>%
  mutate(mutation_type = factor(mutation_type, levels = c("R", "P", "S")))

## Pre-process mutation data

mut_heatmap <- hmapdata %>%
  complete(gene_symbol,case_barcode) %>%
  left_join(idh_codel_data) %>%
  filter(case_barcode %in% common_cases) %>%
  group_by(gene_symbol) %>% 
  mutate(n_mut_gene = sum(!is.na(variant_call))) %>% 
  ungroup() %>%
  arrange(n_mut_gene) %>%
  mutate(gene_symbol = factor(gene_symbol, levels = unique(gene_symbol))) %>%
  group_by(case_barcode) %>%
  mutate(n_mut_case = sum(!is.na(variant_call))) %>% 
  ungroup() %>%
  arrange(desc(n_mut_case)) %>%
  mutate(covered = ifelse(is.na(covered), "Coverage < 5", covered)) %>%
  mutate(case_barcode = factor(case_barcode, levels = unique(case_barcode)))

## Pre-process copy number data 

cnv_heatmap <- cnv_data %>%
  filter(case_barcode %in% common_cases) %>%
  left_join(idh_codel_data)

########################
# Data re-ordering: order plot elements to match one another
########################

case_order <- levels(mut_heatmap$case_barcode)
gene_order <- levels(mut_heatmap$gene_symbol)

idh_codel_data <- idh_codel_data %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
c710_data <- c710_data %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
time_data <- time_data %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
mut_freq_case <- mut_freq_case %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
mut_freq_prop_case <- mut_freq_prop_case %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
cnv_heatmap <- cnv_heatmap %>% mutate(case_barcode = factor(case_barcode, levels = case_order))

mut_freq_prop_gene <- mut_freq_prop_gene %>% mutate(gene_symbol = factor(gene_symbol, levels = gene_order))

######################## 
## Common plotting elements
########################

plot_grid     <- facet_grid(. ~ idh_codel_subtype, scales = "free_x", space = "free")
plot_theme    <- theme_bw() + theme(axis.title = element_text(size = 10), axis.text = element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
null_legend   <- theme(legend.position = 'none')
null_x        <- theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank())
null_y        <- theme(axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())
bottom_x      <- theme(axis.text.x=element_blank())
null_facet    <- theme(strip.background = element_blank(),
                   strip.text.x = element_blank())
top_margin    <- theme(plot.margin= unit(c(1, 1, 0.1, 1), "lines")) ## Top, Right, Bottom, Left
middle_margin <- theme(plot.margin= unit(c(0, 1, 0.1, 1), "lines"))
bottom_margin <- theme(plot.margin= unit(c(0, 1, 1, 1), "lines"))

########################
## View test plot
########################

testPlot <- function(gg, grid = TRUE) {
  if(grid)
    gg + plot_theme + plot_grid
  else
    gg + plot_theme
}

gg_cbind <- function(..., widths = NULL) {
  if(length(match.call()) - 2 != length(widths))
    message("Number of widths does not match number of columns")
  gg <- gtable_cbind(...)
  panels = gg$layout$t[grep("panel", gg$layout$name)]
  gg$widths[panels] <- unit(widths, "null")
  return(gg)
}


gg_rbind <- function(..., heights = NULL, ncol = 2) {
  if(length(match.call()) - 3 != length(heights))
    message("Number of heights does not match number of rows")
  gg <- gtable_rbind(...)
  panels = gg$layout$t[grep("panel", gg$layout$name)]
  gg$heights[panels] <- unit(rep(heights,each = ncol), "null")
  return(gg)
}

## Extract legend
gg_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

########################
## Blank plot
########################

gg_blank <-
  ggplot(cnvdf, aes(y=gene_symbol)) +
  geom_blank()

########################
## Plot mutation frequencies
########################

gg_mut_freq_case <-
  ggplot(mut_freq_case, aes(x=case_barcode)) +
  geom_hline(yintercept = 10, alpha=0.8, linetype=2) +
  geom_point(aes(y=mf_a), color = "#a6611a") +
  geom_linerange(aes(ymin = mf_a, ymax = mf_b, color = mfchange), linetype = 2) +
  geom_point(aes(y=mf_b), color = "#018571") +
  scale_color_manual(values = c("+" = "#bebada", "-" = "#fb8072")) +
  scale_y_log10() +
  labs(y = "Mutation Frequency (mutations / Mb)") 

testPlot(gg_mut_freq_case)


########################
## Plot clinical data
########################

gg_time_data <-
  ggplot(time_data, aes(x=case_barcode)) +
  geom_line(aes(y = surgical_interval_mo), linetype = 1, alpha = 0.5, color = "gray") +
  geom_point(aes(y=surgical_interval_mo, color=factor(surgery_number), shape = factor(event_type), fill= factor(surgery_number))) +
  scale_shape_manual(values = c("surgery" = 1, "alive" = 24, "dead" = 25)) +
  #scale_color_manual(values = c("+" = "#bebada", "-" = "#fb8072")) +
  labs(y = "Surgical Interval (mo)", color = "Surgery Number") 

testPlot(gg_time_data)


########################
## Plot co-variates
########################

gg_idh_codel <-
  ggplot(idh_codel_data, aes(x=case_barcode)) +
  geom_tile(aes(fill = idh_codel_subtype, y = 1)) +
  scale_fill_manual(values = c("#1b9e77","#d95f02", "#7570b3")) 

testPlot(gg_idh_codel)

gg_c710 <-
  ggplot(c710_data, aes(x=case_barcode)) +
  geom_tile(aes(fill = chr7_10_status, y = 1)) +
  scale_fill_manual(values = c("acquired" = "#fc8d59","shared" = "#ffffbf", "shed" = "#91bfdb", "no" = "white")) 

testPlot(gg_c710)

########################
## Plot aneuploidy
########################

gg_anpl_data <-
  ggplot(anpl_data, aes(x=case_barcode)) +
  geom_point(aes(y = prop_het_a), color="red") + 
  geom_point(aes(y = prop_het_b), color="blue") +
  geom_linerange(aes(ymin = prop_het_a, ymax = prop_het_b), linetype = 2)
  #geom_bar(aes(y=prop_delta_eq), stat="identity")
  
  #geom_line(aes(y = surgical_interval_mo), linetype = 1, alpha = 0.5, color = "gray") +
  #geom_point(aes(y=surgical_interval_mo, color=factor(surgery_number), shape = factor(event_type), fill= factor(surgery_number))) +
  #scale_shape_manual(values = c("surgery" = 1, "alive" = 24, "dead" = 25)) +
  #scale_color_manual(values = c("+" = "#bebada", "-" = "#fb8072")) +
  #labs(y = "Surgical Interval (mo)", color = "Surgery Number") 

testPlot(gg_anpl_data)

########################
## Plot proportions shared/private per patient
########################

gg_mut_freq_prop_case <-
  ggplot(mut_freq_prop_case, aes(x = case_barcode, y = mutation_percent, fill=mutation_type)) +
  geom_bar(stat="identity") + 
  labs(y = "% SNVs and Indels") +
  scale_fill_manual(values=c("#2FB3CA", "#CA2F66", "#CA932F")) +
  plot_grid

testPlot(gg_mut_freq_prop_case)

########################
## Plot proportions shared/private per gene
########################

gg_mut_freq_prop_gene <-
  ggplot(mut_freq_prop_gene, aes(x = gene_symbol, y = mutation_percent, fill=mutation_type)) +
  geom_bar(stat="identity") + 
  labs(y = "% SNVs and Indels") +
  scale_fill_manual(values=c("#2FB3CA", "#CA2F66", "#CA932F")) +
  coord_flip()

testPlot(gg_mut_freq_prop_gene, grid = FALSE)

########################
## Plot mutations
########################

gg_mut_heatmap <- 
  ggplot(mut_heatmap, aes(x=case_barcode, y=gene_symbol)) + 
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
  labs(y = "Gene Symbol", color = "Variant Classification", fill = "Coverage", shape = "Variant Retention")

testPlot(gg_mut_heatmap)

########################
## Plot CNV
########################

gg_cnv_heatmap <-
  ggplot(cnv_heatmap, aes(x=case_barcode, y=gene_symbol)) + 
  geom_tile(fill = "#f7f7f7") +
  geom_point(aes(color = cnv_class, shape = cnv_retention), size = 2) +
  scale_color_manual(values=c("Deletion" = "#0571b0",
                              "Heterozygous" = "#f7f7f7",
                              "Amplification" = "#ca0020")) +
  labs(y = "Gene Symbol", color = "Variant Classification", fill = "Coverage", shape = "Variant Retention", x = "Patient")

########################
## Layout plot elements
########################

g1_left <- ggplotGrob(gg_mut_freq_case + plot_grid + plot_theme + null_legend + null_x + null_facet + top_margin) %>% gtable_frame()
g1_right <- ggplotGrob(gg_blank + plot_theme + null_legend + null_x + null_y + null_facet + top_margin) %>% gtable_frame()

g2_left <- ggplotGrob(gg_time_data + plot_grid + plot_theme + null_legend + null_x + null_facet + top_margin) %>% gtable_frame()
g2_right <- ggplotGrob(gg_blank + plot_theme + null_legend + null_x + null_y + null_facet + top_margin) %>% gtable_frame()

g3_left <- ggplotGrob(gg_idh_codel + plot_grid + plot_theme + null_legend + null_x + null_facet + top_margin) %>% gtable_frame()
g3_right <- ggplotGrob(gg_blank + plot_theme + null_legend + null_x + null_y + null_facet + top_margin) %>% gtable_frame()

g4_left <- ggplotGrob(gg_c710 + plot_grid + plot_theme + null_legend + null_x + null_facet + top_margin) %>% gtable_frame()
g4_right <- ggplotGrob(gg_blank + plot_theme + null_legend + null_x + null_y + null_facet + top_margin) %>% gtable_frame()

g5_left <- ggplotGrob(gg_mut_freq_prop_case + plot_grid + plot_theme + null_legend + null_x + null_facet + bottom_margin) %>% gtable_frame()
g5_right <- ggplotGrob(gg_blank + plot_theme + null_legend + null_x + null_y + null_facet + top_margin) %>% gtable_frame()

g6_left <- ggplotGrob(gg_mut_heatmap + plot_grid + plot_theme + null_legend + null_x + null_facet + bottom_margin) %>% gtable_frame()
g6_right <- ggplotGrob(gg_mut_freq_prop_gene + plot_theme + null_legend + null_x + null_y + null_facet + top_margin) %>% gtable_frame()

g7_left <- ggplotGrob(gg_cnv_heatmap + plot_grid + plot_theme + null_legend + bottom_x + null_facet + bottom_margin) %>% gtable_frame()
g7_right <- ggplotGrob(gg_blank + plot_theme + null_legend + null_x + null_y + null_facet + top_margin) %>% gtable_frame()

g1 <- gg_cbind(g1_left, g1_right, widths = c(0.1,8))
g2 <- gg_cbind(g2_left, g2_right, widths = c(0.1,8))
g3 <- gg_cbind(g3_left, g3_right, widths = c(0.1,8))
g4 <- gg_cbind(g4_left, g4_right, widths = c(0.1,8))
g5 <- gg_cbind(g5_left, g5_right, widths = c(0.1,8))
g6 <- gg_cbind(g6_left, g6_right, widths = c(0.1,8))
g7 <- gg_cbind(g7_left, g7_right, widths = c(0.1,8))

g <- gg_rbind(g1, g2, g3, g4, g5, g6, g7, heights = c(1, 1, 0.25, 0.25, 1, 4, 2), ncol = 2)
plot(g)

gleg1 <- gg_rbind(gg_legend(gg_mut_freq_case),
                 gg_legend(gg_time_data),
                 gg_legend(gg_mut_freq_prop_case),
                 heights = rep(1,1),
                 ncol = 1)

gleg2 <- gg_rbind(gg_legend(gg_mut_heatmap),
                  gg_legend(gg_cnv_heatmap),
                  heights = rep(1,1),
                  ncol = 1)

gleg3 <- gg_rbind(gg_legend(gg_idh_codel),
                  gg_legend(gg_c710),
                  heights = rep(1,1),
                  ncol = 1)
plot(gleg1)
plot(gleg2)
plot(gleg3)
