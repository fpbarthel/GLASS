library(DBI)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(egg)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

########################
## Load heatmap data
########################

hmapdata <- dbGetQuery(con, read_file("sql/heatmap/heatmap_snv.sql"))
snvgdata <- dbGetQuery(con, read_file("sql/heatmap/heatmap_snv_by_gene.sql"))
cnv_data <- dbGetQuery(con, read_file("sql/heatmap/heatmap_cnv.sql"))
cnvgdata <- dbGetQuery(con, read_file("sql/heatmap/heatmap_cnv_by_gene.sql"))

mutfdata <- dbGetQuery(con, read_file("sql/heatmap/heatmap_mf.sql"))
clindata <- dbGetQuery(con, read_file("sql/heatmap/heatmap_clinical.sql"))
timedata <- dbGetQuery(con, read_file("sql/heatmap/heatmap_time.sql"))
anpldata <- dbGetQuery(con, read_file("sql/heatmap/heatmap_aneuploidy.sql"))
c710data <- dbGetQuery(con, read_file("sql/heatmap/heatmap_c710.sql"))
puridata <- dbGetQuery(con, read_file("sql/heatmap/heatmap_purity.sql"))
drivdata <- dbGetQuery(con, read_file("sql/heatmap/heatmap_drivers.sql"))

neutrdata = read_tsv("sandbox/glass-evolution-neutrality-20190129.txt")
all_drivers = dbGetQuery(con, "SELECT * FROM analysis.driver_status")
neutrdata = all_drivers %>% 
  left_join(neutrdata, by="tumor_pair_barcode") %>% 
  mutate(case_barcode = substr(tumor_pair_barcode, 1, 12),
         evolution = paste(bayesian_evo_p,bayesian_evo_r, sep="-" ),
         evolution = ifelse(evolution=="NA-NA", NA, evolution)) %>% 
  select(case_barcode, evolution, idh_codel_subtype = idh_codel_subtype.x)

########################
# Data pre-processing: transform data into a form that can be ingested and used by ggplot plotting engine
########################

time_data <- timedata %>%
  mutate(event_type = ifelse(surgery_number %in% as.character(1:5), "surgery", surgery_number))

## Generate a data frame with information about the proportion of unique and shared mutations.
mut_freq_case <- mutfdata %>%
  mutate(mfchange = factor(ifelse(mf_b > mf_a, "+", ifelse(mf_a > mf_b, "-", NA)))) %>%
  arrange(desc(mf_b))

mut_freq_prop_case = mutfdata %>% 
  mutate(P = (count_a-intersection_ab)/union_ab,
         R =  (count_b-intersection_ab)/union_ab,
         S = intersection_ab/union_ab) %>% 
  select(P, R, S, case_barcode, union_ab, idh_codel_subtype) %>% 
  gather(mutation_type, mutation_percent, c(P, R, S), -case_barcode, -union_ab, -idh_codel_subtype) %>%
  mutate(mutation_type = factor(mutation_type, levels = c("R", "P", "S")))

sort_df <-
  mut_freq_case %>%
  left_join(drivdata) %>%
  left_join(clindata) %>% 
  left_join(anpldata) %>%
  arrange(desc(mf_b)) # snv_driver_count, cnv_driver_count)#, aneuploidy_b - aneuploidy_a, 

#mut_freq_prop_gene = genedata %>% 
#  mutate(P = (count_a-shared)/total,
#         R =  (count_b-shared)/total,
#         S = shared/total) %>% 
#  select(P, R, S, gene_symbol) %>% 
#  gather(mutation_type, mutation_percent, c(P, R, S), -gene_symbol) %>%
#  mutate(mutation_type = factor(mutation_type, levels = c("R", "P", "S")))

########################
# Data re-ordering: order plot elements to match one another
########################

case_order <- unique(sort_df$case_barcode)
snv_gene_order <- rev(unique(snvgdata$gene_symbol))
cnv_gene_order <- rev(unique(cnvgdata$gene_symbol))

#idh_codel_data <- idh_codel_data %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
c710_data <- c710data %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
time_data <- time_data %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
clin_data <- clindata %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
mut_freq_case <- mut_freq_case %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
anpl_data <- anpldata %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
puri_data <- puridata %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
mut_freq_prop_case <- mut_freq_prop_case %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
cnv_heatmap <- cnv_data %>% mutate(case_barcode = factor(case_barcode, levels = case_order),
                                   gene_symbol = factor(gene_symbol, levels = cnv_gene_order))
mut_heatmap <- hmapdata %>% mutate(case_barcode = factor(case_barcode, levels = case_order),
                                   gene_symbol = factor(gene_symbol, levels = snv_gene_order))
driv_hm <- drivdata %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
snvgdata <- snvgdata %>% mutate(gene_symbol = factor(gene_symbol, levels = snv_gene_order))
cnvgdata <- cnvgdata %>% mutate(gene_symbol = factor(gene_symbol, levels = cnv_gene_order))
neutrdata <- neutrdata %>% mutate(case_barcode = factor(case_barcode, levels = case_order))

#mut_freq_prop_gene <- mut_freq_prop_gene %>% mutate(gene_symbol = factor(gene_symbol, levels = gene_order))

######################## 
## Common plotting elements
########################

plot_grid     <- facet_grid(. ~ idh_codel_subtype, scales = "free_x", space = "free")
plot_theme    <- theme_bw(base_size = 10) + theme(axis.title = element_text(size = 10),
                                                  axis.text = element_text(size=10),
                                                  panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(),
                                                  panel.background = element_rect(fill = "transparent"),
                                                  axis.line = element_blank())
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
  ggplot(data.frame()) +
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
  labs(y = "Mut. Freq.\n(mutations / Mb)") 

testPlot(gg_mut_freq_case)


########################
## Plot time/interval data
########################

gg_time_data <-
  ggplot(time_data, aes(x=case_barcode)) +
  geom_line(aes(y = time_mo), linetype = 1, alpha = 0.5, color = "gray") +
  geom_point(aes(y=time_mo, color=factor(surgery_number), shape = factor(event_type), fill= factor(surgery_number))) +
  scale_shape_manual(values = c("surgery" = 1, "alive" = 24, "dead" = 25)) +
  #scale_color_manual(values = c("+" = "#bebada", "-" = "#fb8072")) +
  labs(y = "Surg. Int. (mo)", color = "Surgery Number") 

testPlot(gg_time_data)

########################
## Plot clinical
########################

gg_clinical <-
  clin_data %>% 
  gather(key = "type", value = "value", location_distal, grade_change, received_tmz, received_rt, is_hypermutator) %>%
  mutate(type = factor(type,
                       levels = c("location_distal", "grade_change", "received_tmz", "received_rt", "is_hypermutator"),
                       labels = c("Distal location", "Grade increase", "Received TMZ", "Received RT", "Hypermutator"))) %>%
  ggplot(aes(x=case_barcode)) +
  geom_tile(aes(fill = factor(value), y = type)) +
  scale_fill_manual(values=c("white", "#377eb8"), na.value = "gray90") +
  labs(y="", fill = "Event")

testPlot(gg_clinical)


########################
## Plot co-variates
########################
gg_evolution <-
  ggplot(neutrdata, aes(x=case_barcode)) +
  geom_tile(aes(fill = evolution, y = 1)) +
  scale_fill_manual(values = c("N-N" = "#67A3BD", "N-S" = "#8167BD", "S-N" = "#A3BD67", "S-S" = "#BD8167")) +
  labs(y = "Evolution \nmode")

testPlot(gg_evolution)

gg_c710 <-
  ggplot(c710_data, aes(x=case_barcode)) +
  geom_tile(aes(fill = c710_status, y = 1)) +
  scale_fill_manual(values = c("R" = "#fc8d59","S" = "#ffffbf", "P" = "#91bfdb", "WT" = "white"), na.value = "gray90") +
  labs(y="") + null_y

testPlot(gg_c710)

########################
## Plot aneuploidy
########################

gg_anpl_data <-
  ggplot(anpl_data, aes(x=case_barcode)) +
  geom_bar(aes(y = aneuploidy_b - aneuploidy_a), fill = "#beaed4", alpha = 1, stat = "identity") +
  #geom_point(aes(y = 0, color = factor(qc_fail))) +
  coord_cartesian(ylim=c(-1,1)) +
  labs(y = "Aneuploidy \ndifference")
  #scale_color_manual(values = c("1" = "red", "0" = NA))

testPlot(gg_anpl_data)

# gg_anpl_data1 <-
#   ggplot(anpl_data, aes(x=case_barcode)) +
#   geom_bar(aes(y = aneuploidy_amp_score_b - aneuploidy_amp_score_a), fill = "#de2d26", stat = "identity")
# 
# testPlot(gg_anpl_data1)
# 
# gg_anpl_data2 <-
#   ggplot(anpl_data, aes(x=case_barcode)) +
#   geom_bar(aes(y = aneuploidy_del_score_b - aneuploidy_del_score_a), fill = "#3182bd", stat = "identity")
# 
# testPlot(gg_anpl_data2)
# 
# gg_puri_data <-
#   ggplot(puri_data, aes(x=case_barcode)) +
#   geom_bar(aes(y=purity_b-purity_a), fill = "#a1d99b", stat = "identity")
# 
# testPlot(gg_puri_data)

########################
## Plot driver data
########################

# gg_snv_driver_count <- 
#   driv_hm %>% 
#   gather(key = "type", value = "value", snv_driver_count_shared, snv_driver_count_private_a, snv_driver_count_private_b) %>%
#   ggplot(aes(x=case_barcode)) +
#   geom_bar(aes(y=as.integer(value), fill = type), stat = "identity") + 
#   labs(y = "# of drivers") #+
#   #scale_fill_manual(values = c("#ffff33", "#ff7f00"))
# 
# testPlot(gg_snv_driver_count)
# 
# gg_cnv_driver_count <- 
#   driv_hm %>% 
#   gather(key = "type", value = "value", cnv_driver_count_shared, cnv_driver_count_private_a, cnv_driver_count_private_b) %>%
#   ggplot(aes(x=case_barcode)) +
#   geom_bar(aes(y=as.integer(value), fill = type), stat = "identity") + 
#   labs(y = "# of drivers") #+
# #scale_fill_manual(values = c("#ffff33", "#ff7f00"))

# testPlot(gg_cnv_driver_count)

gg_driver_count <- 
  driv_hm %>% 
  gather(key = "type", value = "value", snv_driver_count_shared, snv_driver_count_private_a, snv_driver_count_private_b, cnv_driver_count_shared, cnv_driver_count_private_a, cnv_driver_count_private_b) %>%
  mutate(type = factor(type,
                       levels = c("snv_driver_count_shared", "cnv_driver_count_shared", "snv_driver_count_private_a", "cnv_driver_count_private_a", "snv_driver_count_private_b", "cnv_driver_count_private_b"),
                       labels = c("S - SNV", "S - CNV", "P - SNV", "P - CNV", "R - SNV", "R - CNV"))) %>%
  ggplot(aes(x=case_barcode)) +
  geom_bar(aes(y=as.integer(value), fill = type), stat = "identity") + 
  labs(y = "# of drivers") +
  scale_fill_brewer(palette = "Paired", type = "qual")
#scale_fill_manual(values = c("#ffff33", "#ff7f00"))

testPlot(gg_driver_count)

gg_driver_context <-
  driv_hm %>% 
  gather(key = "type", value = "value", snv_driver_evolution, snv_driver_stability, cnv_driver_stability) %>%
  mutate(type = factor(type,
                       levels = c("snv_driver_evolution", "cnv_driver_stability", "snv_driver_stability"),
                       labels = c("Convergence - SNV", "Driver stability - CNV", "Driver stability - SNV"))) %>%
  ggplot(aes(x=case_barcode)) +
  geom_tile(aes(fill = value, y = type)) +
  scale_fill_manual(values=c("#007c80", "#7b0a6b", "#ce8014"), na.value = "gray95") +
  labs(y="")

testPlot(gg_driver_context)

########################
## Plot proportions shared/private per patient
########################

gg_mut_freq_prop_case <-
  ggplot(mut_freq_prop_case, aes(x = case_barcode, y = mutation_percent, fill=mutation_type)) +
  geom_bar(stat="identity") +
  labs(y = "% mutation") +
  scale_fill_manual(values=c("R" = "#2FB3CA", "P" ="#CA2F66", "S"="#CA932F")) +
  plot_grid

testPlot(gg_mut_freq_prop_case)

########################
## Plot proportions shared/private per gene
########################

gg_snv_prop_gene <-
  snvgdata %>%
  filter(idh_codel_subtype == "all") %>%
  group_by(idh_codel_subtype) %>% 
  gather(key = "type", value = "value", prop_shared, prop_private_a, prop_private_b) %>%
  mutate(type = factor(type,
                       levels = c("prop_private_b", "prop_private_a", "prop_shared"),
                       labels = c("R", "P", "S"))) %>%
  ggplot(aes(x = gene_symbol, y = value, fill = type)) +
  geom_bar(stat="identity") +
  labs(y = "% mutation", x = "Gene Symbol") +
  scale_fill_manual(values=c("R" = "#2FB3CA", "P" ="#CA2F66", "S"="#CA932F")) +
  coord_flip()

testPlot(gg_snv_prop_gene)

pdf(file = "/Users/johnsk/Documents/snv_prop_genes.pdf", width = 12, height = 8)
testPlot(gg_snv_prop_gene)
dev.off()

pdf(file = "~/The Jackson Laboratory/GLASS - Documents/Figure 1/snv_prop_gene_all.pdf", width = 12, height = 8)
testPlot(gg_snv_prop_gene)
dev.off()

gg_cnv_prop_gene <-
  cnvgdata %>%
  filter(idh_codel_subtype == "all") %>%
  group_by(idh_codel_subtype) %>% 
  gather(key = "type", value = "value", prop_shared, prop_private_a, prop_private_b) %>%
  mutate(type = factor(type,
                       levels = c("prop_private_b", "prop_private_a", "prop_shared"),
                       labels = c("R", "P", "S"))) %>%
  ggplot(aes(x = gene_symbol, y = value, fill = type)) +
  geom_bar(stat="identity") +
  labs(y = "% cnv", x = "Gene Symbol") +
  scale_fill_manual(values=c("R" = "#2FB3CA", "P" ="#CA2F66", "S"="#CA932F")) +
  coord_flip()

testPlot(gg_cnv_prop_gene)

pdf(file = "/Users/johnsk/Documents/cnv_prop_genes.pdf", width = 12, height = 8)
testPlot(gg_cnv_prop_gene)
dev.off()


pdf(file = "~/The Jackson Laboratory/GLASS - Documents/Figure 1/cnv_prop_gene_all.pdf", width = 12, height = 8)
testPlot(gg_cnv_prop_gene)
dev.off()

########################
## Plot mutations
########################

gg_mut_heatmap <-
  ggplot(mut_heatmap, aes(x=case_barcode, y=gene_symbol)) +
  geom_tile(aes(fill = covered)) +
  geom_point(aes(shape = variant_call, color = variant_classification), size=1) +
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
  labs(y = "Gene Symbol", color = "Variant Classification", fill = "Coverage", shape = "SNV evolution")

testPlot(gg_mut_heatmap)

########################
## Plot CNV
########################

gg_cnv_heatmap <-
  ggplot(cnv_heatmap, aes(x=case_barcode, y=gene_symbol)) + 
  geom_tile(fill = "#f7f7f7") +
  geom_point(aes(color = cnv_state, shape = cnv_change), size=1) + #
  scale_color_manual(values=c("HLDEL" = "#0571b0", #0 = "#f7f7f7",
                              "HLAMP" = "#ca0020")) +
  labs(y = "Gene Symbol", x= "", color = "CNV type", shape = "CNV evolution") 

testPlot(gg_cnv_heatmap)   

########################
## Layout plot elements
########################

gleg1 <- gg_rbind(gg_legend(gg_clinical),
                  gg_legend(gg_time_data),
                  gg_legend(gg_evolution),
                  heights = rep(1,3),
                  ncol = 1)

gleg2 <- gg_rbind(gg_legend(gg_driver_count),
                  gg_legend(gg_driver_context),
                  gg_legend(gg_mut_freq_prop_case),
                  heights = rep(1,3),
                  ncol = 1)

gleg3 <- gg_rbind(gg_legend(gg_mut_heatmap),
                  gg_legend(gg_cnv_heatmap),
                  gg_legend(gg_mut_freq_case),
                  heights = rep(1,3),
                  ncol = 1)
pdf(file = "~/The Jackson Laboratory/GLASS - Documents/Figures/Figure1/legends-new.pdf", height = 12, width = 4)
plot(gleg1)
plot(gleg2)
plot(gleg3)
dev.off()

g0 <- ggplotGrob(gg_mut_freq_case + plot_grid + plot_theme + null_legend + null_x + null_facet + top_margin) %>% gtable_frame()
g1 <- ggplotGrob(gg_time_data + plot_grid + plot_theme + null_legend + null_x + null_facet + top_margin) %>% gtable_frame()
g2 <- ggplotGrob(gg_clinical + plot_grid + plot_theme + null_legend + null_x + null_facet + top_margin) %>% gtable_frame()
g3 <- ggplotGrob(gg_c710 + plot_grid + plot_theme + null_legend + null_x + null_y + null_facet + top_margin) %>% gtable_frame()
g4 <- ggplotGrob(gg_anpl_data + plot_grid + plot_theme + null_legend + null_x + null_facet + bottom_margin) %>% gtable_frame()
g5 <- ggplotGrob(gg_driver_count + plot_grid + plot_theme + null_legend + null_x + null_facet + bottom_margin) %>% gtable_frame()
g6 <- ggplotGrob(gg_driver_context + plot_grid + plot_theme + null_legend + null_x + null_facet ) %>% gtable_frame()
g7 <- ggplotGrob(gg_mut_freq_prop_case + plot_grid + plot_theme + null_legend + null_x + null_facet + bottom_margin) %>% gtable_frame()
g8 <- ggplotGrob(gg_mut_heatmap + plot_grid + plot_theme + null_legend + null_x + null_facet ) %>% gtable_frame()
g9 <- ggplotGrob(gg_cnv_heatmap + plot_grid + plot_theme + null_legend + bottom_x + null_facet + bottom_margin) %>% gtable_frame()
g10 <- ggplotGrob(gg_evolution + plot_grid + plot_theme + null_legend + null_x + null_y + null_facet) %>% gtable_frame()

g <- gg_rbind(g0, g1, g2, g3, g4, g5, g6, g7, g8, g9, heights = c(5, 4, 5, 1, 4, 4, 7, 4, 12, 12), ncol = 1)
plot(g)

pdf(file = "~/The Jackson Laboratory/GLASS - Documents/Figure 1/fig_cnv_drivers.pdf", height = 16, width = 12)
plot(g)
dev.off()

g_2 <- gg_rbind(g0, g7, g2 , g1 , heights = c(5, 5, 5, 4), ncol = 1)
plot(g_2)

pdf(file = "~/The Jackson Laboratory/GLASS - Documents/Figures/Figure1/f1a-fpb.pdf", height = 12, width = 16, bg = "transparent", useDingbats = FALSE)
plot(g_2) 
dev.off()

## plot with sidebars
gb <- ggplotGrob(gg_blank + plot_theme + null_legend + null_x + null_y + null_facet + top_margin) %>% gtable_frame()
g8b <- ggplotGrob(gg_snv_prop_gene + plot_grid + plot_theme + null_legend + null_x + null_y + null_facet) %>% gtable_frame()
g9b <- ggplotGrob(gg_cnv_prop_gene + plot_grid + plot_theme + null_legend + bottom_x + null_y + null_facet + bottom_margin) %>% gtable_frame()

g0sb <- gg_cbind(g0, gb, widths = c(0.1,8))
g1sb <- gg_cbind(g1, gb, widths = c(0.1,8))
g2sb <- gg_cbind(g2, gb, widths = c(0.1,8))
g3sb <- gg_cbind(g3, gb, widths = c(0.1,8))
g4sb <- gg_cbind(g4, gb, widths = c(0.1,8))
g5sb <- gg_cbind(g5, gb, widths = c(0.1,8))
g6sb <- gg_cbind(g6, gb, widths = c(0.1,8))
g7sb <- gg_cbind(g7, gb, widths = c(0.1,8))
g8sb <- gg_cbind(g8, g8b, widths = c(0.1,8))
g9sb <- gg_cbind(g9, g9b, widths = c(0.1,8))
g10sb <- gg_cbind(g10, gb, widths = c(0.1,8))

gsb <- gg_rbind(g0sb, g1sb, g2sb, g3sb, g4sb, g5sb, g6sb, g7sb, g8sb, g9sb, heights = c(5, 4, 5, 1, 4, 4, 7, 4, 12, 12), ncol = 2)
plot(gsb)

pdf(file = "~/The Jackson Laboratory/GLASS - Documents/Figure 1/fig_cnv_drivers_w_sidebar.pdf", height = 16, width = 12)
plot(gsb)
dev.off()

gsb <- gg_rbind(g2sb, g6sb, g10sb, g8sb, g4sb, g9sb, heights = c(5, 3, 1, 12, 3, 12), ncol = 2)
plot(gsb)

pdf(file = "~/The Jackson Laboratory/GLASS - Documents/Figures/Figure3/f3a-fpb.pdf", height = 12, width = 16, bg = "transparent", useDingbats = FALSE)
plot(gsb)
dev.off()


