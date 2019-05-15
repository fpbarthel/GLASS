library(DBI)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(egg)

con <- DBI::dbConnect(odbc::odbc(), "GLASSv2")

########################
## Load heatmap data
########################

hmapdata <- dbGetQuery(con, read_file("sql/heatmap/heatmap_snv.sql"))
snvgdata <- dbGetQuery(con, read_file("sql/heatmap/heatmap_snv_by_gene.sql"))
cnv_data <- dbGetQuery(con, read_file("sql/heatmap/heatmap_cnv.sql"))
cnvgdata <- dbGetQuery(con, read_file("sql/heatmap/heatmap_cnv_by_gene.sql"))

arm_data <- dbGetQuery(con, read_file("sql/heatmap/heatmap_arm.sql"))
armadata <- dbGetQuery(con, read_file("sql/heatmap/heatmap_arm_by_arm.sql"))

mutfdata <- dbGetQuery(con, read_file("sql/heatmap/heatmap_mf.sql"))
clindata <- dbGetQuery(con, read_file("sql/heatmap/heatmap_clinical.sql"))
timedata <- dbGetQuery(con, read_file("sql/heatmap/heatmap_time.sql"))
anpldata <- dbGetQuery(con, read_file("sql/heatmap/heatmap_aneuploidy.sql"))
c710data <- dbGetQuery(con, read_file("sql/heatmap/heatmap_c710.sql"))
puridata <- dbGetQuery(con, read_file("sql/heatmap/heatmap_purity.sql"))
drivdata <- dbGetQuery(con, read_file("sql/heatmap/heatmap_drivers.sql"))

neutrdata <- dbGetQuery(con, read_file("sql/heatmap/heatmap_evolution.sql"))
pycldata  <- dbGetQuery(con, read_file("sql/heatmap/heatmap_pyclone_clusters.sql"))

## Check counts (cases)

n_distinct(hmapdata$case_barcode)
n_distinct(cnv_data$case_barcode)
n_distinct(arm_data$case_barcode)
n_distinct(mutfdata$case_barcode)
n_distinct(clindata$case_barcode)
n_distinct(timedata$case_barcode)
n_distinct(anpldata$case_barcode)
n_distinct(c710data$case_barcode)
n_distinct(puridata$case_barcode)
n_distinct(drivdata$case_barcode)
n_distinct(neutrdata$case_barcode)
n_distinct(pycldata$case_barcode)

## Check counts (genes)

n_distinct(hmapdata$gene_symbol)
n_distinct(snvgdata$gene_symbol)
n_distinct(cnv_data$gene_symbol)
n_distinct(cnvgdata$gene_symbol)

## Check counts (arms)
n_distinct(armadata$arm)

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
  left_join(neutrdata) %>%
  #arrange(desc(cnv_driver_count + snv_driver_count))
  #arrange(evolution_a,evolution_b)
  arrange(desc(mf_b)) # snv_driver_count, cnv_driver_count)#, aneuploidy_b - aneuploidy_a, 

########################
# Data re-ordering: order plot elements to match one another
########################

case_order <- unique(sort_df$case_barcode)
snv_gene_order <- rev(unique(snvgdata$gene_symbol))
cnv_gene_order <- rev(unique(cnvgdata$gene_symbol))
arm_order <- rev(c("1p","19q","7p","7q","10p","10q"))

c710_data <- c710data %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
time_data <- time_data %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
clin_data <- clindata %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
mut_freq_case <- mut_freq_case %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
anpl_data <- anpldata %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
puri_data <- puridata %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
mut_freq_prop_case <- mut_freq_prop_case %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
arm_heatmap <- arm_data %>% mutate(case_barcode = factor(case_barcode, levels = case_order),
                                   segment_qual = ifelse(!is.na(cnv_state), "Contiguous", "Disrupted"),
                                   cnv_state = ifelse(cnv_state == "neut", NA, cnv_state),
                                   arm = factor(arm, levels = arm_order))
armadata <- armadata %>% mutate(arm = factor(arm, levels = arm_order))
cnv_heatmap <- cnv_data %>% mutate(case_barcode = factor(case_barcode, levels = case_order),
                                   gene_symbol = factor(gene_symbol, levels = cnv_gene_order))
mut_heatmap <- hmapdata %>% mutate(case_barcode = factor(case_barcode, levels = case_order),
                                   gene_symbol = factor(gene_symbol, levels = snv_gene_order),
                                   covered = factor(covered, levels = c("Coverage < 5x", "Coverage >= 5x", "Coverage >= 15x", "Coverage >= 30x")))
driv_hm <- drivdata %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
snvgdata <- snvgdata %>% mutate(gene_symbol = factor(gene_symbol, levels = snv_gene_order))
cnvgdata <- cnvgdata %>% mutate(gene_symbol = factor(gene_symbol, levels = cnv_gene_order))
neutrdata <- neutrdata %>% mutate(case_barcode = factor(case_barcode, levels = case_order))
pycldata <- pycldata %>% mutate(case_barcode = factor(case_barcode, levels = case_order))

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
## Plot PyCle clusters
########################

gg_pyclone_clusters <-
  ggplot(pycldata, aes(x=case_barcode, fill = cut(pycldata$cluster_size,c(0,2,10,30,Inf)))) + #fill = cut(min_ccf,breaks = c(-Inf,0.25,0.50,0.75,Inf)))) + #, fill = cut(pycldata$cluster_size,c(0,2,10,30,Inf)))) +
  geom_bar(position = "stack") +
  labs(y = "No. of PyClone Clusters", fill = "Cluster Size", x = "Patient") +
  scale_fill_manual(values = c("#fef0d9", "#fdcc8a", "#fc8d59", "#d7301f")) +
  scale_y_continuous(breaks = c(0,4,8,12))

testPlot(gg_pyclone_clusters)

########################
## Plot PyClone dominant cluster 
########################

tmp <- pycldata %>% 
  group_by(case_barcode) %>% 
  summarize(ccf=max(min_ccf)) %>% 
  ungroup() %>%
  filter(complete.cases(ccf), !duplicated(ccf)) %>%
  arrange(desc(ccf)) %>%
  mutate(x=row_number()/n())

ggtmp <- ggplot(tmp, aes(x, ccf)) + 
  geom_hline(yintercept = max(tmp$x[tmp$ccf>0.25]), linetype = 2, color = "blue") +
  geom_vline(xintercept = 0.25, linetype = 2, color = "blue") +
  geom_hline(yintercept = max(tmp$x[tmp$ccf>0.50]), linetype = 2, color = "green") +
  geom_vline(xintercept = 0.50, linetype = 2, color = "green") +
  geom_text(x=0.25-0.025, y = max(tmp$x[tmp$ccf>0.25])-0.025, label = sprintf("%s%%", round(max(tmp$x[tmp$ccf>0.25])*100,1)), color = "blue", hjust = 1, vjust = 1) +
  geom_text(x=0.50-0.025, y = max(tmp$x[tmp$ccf>0.50])-0.025, label = sprintf("%s%%", round(max(tmp$x[tmp$ccf>0.50])*100,1)), color = "green", hjust = 1, vjust = 1) +
  geom_point() + 
  geom_line() + 
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) + 
  labs(x= "Mininum Dominant Persistant Cluster Cancer Cell Fraction", y = "Proportion of Patients") + 
  theme_bw(base_size = 12)

ggtmp

########################
## Plot clinical
########################

gg_clinical <-
  clin_data %>% 
  gather(key = "type", value = "value", grade_change, received_rt,received_alk, is_hypermutator) %>%
  mutate(type = factor(type,
                       levels = c("grade_change", "received_rt","received_alk",  "is_hypermutator"),
                       labels = c("Grade increase", "Received RT","Received Alkyl.", "Hypermutator"))) %>%
  ggplot(aes(x=case_barcode)) +
  geom_tile(aes(fill = factor(value), y = type)) +
  scale_fill_manual(values=c("white", "#377eb8"), na.value = "gray90") +
  labs(y="", fill = "Event")

testPlot(gg_clinical)

########################
## Plot co-variates
########################

gg_evolution <-
  neutrdata %>%
  gather(key = "type", value = "value", evolution_a, evolution_b) %>%
  mutate(type = factor(type,
                       levels = c("evolution_a", "evolution_b"),
                       labels = c("Primary", "Recurrence"))) %>%
  ggplot(aes(x=case_barcode)) +
  geom_tile(aes(fill = factor(value), y = type)) +
  scale_fill_manual(values = c("N" = "#A3BD67", "S" = "#BD8167")) +
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
  coord_cartesian(ylim=c(-0.5,0.5)) +
  labs(y = "Aneuploidy \ndifference")

testPlot(gg_anpl_data)

########################
## Plot driver data
########################

gg_driver_count <- 
  driv_hm %>% 
  gather(key = "type", value = "value", snv_driver_count_shared, snv_driver_count_private_a, snv_driver_count_private_b, cnv_driver_count_shared, cnv_driver_count_private_a, cnv_driver_count_private_b) %>%
  mutate(type = factor(type,
                       levels = c("snv_driver_count_shared", "cnv_driver_count_shared", "snv_driver_count_private_a", "cnv_driver_count_private_a", "snv_driver_count_private_b", "cnv_driver_count_private_b"),
                       labels = c("S - SNV", "S - CNV", "P - SNV", "P - CNV", "R - SNV", "R - CNV"))) %>%
  ggplot(aes(x=case_barcode)) +
  geom_bar(aes(y=as.integer(value), fill = type), stat = "identity") + 
  labs(y = "# of drivers") +
  scale_fill_brewer(palette = "Paired", type = "qual") +
  scale_y_continuous(breaks = c(0,4,8,12)) +
  coord_cartesian(ylim=c(0,12))

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
## Plot proportions shared/private per gene (SNVs)
########################

tmp1 <- snvgdata %>%
  group_by(idh_codel_subtype) %>% 
  gather(key = "type", value = "value", prop_shared, prop_private_a, prop_private_b) %>%
  mutate(type = factor(type,
                       levels = c("prop_private_b", "prop_private_a", "prop_shared"),
                       labels = c("R", "P", "S"))) %>%
  ungroup()

gg_snv_prop_gene_all <- ggplot(tmp1, aes(x = gene_symbol, y = value, fill = type)) +
  geom_bar(stat="identity") + 
  geom_text(aes(label = sprintf("n=%s", total)), y = 0.90, stat = "identity", check_overlap = TRUE) +
  labs(y = "% mutation", x = "Gene Symbol") +
  scale_fill_manual(values=c("R" = "#2FB3CA", "P" ="#CA2F66", "S"="#CA932F")) +
  coord_flip()

testPlot(gg_snv_prop_gene_all)

gg_snv_prop_gene <- tmp1 %>% filter(idh_codel_subtype == "all") %>% ggplot(aes(x = gene_symbol, y = value, fill = type)) +
  geom_bar(stat="identity") +
  labs(y = "% mutation", x = "Gene Symbol") +
  scale_fill_manual(values=c("R" = "#2FB3CA", "P" ="#CA2F66", "S"="#CA932F")) +
  coord_flip()

testPlot(gg_snv_prop_gene)

########################
## Plot proportions shared/private per gene (CNV)
########################

tmp2 <- cnvgdata %>%
  group_by(idh_codel_subtype) %>% 
  gather(key = "type", value = "value", prop_shared, prop_private_a, prop_private_b) %>%
  mutate(type = factor(type,
                       levels = c("prop_private_b", "prop_private_a", "prop_shared"),
                       labels = c("R", "P", "S"))) %>%
  ungroup()

gg_cnv_prop_gene_all <- ggplot(tmp2, aes(x = gene_symbol, y = value, fill = type)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = sprintf("n=%s", total)), y = 0.90, stat = "identity", check_overlap = TRUE) +
  labs(y = "% cnv", x = "Gene Symbol") +
  scale_fill_manual(values=c("R" = "#2FB3CA", "P" ="#CA2F66", "S"="#CA932F")) +
  coord_flip()

testPlot(gg_cnv_prop_gene_all)

gg_cnv_prop_gene <- tmp2 %>% filter(idh_codel_subtype == "all") %>% ggplot(aes(x = gene_symbol, y = value, fill = type)) +
  geom_bar(stat="identity") +
  labs(y = "% cnv", x = "Gene Symbol") +
  scale_fill_manual(values=c("R" = "#2FB3CA", "P" ="#CA2F66", "S"="#CA932F")) +
  coord_flip()

testPlot(gg_cnv_prop_gene)

########################
## Plot proportions shared/private per gene (arm-level)
########################

tmp3 <- armadata %>%
  group_by(idh_codel_subtype) %>% 
  gather(key = "type", value = "value", prop_shared, prop_private_a, prop_private_b) %>%
  mutate(type = factor(type,
                       levels = c("prop_private_b", "prop_private_a", "prop_shared"),
                       labels = c("R", "P", "S"))) %>%
  ungroup()

gg_arm_prop_arm_all <- ggplot(tmp3, aes(x = arm, y = value, fill = type)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = sprintf("n=%s", total)), y = 0.90, stat = "identity", check_overlap = TRUE) +
  labs(y = "% cnv", x = "Arm") +
  scale_fill_manual(values=c("R" = "#2FB3CA", "P" ="#CA2F66", "S"="#CA932F")) +
  coord_flip()

testPlot(gg_arm_prop_arm_all)

gg_arm_prop_arm <- tmp3 %>% filter(idh_codel_subtype == "all") %>% ggplot(aes(x = arm, y = value, fill = type)) +
  geom_bar(stat="identity") +
  labs(y = "% cnv", x = "Arm") +
  scale_fill_manual(values=c("R" = "#2FB3CA", "P" ="#CA2F66", "S"="#CA932F")) +
  coord_flip()

testPlot(gg_arm_prop_arm)

########################
## Plot mutations
########################

gg_mut_heatmap <-
  ggplot(mut_heatmap, aes(x=case_barcode, y=gene_symbol)) +
  geom_tile(aes(fill = covered)) +
  geom_point(aes(shape = variant_call, color = variant_classification), size=1) +
  scale_fill_manual(values=c("gray90", "gray95", "gray99", "white"), na.value = "gray90") +
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
## Plot arm cnv
########################

gg_arm_heatmap <-
  ggplot(arm_heatmap, aes(x=case_barcode, y=arm)) + 
  geom_tile(aes(fill=segment_qual)) +
  scale_fill_manual(values=c("white", "#f2f0f7"), na.value = "gray90") +
  geom_point(aes(color = cnv_state, shape = cnv_change), size=1) + #
  scale_color_manual(values=c("del" = "#0571b0", #0 = "#f7f7f7",
                              "amp" = "#ca0020")) +
  labs(y = "Chr. Arm", color = "CNV type", shape = "Arm evolution", fill = "Segment Quality")

testPlot(gg_arm_heatmap)   

########################
## Plot CNV
########################

gg_cnv_heatmap <-
  ggplot(cnv_heatmap, aes(x=case_barcode, y=gene_symbol)) + 
  geom_tile(aes(fill=gold_set)) +
  scale_fill_manual(values=c("white", "gray90"), na.value = "gray90") +
  geom_point(aes(color = cnv_state, shape = cnv_change), size=1) + #
  scale_color_manual(values=c("HLDEL" = "#0571b0", #0 = "#f7f7f7",
                              "HLAMP" = "#ca0020")) +
  labs(y = "Gene Symbol", color = "CNV type", shape = "CNV evolution", x = "Patient")

testPlot(gg_cnv_heatmap)   


########################
## Plot figure legends
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

gleg4 <- gg_rbind(gg_legend(gg_arm_heatmap),
                  gg_legend(gg_pyclone_clusters),
                  heights = rep(1,2),
                  ncol = 1)

pdf(file = "~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/legends-fig1-fig3.pdf", height = 12, width = 4)
plot(gleg1)
plot(gleg2)
plot(gleg3)
plot(gleg4)
dev.off()

########################
## Plot exploratory/overview heatmap
########################

g <- gg_rbind(gtable_frame(ggplotGrob(gg_mut_freq_case      + plot_grid + plot_theme + null_legend + null_x + null_facet)),
              gtable_frame(ggplotGrob(gg_pyclone_clusters   + plot_grid + plot_theme + null_legend + null_x + null_facet)),
              gtable_frame(ggplotGrob(gg_clinical           + plot_grid + plot_theme + null_legend + null_x + null_facet)),
              gtable_frame(ggplotGrob(gg_c710               + plot_grid + plot_theme + null_legend + null_x + null_y + null_facet)),
              gtable_frame(ggplotGrob(gg_anpl_data          + plot_grid + plot_theme + null_legend + null_x + null_facet)),
              gtable_frame(ggplotGrob(gg_driver_count       + plot_grid + plot_theme + null_legend + null_x + null_facet)),
              gtable_frame(ggplotGrob(gg_driver_context     + plot_grid + plot_theme + null_legend + null_x + null_facet)),
              gtable_frame(ggplotGrob(gg_mut_freq_prop_case + plot_grid + plot_theme + null_legend + null_x + null_facet)),
              gtable_frame(ggplotGrob(gg_evolution          + plot_grid + plot_theme + null_legend + null_x + null_facet)),
              gtable_frame(ggplotGrob(gg_arm_heatmap        + plot_grid + plot_theme + null_legend + null_x + null_facet)),
              gtable_frame(ggplotGrob(gg_mut_heatmap        + plot_grid + plot_theme + null_legend + null_x + null_facet)),
              gtable_frame(ggplotGrob(gg_cnv_heatmap        + plot_grid + plot_theme + null_legend + bottom_x + null_facet)),
              heights = c(4, 4, 4, 1, 4, 4, 3, 4, 2, 6, 12, 12),
              ncol = 1)
plot(g)

pdf(file = "~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/heatmap_full.pdf", height = 16, width = 12)
plot(g)
dev.off()

########################
## Plot figure 1 heatmap
########################

fig1a <- gg_rbind(gtable_frame(ggplotGrob(gg_mut_freq_case + plot_grid + plot_theme + null_legend + null_x + null_facet)),
                  gtable_frame(ggplotGrob(gg_mut_freq_prop_case + plot_grid + plot_theme + null_legend + null_x + null_facet)),
                  gtable_frame(ggplotGrob(gg_clinical + plot_grid + plot_theme + null_legend + null_x + null_facet)),
                  gtable_frame(ggplotGrob(gg_pyclone_clusters + plot_grid + plot_theme + null_legend + bottom_x + null_facet + bottom_margin)),
                  heights = c(5, 5, 4, 5), ncol = 1)
plot(fig1a)

pdf(file = "~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/Figure 1/f1a-fpb.pdf", height = 12, width = 16, bg = "transparent", useDingbats = FALSE)
plot(fig1a) 
dev.off()

########################
## Plot Figure 3 heatmap
########################

gb <- gtable_frame(ggplotGrob(gg_blank + plot_theme + null_legend + null_x + null_y + null_facet))  

fig3a <- gg_rbind(gg_cbind(gtable_frame(ggplotGrob(gg_clinical + plot_grid + plot_theme + null_legend + null_x + null_facet)), 
                           gb,
                           widths = c(0.1,8)),
                  gg_cbind(gtable_frame(ggplotGrob(gg_evolution     + plot_grid + plot_theme + null_legend + null_x + null_facet)),
                           gb,
                           widths = c(0.1,8)),
                  gg_cbind(gtable_frame(ggplotGrob(gg_anpl_data     + plot_grid + plot_theme + null_legend + null_x + null_facet)),
                           gb,
                           widths = c(0.1,8)),
                  gg_cbind(gtable_frame(ggplotGrob(gg_driver_count  + plot_grid + plot_theme + null_legend + null_x + null_facet)),
                           gb,
                           widths = c(0.1,8)),
                  gg_cbind(gtable_frame(ggplotGrob(gg_mut_heatmap   + plot_grid + plot_theme + null_legend + null_x + null_facet)),
                           gtable_frame(ggplotGrob(gg_snv_prop_gene + plot_grid + plot_theme + null_legend + null_x + null_y + null_facet)),
                           widths = c(0.1,8)),
                  gg_cbind(gtable_frame(ggplotGrob(gg_arm_heatmap   + plot_grid + plot_theme + null_legend + null_x + null_facet)),
                           gtable_frame(ggplotGrob(gg_arm_prop_arm  + plot_grid + plot_theme + null_legend + null_x + null_y + null_facet)),
                           widths = c(0.1,8)),
                  gg_cbind(gtable_frame(ggplotGrob(gg_cnv_heatmap   + plot_grid + plot_theme + null_legend + bottom_x + null_facet)),
                           gtable_frame(ggplotGrob(gg_cnv_prop_gene + plot_grid + plot_theme + null_legend + bottom_x + null_y + null_facet)),
                           widths = c(0.1,8)),
                  heights = c(4, 2, 3, 3, 12, 6, 12),
                  ncol = 2)
plot(fig3a)

pdf(file = "~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/Figure 3/f3a-fpb.pdf", height = 16, width = 16, bg = "transparent", useDingbats = FALSE)
plot(fig3a)
dev.off()

########################
## Plot PyClone minimum dominant cluster CCF
########################

pdf(file = "~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/min_dominant_persistant_cluster_ccf.pdf", width = 5, height = 5, useDingbats = FALSE)
plot(ggtmp)
dev.off()

########################
## Plot supplementary figure proportions (SNV)
########################

pdf(file = "~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/snv_prop_gene_all.pdf", width = 12, height = 8, useDingbats = FALSE)
testPlot(gg_snv_prop_gene_all)
dev.off()

########################
## Plot supplementary figure proportions (CNV)
########################

pdf(file = "~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/cnv_prop_genes_all.pdf", width = 12, height = 8, useDingbats = FALSE)
testPlot(gg_cnv_prop_gene_all)
dev.off()

########################
## Plot supplementary figure proportions (arm)
########################

pdf(file = "~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/arm_prop_arm_all.pdf", width = 12, height = 8, useDingbats = FALSE)
testPlot(gg_arm_prop_arm_all)
dev.off()
