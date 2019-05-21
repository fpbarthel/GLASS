## Data cleanup
library(tidyverse)

## Database
library(odbc)
library(DBI)

## Plotting
library(ggthemes)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
#library(egg)
library(RColorBrewer)
library(wesanderson)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

q = "SELECT ts.case_barcode, case_sex, s.idh_status, s.codel_status, s.idh_codel_subtype, v.gene_symbol, ts.sample_type, ts.sample_barcode, 
v.chrom, v.start, v.end, v.alt, v.variant_classification, v.variant_type, 
gt.aliquot_barcode, v.hgvs_p, gt.ref_count, gt.alt_count, gt.read_depth, sift, polyphen, gt.called, mf.coverage_adj_mut_freq
FROM analysis.snvs v
FULL JOIN analysis.snv_genotypes gt ON v.chrom = gt.chrom AND v.start = gt.start AND v.end = gt.end AND v.alt = gt.alt
INNER JOIN biospecimen.aliquots ta ON ta.aliquot_barcode = gt.aliquot_barcode
INNER JOIN biospecimen.samples ts ON ts.sample_barcode = ta.sample_barcode
LEFT JOIN analysis.mutation_freq mf ON mf.aliquot_barcode = gt.aliquot_barcode
LEFT JOIN clinical.surgeries s ON s.sample_barcode = ts.sample_barcode
INNER JOIN clinical.cases ca ON ts.case_barcode = ca.case_barcode
WHERE ts.sample_type IN ('TP','R1','R2','R3') AND ((v.variant_classification IN ('In_Frame_Ins','In_Frame_Del','Frame_Shift_Del','Frame_Shift_Ins','Missense_Mutation','Nonstop_Mutation','Nonsense_Mutation')
AND v.gene_symbol IN ('TP53','ATRX','RB1','EGFR','PTEN','NF1','CIC','FUBP1','PDGFRA','PIK3CA','PIK3R1','NOTCH1')) 
OR (v.variant_classification = '5''Flank' AND v.gene_symbol = 'TERT' AND (v.start = 1295250 OR v.start = 1295228)))"

qres <- dbGetQuery(con, q)

q2 <- "SELECT case_barcode, sample_barcode, surgical_interval_mo, histology, grade, idh_status, codel_status, who_classification, surgery_location FROM clinical.surgeries"

q2res <- dbGetQuery(con, q2)

#q3 <- "SELECT * FROM analysis.mutation_freq"
#
#q3res <- dbGetQuery(con, q3)

df = qres %>% 
  mutate(var = sprintf("%s:%s-%s_%s", chrom, start, end, alt),
         called = called == "1",
         severity_score = case_when(variant_classification == "Nonsense_Mutation" ~ 0,
                                    variant_classification %in% c("Frame_Shift_Del","Frame_Shift_Del") ~ 1,
                                    variant_classification %in% c("In_Frame_Del","In_Frame_Ins") ~ 2,
                                    variant_classification == "Missense_Mutation" ~ 3,
                                    variant_classification == "5'Flank" ~ 4)) %>%
  group_by(sample_barcode, var) %>%
  mutate(optimal_variant = order(read_depth, decreasing = T)) %>%
  ungroup() %>%
  filter(optimal_variant == 1) %>% ## (2) For each sample/variant combination, select the variant with the highest read depth, eg. when a variant was profiled across multple sectors or with both WGS and WES
  group_by(case_barcode, gene_symbol) %>%
  mutate(any_called = any(called, na.rm=T),
         num_samples = n_distinct(sample_barcode)) %>%
  ungroup() %>%
  filter(num_samples > 1) %>% ## (4) filter out singletons (unpaired samples)
  group_by(gene_symbol) %>%
  mutate(num_patient = n_distinct(case_barcode),
         gene_symbol_label = sprintf("%s (n=%s)", gene_symbol, num_patient)) %>%
  ungroup() %>%
  arrange(desc(called), severity_score) %>%
  group_by(sample_barcode, gene_symbol) %>%
  mutate(avg_cov = mean(read_depth)) %>%
  filter(row_number()==1) %>%
  ungroup() %>%
  left_join(q2res)

# theme_base <- theme(plot.title = element_text(size = 10), #plot.margin = unit(c(0,0.5,0,0.5),"lines"),
#                     panel.spacing = unit(0.2, "lines"),
#                     axis.title.x = element_blank(),#element_text(size = 12, angle = 0),
#                     axis.title.y = element_blank(),#element_text(size = 12, angle = 90),
#                     panel.background = element_rect(fill = "white",colour=NA),
#                     panel.border=element_blank(),
#                     plot.background = element_rect(fill = NA),
#                     axis.text.x = element_text(size=12, angle=0,colour = 'black', hjust = 0.5, vjust = 1),
#                     axis.text.y = element_text(size=18, angle=0,colour = 'black',  hjust = 1, vjust = 0.5, face = 'italic'),
#                     legend.text = element_text(size=16),
#                     legend.key.size = unit(1.5, "line"),
#                     legend.key = element_rect(fill="white",colour=NA),
#                     legend.background = element_blank(),
#                     legend.title = element_text(size=18,face='bold'),
#                     axis.ticks.x=element_blank(),
#                     axis.ticks.y=element_blank(),
#                     axis.ticks.length = unit(0,"null"),
#                     axis.ticks=element_line(colour="white"),
#                     strip.background = element_rect(fill="white",colour=NA),
#                     strip.text.x = element_text(size=14, colour="black",face="bold",angle = 90),
#                     strip.text.y = element_text(size=18, colour="black",face="bold.italic",angle = 0))

g_theme = theme_minimal(base_size = 12, base_family = "sans") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
null_legend = theme(legend.position = 'none')
null_x = theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank()) 
bottom_x = theme(axis.text.x=element_blank()) 
null_facet = theme(strip.background = element_blank(),
                   strip.text.x = element_blank())
top_margin = theme(plot.margin= unit(c(1, 1, 0.1, 1), "lines")) ## Top, Right, Bottom, Left
middle_margin = theme(plot.margin= unit(c(0, 1, 0.1, 1), "lines"))
bottom_margin = theme(plot.margin= unit(c(0, 1, 1, 1), "lines"))
plot_grid = facet_grid(. ~ case_barcode, scales = "free_x", space = "free")

codel_subtype = "IDHwt_noncodel"

for (codel_subtype in unique(df$idh_codel_subtype)) {

  ggdat <- df %>%
    filter(idh_codel_subtype == codel_subtype) %>%
    complete(nesting(case_barcode, sample_barcode, sample_type), gene_symbol) %>%
    select(case_barcode, sample_barcode, idh_codel_subtype, sample_type, gene_symbol, variant_classification, called, avg_cov, mf = coverage_adj_mut_freq) %>%
    left_join(q2res) %>%
    mutate(variant = ifelse(called, gsub("_","-",variant_classification), ifelse(avg_cov < 15, NA, "Wild-Type")),
           variant_classification = factor(variant_classification),
           gene_symbol = factor(gene_symbol),
           sample_type = factor(sample_type, levels = c("TP","R1","R2","R3"), labels = c("P","2","3","4")),
           sample_barcode = fct_reorder2(sample_barcode, as.numeric(variant_classification), as.numeric(gene_symbol), .desc = FALSE),
           case_barcode = fct_reorder2(case_barcode, as.numeric(variant_classification), as.numeric(gene_symbol), .desc = FALSE))
  
  p1 <- ggplot() + 
    geom_tile(data = ggdat, aes(x = sample_type, y = "WHO Classification", fill = who_classification), color = "black") +
    labs(y="Gene", fill = "WHO Classification")
  
  #p1 <- p1 + facet_grid(. ~ case_barcode, scales = "free_x", space = "free_x", drop = F) + theme_base
  
  p2 <- ggplot() + 
    geom_tile(data = ggdat, aes(x = sample_type, y = gene_symbol, fill = variant), color = "black") +
    labs(y="Gene", fill = "Variant Classification") +
    scale_fill_manual(values=c("5'Flank" = brewer.pal(7, "Paired")[7],
                               "Frame-Shift-Del" = brewer.pal(7, "Paired")[1],
                               "Frame-Shift-Ins" = brewer.pal(7, "Paired")[2],
                               "In-Frame-Del" = brewer.pal(7, "Paired")[3],
                               "In-Frame-Ins" = brewer.pal(7, "Paired")[4],
                               "Missense-Mutation" = brewer.pal(7, "Paired")[5],
                               "Nonsense-Mutation" = brewer.pal(7, "Paired")[6],
                               "Wild-Type" = "white"), na.value = "gray75")
  
  g1 = ggplotGrob(p1 + plot_grid + g_theme + null_legend + null_x + top_margin)                  %>% gtable_frame()
  g2 = ggplotGrob(p2 + plot_grid + g_theme + null_legend + null_x + null_facet + middle_margin + theme(axis.title = element_text(size = 10), axis.text = element_text(size=10)))  %>% gtable_frame()
  #g2 = ggplotGrob(p2 + plot_grid + g_theme + null_legend + null_x + null_facet + middle_margin + theme(axis.title = element_text(size = 10), axis.text = element_text(size=10)))  %>% gtable_frame()
  #g3 = ggplotGrob(p3 + plot_grid + g_theme + null_legend + null_x + null_facet + middle_margin)  %>% gtable_frame()
  #g4 = ggplotGrob(p4 + plot_grid + g_theme + null_legend + null_x + null_facet + middle_margin)  %>% gtable_frame()
  #g5 = ggplotGrob(p5 + plot_grid + g_theme + null_legend + null_x + null_facet + middle_margin)  %>% gtable_frame()
  #g6 = ggplotGrob(p6 + plot_grid + g_theme + null_legend + null_x + null_facet + middle_margin)  %>% gtable_frame()
  #g7 = ggplotGrob(p7 + plot_grid + g_theme + null_legend + null_x + null_facet + middle_margin)  %>% gtable_frame()
  #g8 = ggplotGrob(p8 + gene_grid + g_theme + null_legend + bottom_x + null_facet + bottom_margin + theme(strip.text.y = element_text(size = 12)))      %>% gtable_frame()
  
  g = gtable_rbind(g1, g2)#, g2, g3, g4, g5, g6, g7, g8)
  #gleg = gtable_rbind(gleg0, gleg2, gleg3, gleg4, gleg5, gleg6, gleg7, gleg8)
  
  ## Adjust relative height of panels
  panels = g$layout$t[grep("panel", g$layout$name)]
  g$heights[panels] <- unit(c(0.1,2), "null")
  
  #pdf(file = "figures/HM_v7_Verhaak2010.pdf", width = 8.5, height = 11)
  #grid.newpage()
  #grid.draw(g)
  #dev.off()
  
  
  #p2 <- p2 + facet_grid(. ~ case_barcode, scales = "free_x", space = "free_x", drop = F) + theme_base +
  #  theme(strip.text.x = element_blank())
  
  #p <- grid.arrange(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "first"),
  #                  nrow=2, ncol=1, top=NULL,
   #                 widths = c(10),
  #                  heights = c(50, 350),
  #                  layout_matrix = rbind(c(1),
  #                                        c(2)))
    
  #grid.newpage()
  #grid.draw(p)
  
  ggsave(g,
         file=sprintf("%s-test.pdf",codel_subtype),
         width=40,
         height=10)
  
}

# ## plot all subtypes
# 
# ggdat = df %>%
#   complete(nesting(case_barcode, sample_barcode, sample_type), gene_symbol) %>%
#   select(case_barcode, sample_barcode, idh_codel_subtype, sample_type, gene_symbol, variant_classification, called, avg_cov) %>%
#   left_join(q2res) %>%
#   mutate(variant = ifelse(called, gsub("_","-",variant_classification), ifelse(avg_cov < 15, NA, "Wild-Type")),
#          variant_classification = factor(variant_classification),
#          gene_symbol = factor(gene_symbol),
#          sample_type = factor(sample_type, levels = c("TP","R1","R2","R3"), labels = c("P","2","3","4")),
#          gene_symbol = fct_reorder(gene_symbol, as.numeric(variant_classification), .desc = FALSE),
#          sample_barcode = fct_reorder2(sample_barcode, as.numeric(variant_classification), as.numeric(gene_symbol), .desc = FALSE),
#          case_barcode = fct_reorder2(case_barcode, as.numeric(variant_classification), as.numeric(gene_symbol), .desc = FALSE))
# 
# p2 = ggplot() + 
#   geom_tile(data = ggdat, aes(x = gene_symbol, y = sample_type , fill = variant), color = "black") +
#   labs(y="Gene", fill = "Variant Classification") +
#   scale_fill_manual(values=c("5'Flank" = brewer.pal(7, "Paired")[7],
#                              "Frame-Shift-Del" = brewer.pal(7, "Paired")[1],
#                              "Frame-Shift-Ins" = brewer.pal(7, "Paired")[2],
#                              "In-Frame-Del" = brewer.pal(7, "Paired")[3],
#                              "In-Frame-Ins" = brewer.pal(7, "Paired")[4],
#                              "Missense-Mutation" = brewer.pal(7, "Paired")[5],
#                              "Nonsense-Mutation" = brewer.pal(7, "Paired")[6],
#                              "Wild-Type" = "white"), na.value = "gray75")
# 
# p2 = p2 + facet_grid(case_barcode ~ ., scales = "free_y", space = "free_y", drop = F) + 
#   scale_x_discrete(position = "top") +
#   theme(plot.title = element_text(size = 1), #plot.margin = unit(c(0,0.5,0,0.5),"lines"),
#         panel.spacing = unit(0.05, "lines"),
#         axis.title.x = element_blank(),#element_text(size = 12, angle = 0),
#         axis.title.y = element_blank(),#element_text(size = 12, angle = 90),
#         panel.background = element_rect(fill = "white",colour=NA),
#         panel.border=element_blank(),
#         plot.background = element_rect(fill = NA),
#         axis.text.x = element_text(size=10, angle=0,colour = 'black', hjust = 0.5, vjust = 0.5, face = 'italic'),
#         axis.text.y = element_text(size=6, angle=0, colour = 'black', hjust = 0.5, vjust = 0.5),
#         legend.text = element_text(size=1),
#         legend.key.size = unit(1, "line"),
#         legend.key = element_rect(fill="white",colour=NA),
#         legend.background = element_blank(),
#         legend.title = element_text(size=1,face='bold'),
#         axis.ticks.x=element_blank(),
#         axis.ticks.y=element_blank(),
#         axis.ticks.length = unit(0,"null"),
#         axis.ticks=element_line(colour="white"),
#         strip.background = element_rect(fill="white",colour=NA),
#         strip.text.x = element_text(size=8, colour="black",face="bold",angle = 90),
#         strip.text.y = element_text(size=8, colour="black",face="bold",angle = 0))
#   
# ggsave(p2,
#        file="all.pdf",
#        width=10,
#        height=30)
