
library(DBI)
library(tidyverse)
library(ggplot2)
library(ggforce)
library(egg)
library(directlabels)

con <- DBI::dbConnect(odbc::odbc(), "GLASSv2")
res <- dbGetQuery(con, read_file("sql/pyclone/pyclone_cluster_pairs_anno_drivers.sql"))

df <- res %>% 
  select(case_barcode, idh_codel_subtype, cluster_id, size, ccf_a, ccf_b, rank_a, rank_b, drivers) %>%
  mutate(cat = case_when(ccf_a > 0.5 & ccf_b > 0.5 ~ "Neutral (Clonal)",
                         ccf_a < 0.5 & ccf_b < 0.5 ~ "Neutral (Subclonal)",
                         ccf_a < 0.5 & ccf_b > 0.5 ~ "Positive Selection",
                         ccf_a > 0.5 & ccf_b < 0.5 ~ "Negative Selection")) %>%
  gather(key = "k", value = "v", ccf_a, ccf_b, rank_a, rank_b) %>% 
  separate(k, into = c("var", "sample_type")) %>%
  spread(var, v) %>%
  mutate(id = sprintf("%s%s", case_barcode, cluster_id),
         sample_type = ifelse(sample_type == "a", "P", "R")) 

df_driver <- df %>%
  group_by(case_barcode) %>%
  mutate(has_driver = any(!is.na(drivers)), has_ccf = any(ccf > 0.90 & sample_type == "P") & any(ccf > 0.90 & sample_type == "R")) %>%
  ungroup() %>%
  filter(has_ccf) %>%
  group_by(case_barcode, cluster_id) %>%
  mutate(n = sum(ccf > 0.90), is_truncal = sum(ccf > 0.90) == 2) %>%
  ungroup() %>%
  #filter(n==2) %>%
  group_by(case_barcode) %>%
  summarize(has_trunk = any(n==2), trunk_driver = any(n==2 & !is.na(drivers)), 
            trunk_size = sum(size[is_truncal]), all_size = sum(size),
            rel_size = trunk_size / all_size) %>%
  ungroup()

tmp = df_driver %>% filter(all_size < 1e4, has_trunk)

# df2 <- res %>%
#   group_by(case_barcode) %>%
#   summarize(w_selected_ps = sum(size[ccf_a < 0.5 & ccf_b > 0.5]),
#             w_selected_ns = sum(size[ccf_a > 0.5 & ccf_b < 0.5]),
#             w_total = sum(size)) %>%
#   mutate(p_selection_pos = w_selected_ps / (w_total - w_selected_ns),
#          p_selection_all = (w_selected_ps + w_selected_ns) / w_total) %>%
#   ungroup()

g <- ggplot(df, aes(x=sample_type, y=ccf, size = cut(size,breaks = c(-Inf,2,10,30,Inf)), group = id, color = cat)) + 
  geom_point() + 
  geom_line() + 
  scale_size_manual(values = c(1,1.5,2,2.5)) +
  labs(x = "Sample Type", y = "Cancer Cell Fraction", size = "Cluster Size", color = "Selection") +
  theme_bw(base_size = 12) +
  geom_dl(aes(label = drivers), method = list(dl.trans(x = x + 0.2), "last.points", cex = 0.8)) +
  geom_dl(aes(label = drivers), method = list(dl.trans(x = x - 0.2), "first.points", cex = 0.8))

g + facet_wrap_paginate(~idh_codel_subtype+case_barcode, nrow = 4, ncol = 4, page = 1)

n_pages(g + facet_wrap_paginate(~idh_codel_subtype+case_barcode, nrow = 5, ncol = 4, page = 1))

pdf(file = "~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/pyclone_plots.pdf", width = 12, height = 12)
for(i in 1:11)
  plot(g + facet_wrap_paginate(~idh_codel_subtype+case_barcode, nrow = 5, ncol = 4, page = i))
dev.off()

## Plot rank tile plots

## First for all samples

df <- res %>%
  count(rank_a, rank_b) %>%
  complete(rank_a, rank_b, fill = list(n=0))

g <- ggplot(df, aes(x=rank_a, y=rank_b, fill = n)) + 
  geom_tile() +
  labs(x = "Cluster CCF Rank in Primary", y = "Cluster CCF Rank in Recurrence", fill = "Number of Samples") +
  scale_fill_distiller(palette = "PuRd", direction = 1) +
  coord_cartesian(ylim = c(1,10), xlim = c(1,10)) + 
  theme_bw(base_size = 12) + theme(axis.title = element_text(size = 12),
                                   axis.text = element_text(size=12),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   panel.background = element_rect(fill = "transparent"),
                                   axis.line = element_blank()) +
  scale_y_continuous(breaks = seq(2,10,2)) +
  scale_x_continuous(breaks = seq(2,10,2))

g

pdf(file = "~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/pyclone_rank_ccf_comparison_all.pdf", width = 5.5, height = 4)
plot(g)
dev.off()

## Next seperately by subtype

myplots <- lapply(split(res, res$idh_codel_subtype), function(df){
  df <- df %>%
    group_by(idh_codel_subtype) %>%
    count(rank_a, rank_b) %>%
    complete(rank_a, rank_b, fill = list(n=0)) %>%
    ungroup()
  
  p <- ggplot(df, aes(x=rank_a, y=rank_b, fill = n)) + 
    geom_tile() +
    facet_wrap(~idh_codel_subtype) +
    labs(x = "Cluster CCF Rank in Primary", y = "Cluster CCF Rank in Recurrence", fill = "Number of Samples") +
    scale_fill_distiller(palette = "PuRd", direction = 1) +
    coord_cartesian(ylim = c(1,10), xlim = c(1,10)) + 
    theme_bw(base_size = 12) + theme(axis.title = element_text(size = 12),
                                                      axis.text = element_text(size=12),
                                                      panel.grid.major = element_blank(),
                                                      panel.grid.minor = element_blank(),
                                                      panel.background = element_rect(fill = "transparent"),
                                                      axis.line = element_blank()) +
    scale_y_continuous(breaks = seq(2,10,2)) +
    scale_x_continuous(breaks = seq(2,10,2))
  
  return(p)
})

do.call(grid.arrange,(c(myplots, ncol=3)))

pdf(file = "~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/pyclone_rank_ccf_comparison_subtype.pdf", width = 17, height = 4)
do.call(grid.arrange,(c(myplots, ncol=3)))
dev.off()


g <- ggplot(res, aes(x=ccf_a, y=ccf_b, color = idh_codel_subtype, size=cut(size,breaks = c(-Inf,10,100,1000,Inf)))) +
  geom_point(alpha = 0.1) +
  labs(x = "Cancer Cell Fraction Rank in Primary", y = "Cancer Cell Fraction Rank in Recurrent", size = "Number of Samples", color = "Subtype") +
  theme_bw(base_size = 12) + theme(axis.title = element_text(size = 10),
                                   axis.text = element_text(size=10),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   panel.background = element_rect(fill = "transparent"),
                                   axis.line = element_blank())
g

pdf(file = "~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/pyclone_ccf_comparison_all.pdf", width = 6, height = 4)
plot(g)
dev.off()
