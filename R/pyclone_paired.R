
library(DBI)
library(tidyverse)
library(ggplot2)
library(ggforce)
library(egg)

con <- DBI::dbConnect(odbc::odbc(), "GLASSv2")
res <- dbGetQuery(con, q)#read_file("sql/pyclone/pyclone_cluster_pairs.sql"))

df <- res %>% 
  select(case_barcode, idh_codel_subtype, cluster_id, size, ccf_a, ccf_b, rank_a, rank_b) %>%
  mutate(cat = case_when(ccf_a > 0.5 & ccf_b > 0.5 ~ "Neutral (Clonal)",
                         ccf_a < 0.5 & ccf_b < 0.5 ~ "Neutral (Subclonal)",
                         ccf_a < 0.5 & ccf_b > 0.5 ~ "Positive Selection",
                         ccf_a > 0.5 & ccf_b < 0.5 ~ "Negative Selection")) %>%
  gather(key = "k", value = "v", ccf_a, ccf_b, rank_a, rank_b) %>% 
  separate(k, into = c("var", "sample_type")) %>%
  spread(var, v) %>%
  mutate(id = sprintf("%s%s", case_barcode, cluster_id),
         sample_type = ifelse(sample_type == "a", "P", "R"))

g <- ggplot(df, aes(x=sample_type, y=ccf, size = cut(size,breaks = c(-Inf,2,10,30,Inf)), group = id, color = cat)) + 
  geom_point() + 
  geom_line() + 
  scale_size_manual(values = c(1,1.5,2,2.5)) +
  labs(x = "Sample Type", y = "Cancer Cell Fraction", size = "Cluster Size", color = "Selection") +
  theme_bw(base_size = 12)

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

ggplot(res, aes(x=ccf_a, y=ccf_b, color = idh_codel_subtype, size=cut(size,breaks = c(-Inf,10,100,1000,Inf)))) +
  geom_point(alpha = 0.1) +
  labs(x = "Cancer Cell Fraction Rank in Primary", y = "Cancer Cell Fraction Rank in Recurrent", size = "Number of Samples", color = "Subtype") +
  theme_bw(base_size = 12) + theme(axis.title = element_text(size = 10),
                                   axis.text = element_text(size=10),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   panel.background = element_rect(fill = "transparent"),
                                   axis.line = element_blank())
  
