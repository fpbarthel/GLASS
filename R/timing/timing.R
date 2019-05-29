library(DBI)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(egg)

con <- DBI::dbConnect(odbc::odbc(), "GLASSv2")
res <- dbGetQuery(con, read_file('sql/timing/timing.sql'))
pairs <- dbGetQuery(con, read_file('sql/timing/timing.sql'))

## CCF compare

testWilcoxGroup <- function(df) {
  wtest = wilcox.test(df$ccf ~ df$evnt, paired = TRUE, conf.int = TRUE)
  data.frame(n = nrow(df)/2,
             statistic = wtest$statistic,
             estimate = wtest$estimate,
             lcl = wtest$conf.int[1],
             ucl = wtest$conf.int[2],
             wilcox_p = wtest$p.value,
             test_str = sprintf("n=%s\n%s",
                                nrow(df)/2,
                                case_when(wtest$p.value < 0.0001 ~ "P<0.0001",
                                          wtest$p.value > 0.05 ~ sprintf("P=%s", format(round(wtest$p.value, 2),scientific=F)),
                                          TRUE ~ sprintf("P=%s", format(round(wtest$p.value, 4),scientific=F)))),
             stringsAsFactors = FALSE)
}

grobs <- lapply(1:nrow(pairs), function(i) {
  tmp <- res %>% filter(subtype == pairs$idh_codel_subtype[i], evnt %in% c(pairs$evnt_a[i], pairs$evnt_b[i])) %>% 
    group_by(aliquot_barcode) %>% 
    mutate(n=n()) %>% 
    ungroup() %>%
    filter(n==2)
  
  if(nrow(tmp)/2 < 4)
    return()
  
  wtest_all <- testWilcoxGroup(tmp)
  wtest_sam <- tmp %>% group_by(sample_type) %>% do(testWilcoxGroup(.)) %>% ungroup()
  
  ggplot(tmp, aes(x=evnt, y=ccf, group = aliquot_barcode, color = sample_type)) + 
    geom_point() +
    geom_line() + 
    theme_bw() +
    coord_cartesian(ylim = c(0,1), xlim = c(1,2)) +
    geom_text(data = wtest_all, aes(x=1.5, y=0.05, label = test_str), group = NA, color = "black", size = 3) +
    geom_text(data = filter(wtest_sam, sample_type == "P"), aes(x=1, y=0.05, label = test_str), group = NA, color = "#007C80", size = 3) +
    geom_text(data = filter(wtest_sam, sample_type == "R"), aes(x=2, y=0.05, label = test_str), group = NA, color = "#CE8014", size = 3) +
    guides(color = FALSE) +
    labs(x = "Event", y = "Cancer Cell Fraction", title = pairs$idh_codel_subtype[i]) +
    scale_color_manual(values = c("P" = "#007C80", "R" = "#CE8014"))
})

grobs = grobs[which(!sapply(grobs,is.null))]

pdf("~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/timing.pdf", width=9, height = 12, useDingbats = FALSE)
marrangeGrob(grobs = grobs, nrow = 4, ncol = 3)
dev.off()


tmp <- res %>% filter(subtype == "IDHwt", evnt %in% c("TP53 mut", "NF1 mut")) %>% 
  group_by(aliquot_barcode) %>% 
  mutate(n=n()) %>% 
  ungroup() %>%
  filter(n==2)

ggplot(tmp, aes(x=evnt, y=ccf, group = aliquot_barcode)) + geom_point() + geom_line() + facet_wrap(~sample_type)

## Proportion clonal radar plot

myplots <- lapply(split(res, res$subtype), function(df){
  df <- df %>% 
    arrange(desc(prop_clonal)) %>%
    group_by(evnt) %>%
    mutate(label = sprintf("%s\n(n=%s/%s)", evnt, n_total[sample_type == 'P'], n_total[sample_type == 'R'])) %>%
    ungroup()
  df$label = factor(df$label, levels = unique(df$label))
  df = df[order(df$label),]

  p <- ggplot(df, aes(x=label, y=prop_clonal, group = sample_type, color = sample_type, fill = sample_type)) +
    geom_point(size=2) + 
    geom_polygon(size = 1, alpha= 0.2) + 
    ylim(0,1) +
    labs(y = "Proportion clonal", x = "Alteration", fill = "Sample Type", color = "Sample Type") +
    scale_x_discrete() +
    theme_bw(base_size = 12) +
    coord_polar() +
    guides(fill = FALSE, color = FALSE)
  
  return(p)
})

do.call(grid.arrange,(c(myplots, ncol=3)))

## Proportion clonal barplot

df <- res %>% 
  mutate(n_subclonal = n_total - n_clonal) %>%
  group_by(subtype, evnt) %>%
  mutate(label = sprintf("%s\n(n=%s/%s)", evnt, n_total[sample_type == 'P'], n_total[sample_type == 'R'])) %>%
  ungroup() %>%
  select(evnt,label,subtype,sample_type,n_clonal,n_subclonal) %>%
  gather(key = "clonality", "count", n_clonal, n_subclonal) %>%
  mutate(clonality = factor(clonality, levels = c("n_clonal", "n_subclonal"), labels = c("Clonal", "Subclonal")))
  
df$label = factor(df$label, levels = unique(df$label))
df = df[order(df$label),]

ggplot(df, aes(x=sample_type, y = count, fill = clonality)) + geom_bar(stat = "identity") + facet_grid(subtype~evnt, scales = "free")

## Primary vs recurrent CCF

ggplot(res, aes(x=evnt, color=sample_type, y=ccf)) +
  geom_point(position=position_jitterdodge()) +
  geom_boxplot(alpha=0.6) +
  facet_wrap(~subtype, scales = "free_x", ncol = 1)

ggplot(res, aes(x=sample_type, color=subtype, y=ccf)) +
  geom_point(position=position_jitterdodge()) +
  geom_boxplot(alpha=0.6) +
  facet_wrap(~evnt) +
  theme_bw() + 
  labs(x = "Sample Type", y = "Cancer Cell Fraction")

## Barplots for rank per gene and subtype (counts)

res2 <- res %>% mutate(subtype = factor(subtype),
                      evnt_rank = factor(evnt_rank, levels = 1:7, labels = c("1st", "2nd", "3rd", sprintf("%sth",4:7))))

myplots <- lapply(split(res2, list(res2$subtype, res2$sample_type)), function(df){
  df <- df %>% 
    group_by(evnt) %>%
    mutate(mean_ccf = mean(ccf),
           num_evnt = n()) %>%
    ungroup() %>%
    arrange(desc(num_evnt), desc(mean_ccf), evnt) %>%
    mutate(evnt = factor(gsub(" ", "\n", evnt), levels = unique(gsub(" ", "\n", evnt))))
  
  mytitle = sprintf("%s (%s)", unique(df$subtype), ifelse(unique(df$sample_type) == "P", "Primary", "Recurrent"))
  
  p <- ggplot(df, aes(x=evnt, fill=evnt_rank)) +
    geom_bar() + 
    labs(y = "Count", x = "Alteration", fill = "Event Rank", title = mytitle) + 
    scale_fill_manual(values = c("#FFEC8B", "#AFFF8B", "#8BFFCB", "#8BD0FF", "#AA8BFF", "#FF8BF1", "#FF8B8D"), drop = F) +
    theme_bw(base_size = 12) +
    guides(fill = FALSE, color = FALSE)
  
  return(p)
})

grid.arrange(myplots[[1]], myplots[[2]], myplots[[3]], myplots[[4]], myplots[[5]], myplots[[6]], ncol=3, widths = c(6,8,17))

## barplots for rank per gene and subtype (median)

res <- res %>% mutate(subtype = factor(subtype))

myplots <- lapply(split(res, list(res$subtype, res$sample_type)), function(df){
  df <- df %>% 
    arrange(med_rank) %>%
    mutate(evnt = factor(gsub(" ", "\n", evnt), levels = unique(gsub(" ", "\n", evnt))))
  
  mytitle = sprintf("%s (%s)", unique(df$subtype), ifelse(unique(df$sample_type) == "P", "Primary", "Recurrent"))
  
  p <- ggplot(df, aes(x=evnt, y=med_rank)) +
    geom_bar(stat="identity", fill = "#FFEC8B") + 
    geom_text(aes(label = sprintf("n=%s", num_events)), position = position_dodge(width=0.9), size=4) +
    ylim(0,5) +
    labs(y = "Median Rank", x = "Alteration", title = mytitle) +
    scale_x_discrete() +
    theme_bw(base_size = 12) +
    guides(fill = FALSE, color = FALSE)
  
  return(p)
})

grid.arrange(myplots[[1]], myplots[[2]], myplots[[3]], myplots[[4]], myplots[[5]], myplots[[6]], ncol=3, widths = c(6,8,17))

## END ##