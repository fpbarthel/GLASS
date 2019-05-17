library(DBI)
library(tidyverse)
library(MutationalPatterns)

con <- DBI::dbConnect(odbc::odbc(), "GLASSv2") 
res <- dbGetQuery(con, read_file("sql/figures/mutsig_boxplot_fig1.sql"))

################ stacked bar plots

res <- res %>% filter(rnk == 1, all_fractions_counts == 3) %>%
  mutate(signature = factor(signature, levels = c(1,3,8,11,16), labels = c("Signature 1\nAging","Signature 3\nDNA DSB repair","Signature 8\nUnknown","Signature 11\nAlkylating Agent","Signature 16\nUnknown")),
         fraction = factor(fraction, levels = c("P","S","R")))

plot_theme <- theme_bw(base_size = 12) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size=12),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        axis.line = element_blank())

gg_sig_stackedbar <- ggplot(res, aes(x = fraction, fill = factor(signature))) +
  geom_bar(width=1, color = "black") +
  facet_wrap(~idh_codel_subtype, scales = "free_y") +
  scale_fill_manual(values = c("#e41a1c","#4daf4a","#377eb8","#984ea3","#ff7f00")) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) +
  plot_theme + 
  labs(x = "Fraction", y = "Number of Patients (n=216)", fill = "Dominant\nMutational\nSignature")

pdf("~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/Figure 1/mutsigs_stacked_bars.pdf", width = 8, height  = 6, useDingbats = FALSE)
plot(gg_sig_stackedbar)
dev.off()




################ ################ ################ ################ ################ ################ ################ 
## New figure
################ ################ ################ ################ ################ ################ ################ 

res <- dbGetQuery(con, read_file("sql/figures/mutsig_boxplot_fig1.sql"))

res <- res %>% mutate(signature = factor(signature, levels = c(16,3,8,11,1), labels = c("Signature 16\nUnknown","Signature 3\nDNA DSB repair","Signature 8\nUnknown","Signature 11\nAlkylating Agent","Signature 1\nAging")),
                      fraction = factor(fraction, levels = c("R","S","P"))) #sd_rel_score = ifelse(avg_rel_score + sd_rel_score > 1, 1 - avg_rel_score, sd_rel_score))

tmp0 <- res %>% filter(hypermutator_status == 0)
tmp1 <- res %>% filter(hypermutator_status == 1)

plot_theme <- theme_bw(base_size = 12) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size=12),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        axis.line = element_blank(),
        legend.position = 'none',
        axis.title.y=element_blank()) + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())

plot_colors <- scale_fill_manual(values=c("#2FB3CA","#CA932F","#CA2F66"))

null_y        <- theme(axis.text.y=element_blank())

g0_box <- ggplot(tmp0, aes(x = signature, y = rel_score, fill = fraction)) + 
  geom_boxplot() +
  scale_y_reverse(expand=c(0,0)) +
  scale_x_discrete(position = "top") +
  coord_flip(ylim=c(0,1.0)) +
  plot_theme + null_y + plot_colors +
  labs(y = "Signature Score (Non-hypermutants)") +
  theme(plot.margin= unit(c(1, 0.2, 1, 1), "lines")) ## Top, Right, Bottom, Left

g1_box <- ggplot(tmp1, aes(x = signature, y = rel_score, fill = fraction)) + 
  geom_boxplot() +
  scale_y_continuous(expand=c(0,0)) + 
  scale_x_discrete() +
  coord_flip(ylim=c(0,1.0)) +
  plot_theme + plot_colors +
  theme(axis.text.y = element_text(angle=0, vjust = 0.5)) +
  labs(y = "Signature Score (Hypermutants)") +
  theme(plot.margin= unit(c(1, 1, 1, 0.2), "lines")) ## Top, Right, Bottom, Left

g0_bar <- ggplot(tmp0, aes(x = signature, y = avg_rel_score, fill = fraction, ymax = avg_rel_score + sd_rel_score, ymin = 0)) + 
  geom_errorbar(position = position_dodge(width = 0.9), width = 0.5) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_y_reverse(expand=c(0,0)) +
  scale_x_discrete(position = "top") +
  coord_flip(ylim=c(0,1.0)) +
  plot_theme + null_y + plot_colors +
  labs(y = "Signature Score (Non-hypermutants)") +
  theme(plot.margin= unit(c(1, 0.2, 1, 1), "lines")) ## Top, Right, Bottom, Left

g1_bar <- ggplot(tmp1, aes(x = signature, y = avg_rel_score, fill = fraction, ymax = avg_rel_score + sd_rel_score, ymin = 0)) + 
  geom_errorbar(position = position_dodge(width = 0.9), width = 0.5) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_y_continuous(expand=c(0,0)) + 
  scale_x_discrete() +
  coord_flip(ylim=c(0,1.0)) +
  plot_theme + plot_colors +
  theme(axis.text.y = element_text(angle=0, vjust = 0.5)) +
  labs(y = "Signature Score (Hypermutants)") +
  theme(plot.margin= unit(c(1, 1, 1, 0.2), "lines")) ## Top, Right, Bottom, Left

pdf("~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/Figure 1/mutsigs_main_boxplot.pdf", width = 8, height  = 6, useDingbats = FALSE)
grid.arrange(g0_box, g1_box, ncol=2, widths = c(0.725,1))
dev.off()

pdf("~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/Figure 1/mutsigs_main_barplot.pdf", width = 8, height  = 6, useDingbats = FALSE)
grid.arrange(g0_bar, g1_bar, ncol=2, widths = c(0.725,1))
dev.off()


message("Stop here")

### Below are former attempts
# res <- dbGetQuery(con, read_file("sql/mut_sig/mut_sig_effect.sql"))
ggplot(res, aes(y=rel_score, x=factor(signature), fill = variant_effect)) + geom_bar(stat = "identity",position="dodge" ) + facet_grid(known_driver_status~sample_type)

mat <- res %>% 
  filter(variant_effect == 'Coding', mut_n_total > 9) %>%
  mutate(subst = sprintf("%s[%s>%s]%s", substr(trinucleotide_context,1,1), substr(trinucleotide_context,2,2), alt, substr(trinucleotide_context,3,3)),
         set = sprintf("%s\n%s\n%s\n(n=%s)", ifelse(variant_effect == "Coding", "Non-synonymous", "Synonymous/intronic"), gene_symbol, idh_codel_subtype, mut_n_total)) %>%
  select(subst,set,mut_n) %>%
  arrange(subst) %>%
  spread(set, mut_n, fill = 0)

mat <- res %>%
  filter(variant_effect == 'Non-coding') %>%
  mutate(subst = sprintf("%s[%s>%s]%s", substr(trinucleotide_context,1,1), substr(trinucleotide_context,2,2), alt, substr(trinucleotide_context,3,3)),
         set = sprintf("%s\n%s\n%s", "Synonymous", "Intronic", idh_codel_subtype)) %>%
  group_by(set) %>%
  mutate(mut_n_total = sum(mut_n)) %>%
  ungroup() %>%
  mutate(set = sprintf("%s\n(n=%s)", set, mut_n_total)) %>%
  group_by(subst,set) %>%
  summarize(mut_n = sum(mut_n)) %>%
  ungroup() %>%
  select(subst,set,mut_n) %>%
  arrange(subst) %>%
  spread(set, mut_n, fill = 0)
  

tmp = mat[,1,drop=T]
mat = as.matrix(mat[,-1])
rownames(mat) = tmp

plot_96_profile(mat)


res <- dbGetQuery(con, read_file("sql/mut_sig/mut_sig_fraction_subtype_hypermutation.sql"))

mat <- res %>%
  filter(complete.cases(hypermutator_status)) %>%
  mutate(subst = sprintf("%s[%s>%s]%s", substr(trinucleotide_context,1,1), substr(trinucleotide_context,2,2), alt, substr(trinucleotide_context,3,3)),
         set = sprintf("%s\n%s\n%s", ifelse(hypermutator_status == "0", "Non-hyperm.", "Hyperm."), idh_codel_subtype, fraction)) %>%
  group_by(set) %>%
  mutate(mut_n_total = sum(mut_n)) %>%
  ungroup() %>%
  mutate(set = sprintf("%s\n(n=%s)", set, mut_n_total)) %>%
  group_by(subst,set) %>%
  summarize(mut_n = sum(mut_n)) %>%
  ungroup() %>%
  select(subst,set,mut_n) %>%
  arrange(subst) %>%
  spread(set, mut_n, fill = 0)

tmp = mat[,1,drop=T]
mat = as.matrix(mat[,-1])
rownames(mat) = tmp

plot_96_profile(mat, condensed = F)

res <- dbGetQuery(con, read_file("sql/mut_sig/mut_sig_fraction_subtype_hypermutation.sql"))

df <- res %>% group_by(signature) %>% mutate(mean_sig = mean(rel_score), max_sig = max(rel_score)) %>% ungroup() %>% filter(max_sig > 0.1) %>% arrange(desc(rel_score)) %>% mutate(signature = factor(signature, levels = unique(signature)))

ggplot(df, aes(x=signature, y=rel_score, color = fraction, group = fraction, fill = fraction)) +
  geom_point(size=2) + 
  geom_polygon(size=1,alpha=0.2) +
  ylim(0,1) +
  scale_x_discrete() + 
  coord_polar() + 
  theme_bw() +
  facet_grid(hypermutator_status ~ idh_codel_subtype)

ggplot(df, aes(x=signature, y=rel_score, color = fraction, group = fraction, fill = fraction)) +
  geom_bar(size=1,alpha=0.2, stat="identity", position = "dodge") +
  ylim(0,1) +
  scale_x_discrete() + #coord_polar() + 
  theme_bw() +
  facet_grid(hypermutator_status ~ idh_codel_subtype)

