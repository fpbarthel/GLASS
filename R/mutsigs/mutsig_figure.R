library(DBI)
library(tidyverse)
library(MutationalPatterns)

con <- DBI::dbConnect(odbc::odbc(), "GLASSv2")  
res <- dbGetQuery(con, read_file("sql/mut_sig/mut_sig_effect.sql"))

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

