
## To get paper figures, run:
##   - Section titled: "DB connection"
##   - Section titled: "Run dNdS for all samples in cohort"
##   - Section titled: "Run dNdS seperately for private/shared variants and IDH status"
##   - Section titled: "Plot dNdS-CV by subtype - by fraction
##   - Section titled: "Plot dNdS-CV genes by subtype - fraction"

library(dndscv)
library(tidyverse)
library(DBI)
library(ggthemes)
library(gridExtra)
library(ggplot2)

data("cancergenes_cgc81", package="dndscv")

###################
## DB connection
###################

con <- DBI::dbConnect(odbc::odbc(), "GLASSv2")

res_fraction <- dbGetQuery(con, read_file("sql/dndscv/dndscv_input_by_fraction.sql"))
res_sample <- dbGetQuery(con, read_file("sql/dndscv/dndscv_input_by_sample.sql"))

###################
## Run dNdS for all samples in cohort
###################

dnds_all <- res_fraction %>%
  select(case_barcode,chrom,pos,ref,mut) %>% 
  distinct() %>%
  dndscv(refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = NULL)

message("Ran dNdSCV for ", n_distinct(res_fraction$case_barcode) - length(dnds_all$exclsamples))

dnds_all_global = dnds_all$globaldnds
dnds_all_sel_cv = dnds_all$sel_cv
dnds_all_gene_ci = geneci(dnds_all, gene_list = known_cancergenes)

###################
## Examine results
###################
print(dnds_all$globaldnds)
print(dnds_all$nbreg$theta)

sel_cv = dnds_all$sel_cv
print(head(sel_cv, sum(sel_cv$qglobal_cv < 0.05),n=25), digits = 3)

###################
## Run dNdS seperately for private/shared variants and IDH status
###################

result_list <- lapply(na.omit(unique(res_fraction$subtype)), function(st) {
  result_list <- lapply(na.omit(unique(res_fraction$fraction)), function(fr) {
    message("Computing dNdS for ", fr, " and ", st)
    qres_subset = res_fraction %>% filter(fraction == fr, subtype == st) %>% select(case_barcode,chrom,pos,ref,mut) %>% distinct()
    dnds_subset = dndscv(qres_subset, refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = NULL)
    message(".. Ran dNdSCV for ", n_distinct(qres_subset$case_barcode) - length(dnds_subset$exclsamples))
    message(".. Found ", nrow(dnds_subset$annotmuts), " mutations")
    globaldnds <- cbind(fraction = fr, subtype = st, dnds_subset$globaldnds)
    sel_cv <- cbind(fraction = fr, subtype = st, dnds_subset$sel_cv)
    gene_ci <- cbind(fraction = fr, subtype = st, geneci(dnds_subset, gene_list = known_cancergenes))
    list(globaldnds, sel_cv, gene_ci)
  })
  res1 <- data.table::rbindlist(lapply(result_list,'[[',1))
  res2 <- data.table::rbindlist(lapply(result_list,'[[',2))
  res3 <- data.table::rbindlist(lapply(result_list,'[[',3))
  list(res1,res2,res3)
})
dnds_fraction_global <- data.table::rbindlist(lapply(result_list,'[[',1))
dnds_fraction_sel_cv <- data.table::rbindlist(lapply(result_list,'[[',2))
dnds_fraction_gene_ci <- data.table::rbindlist(lapply(result_list,'[[',3))
print(dnds_fraction_global)
print(dnds_fraction_sel_cv)

###################
## Run dNdS seperately for all primaries and recurrences and IDH status
###################

result_list <- lapply(na.omit(unique(res_sample$subtype)), function(st) {
  result_list <- lapply(na.omit(unique(res_sample$sample_type)), function(sampt) {
    message("Computing dNdS for ", sampt, " and ", st)
    qres_subset = res_sample %>% filter(sample_type == sampt, subtype == st) %>% select(aliquot_barcode,chrom,pos,ref,mut) %>% distinct()
    dnds_subset = dndscv(qres_subset, refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = NULL)
    globaldnds <- cbind(sample_type = sampt, subtype = st, dnds_subset$globaldnds)
    sel_cv <- cbind(sample_type = sampt, subtype = st, dnds_subset$sel_cv)
    gene_ci <- cbind(sample_type = sampt, subtype = st, geneci(dnds_subset, gene_list = known_cancergenes))
    list(globaldnds, sel_cv, gene_ci)
  })
  res1 <- data.table::rbindlist(lapply(result_list,'[[',1))
  res2 <- data.table::rbindlist(lapply(result_list,'[[',2))
  res3 <- data.table::rbindlist(lapply(result_list,'[[',3))
  list(res1,res2,res3)
})
dnds_sample_global <- data.table::rbindlist(lapply(result_list,'[[',1))
dnds_sample_sel_cv <- data.table::rbindlist(lapply(result_list,'[[',2))
dnds_sample_gene_ci <- data.table::rbindlist(lapply(result_list,'[[',3))
print(dnds_sample_global)
print(dnds_sample_sel_cv)

###################
## Plot dNdS for all data
###################

p <- dnds_all_global %>%
  ggplot(aes(x = name, y = mle, ymin = cilow, ymax = cihigh, color = name)) +
  geom_pointrange(position=position_dodge(width = 0.5)) +
  geom_hline(yintercept = 1) +
  theme_bw(base_size = 18) +
  labs(y = "Global dN/dS", color = "Variant Type")
plot(p)

###################
## Plot all genes
###################

df = dnds_all_sel_cv[1:25,]
df$gene_name = factor(df$gene_name, levels = unique(df$gene_name))
p<-ggplot(df, aes(x=gene_name, y=-log10(qglobal_cv)), width=0.75) +
  geom_bar(aes(y=-log10(pglobal_cv)), fill = "gray90", stat="identity", position = position_dodge()) + 
  geom_bar(stat="identity", position = position_dodge()) + 
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = 2) +
  labs(y="-log10(FDR)", x="") +
  guides(fill = F) +
  coord_flip(ylim = c(0,10)) +
  scale_fill_manual(values = c("P" = "#2fb3ca","S" = "#CA932F", "R" = "#CA2FB4"), drop=F) +
  theme_bw(base_size = 16) +
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.ticks = element_line(size = 0))
plot(p)

###################
## Plot dNdS-CV by subtype - by fraction
###################

p<-dnds_fraction_global %>%
  filter(name == 'wall') %>%
  mutate(subtype = factor(subtype, levels = c("IDHmut-codel","IDHmut-noncodel","IDHwt")),
         fraction = factor(fraction, levels = c("P", "S", "R")),
         fraction = fct_recode(fraction, "Primary only" = "P", "Shared" = "S", "Recurrence only" = "R")) %>%
  ggplot(aes(x = subtype, y = mle, ymin = cilow, ymax = cihigh, color = fraction)) +
  geom_pointrange(position=position_dodge(width = 0.5)) +
  scale_color_manual(values = c("Primary only" = "#CA2FB4", "Shared" = "#CA932F", "Recurrence only" = "#2fb3ca"), drop=F) +
  #facet_wrap( ~ subtype) +
  geom_hline(yintercept = 1) +
  theme_bw(base_size = 10) +
  coord_cartesian(ylim=c(0.75,2)) +
  labs(x = "Subtype", y = "Global dN/dS", color = "Variant Fraction")

pdf('~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/F2/dnds-fraction.pdf', height = 4, width = 6)
plot(p)
dev.off()

 ###################
## Plot dNdS-CV by subtype - by sample
###################

p<-dnds_sample_global %>%
  filter(name == 'wall') %>%
  mutate(subtype = factor(subtype, levels = c("IDHmut-codel","IDHmut-noncodel","IDHwt")),
         sample_type = factor(sample_type, levels = c("P", "R")),
         sample_type = fct_recode(sample_type, "Primary" = "P", "Recurrence" = "R")) %>%
  ggplot(aes(x = subtype, y = mle, ymin = cilow, ymax = cihigh, color = sample_type)) +
  geom_pointrange(position=position_dodge(width = 0.5)) +
  scale_color_manual(values = c("Primary" = "#a6611a", "Recurrence" = "#018571"), drop=F) +
  #facet_wrap( ~ subtype) +
  geom_hline(yintercept = 1) +
  theme_bw(base_size = 10) +
  coord_cartesian(ylim=c(0.75,2)) +
  labs(x = "Subtype", y = "Global dN/dS", color = "Sample Type")

pdf('~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/F2/dnds-sample.pdf', height = 4, width = 6)
plot(p)
dev.off()

###################
## Plot dNdS-CV genes by subtype - fraction
###################

dat <- dnds_fraction_sel_cv %>%
  mutate(fraction = factor(fraction, levels = c("P", "S", "R"))) %>%
  complete(gene_name,fraction,subtype, fill = list(qglobal_cv=1)) %>%
  group_by(fraction,subtype) %>%
  arrange(qglobal_cv, pglobal_cv) %>%
  ungroup() %>%
  arrange(fraction)

myplots <- lapply(split(dat, paste(dat$subtype, dat$fraction))[c(1,3,2,4,6,5,7,9,8)], function(df){
  df$gene_name = factor(df$gene_name, levels = unique(df$gene_name))
  p <- ggplot(df[1:7,], aes(x=gene_name, y=-log10(qglobal_cv), fill=fraction), width=0.75) +
    geom_bar(aes(y=-log10(pglobal_cv)), fill = "gray90", stat="identity", position = position_dodge()) + 
    geom_bar(stat="identity", position = position_dodge()) + 
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = 2) +
    labs(y="-log10(FDR)", x="") +
    guides(fill = F) +
    facet_wrap(~ subtype, scales = "free_y") + 
    coord_flip(ylim = c(0,10)) +
    scale_fill_manual(values = c("P" = "#2fb3ca","S" = "#CA932F", "R" = "#CA2FB4"), drop=F) +
    theme_bw(base_size = 16) +
    theme( panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           axis.ticks = element_line(size = 0))
  
  return(p)
})


pdf('~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/EDF7/dnds-genes.pdf', height = 8, width = 12)
do.call(grid.arrange,(c(myplots, ncol=3)))
dev.off()

###################
## Plot dNdS-CV genes by subtype - samples
###################

dat <- dnds_sample_sel_cv %>%
  mutate(sample_type = factor(sample_type, levels = c("P", "R"))) %>%
  complete(gene_name,sample_type,subtype, fill = list(qglobal_cv=1)) %>%
  group_by(sample_type,subtype) %>%
  arrange(qglobal_cv, pglobal_cv) %>%
  ungroup()

myplots <- lapply(split(dat, paste(dat$subtype, dat$sample_type)), function(df){
  df$gene_name = factor(df$gene_name, levels = unique(df$gene_name))
  p <- ggplot(df[1:9,], aes(x=gene_name, y=-log10(qglobal_cv), fill=sample_type), width=0.75) +
    geom_bar(aes(y=-log10(pglobal_cv)), fill = "gray90", stat="identity", position = position_dodge()) + 
    geom_bar(stat="identity", position = position_dodge()) + 
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = 2) +
    labs(x="", y="-log10(FDR)") +
    guides(fill = F) +
    facet_wrap(~ subtype, scales = "free_y") + 
    coord_flip(ylim = c(0,10)) +
    scale_fill_manual(values = c("P" = "#2fb3ca", "R" = "#CA2FB4"), drop=F) +
    theme_bw(base_size = 16) +
    theme( panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           axis.ticks = element_line(size = 0))
  
  return(p)
})

do.call(grid.arrange,(c(myplots, ncol=2)))

## Save image
save.image(file="results/dndscv/dndscv.RData")

library(scales)

## Join sel_cv and gene_ci for all
tmp <- dnds_all_sel_cv %>% 
  rename(gene = gene_name) %>% 
  inner_join(dnds_all_gene_ci) %>%
  filter(pglobal_cv < 0.05) %>%
  gather("k","v",mis_mle:tru_high) %>%
  separate(k, into = c("variant_type","k")) %>%
  spread(k, v) %>% 
  mutate(n = case_when(variant_type == "mis" ~ n_mis, variant_type == "tru" ~ n_non + n_spl + n_ind, TRUE ~ NA_real_)) %>%
  arrange(mle) %>%
  mutate(gene = factor(gene, levels = unique(gene)),
         variant_type = case_when(variant_type == "mis" ~ "Missense", variant_type == "tru" ~ "Truncating", TRUE ~ NA_character_))

ggplot(tmp, aes(x = gene, y = mle, ymin = low, ymax = high, color = variant_type, alpha = n)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_pointrange(position = position_dodge(width=1)) +
  coord_flip() +
  scale_y_log10(breaks = c(1e-1,1e0,1e1,1e2,1e3),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_few(base_size = 16) +
  labs(alpha = "Non-synonymous\nmutation count", color = "Variant classification", x = "Gene Symbol")

## Join sel_cv and gene_ci for IDHwt R fraction
tmp <- dnds_fraction_sel_cv %>%
  rename(gene = gene_name) %>% 
  inner_join(dnds_fraction_gene_ci) %>%
  gather("k","v",mis_mle:tru_high) %>%
  separate(k, into = c("variant_type","k")) %>%
  spread(k, v) %>% 
  mutate(n = case_when(variant_type == "mis" ~ n_mis, variant_type == "tru" ~ n_non + n_spl + n_ind, TRUE ~ NA_real_),
         p = case_when(variant_type == "mis" ~ pmis_cv, variant_type == "tru" ~ ptrunc_cv, TRUE ~ NA_real_)) %>%
  filter(p < 0.05) %>%
  arrange(mle) %>%
  mutate(gene = factor(gene, levels = unique(gene)),
         variant_type = case_when(variant_type == "mis" ~ "Missense", variant_type == "tru" ~ "Truncating", TRUE ~ NA_character_))

ggplot(tmp, aes(x = gene, y = mle, ymin = low, ymax = high, color = variant_type, alpha = n)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_pointrange(position = position_dodge(width=1)) +
  coord_flip() +
  scale_y_log10(breaks = c(1e-1,1e0,1e1,1e2,1e3),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_few(base_size = 16) +
  labs(alpha = "Non-synonymous\nmutation count", color = "Variant classification", x = "Gene Symbol") + 
  facet_wrap(subtype ~ fraction)



###################
## Write to database
###################

#dbWriteTable(con, Id(schema="analysis",table="dnds_sample_global"), dnds_sample_global, append = FALSE)
#dbWriteTable(con, Id(schema="analysis",table="dnds_fraction_global"), dnds_fraction_global, append = FALSE)
#dbWriteTable(con, Id(schema="analysis",table="dnds_sample_sel_cv"), dnds_sample_sel_cv, append = FALSE)
#dbWriteTable(con, Id(schema="analysis",table="dnds_fraction_sel_cv"), dnds_fraction_sel_cv, append = FALSE)
#dbWriteTable(con, Id(schema="analysis",table="dnds_neutralitytestr_fraction_global"), dnds_neutralitytestr_fraction_global, append = FALSE)
#dbWriteTable(con, Id(schema="analysis",table="dnds_neutralitytestr_sample_global"), dnds_neutralitytestr_sample_global, append = FALSE)


###################
## Re-run for hypermutators
###################

con <- DBI::dbConnect(odbc::odbc(), "GLASSv2")

res_fraction <- dbGetQuery(con, read_file("sql/dndscv/dndscv_input_by_fraction_hyperm.sql"))
res_sample <- dbGetQuery(con, read_file("sql/dndscv/dndscv_input_by_sample_hyperm.sql"))

###################
## Run dNdS for all samples in cohort
###################

dnds_all <- res_fraction %>%
  select(case_barcode,chrom,pos,ref,mut) %>% 
  distinct() %>%
  dndscv(refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = NULL)

message("Ran dNdSCV for ", n_distinct(res_fraction$case_barcode) - length(dnds_all$exclsamples))

dnds_all_global = dnds_all$globaldnds
dnds_all_sel_cv = dnds_all$sel_cv
dnds_all_gene_ci = geneci(dnds_all, gene_list = known_cancergenes)

###################
## Examine results
###################
print(dnds_all$globaldnds)
print(dnds_all$nbreg$theta)

sel_cv = dnds_all$sel_cv
print(head(sel_cv, sum(sel_cv$qglobal_cv < 0.05),n=25), digits = 3)

###################
## Run dNdS seperately for private/shared variants and IDH status
###################

result_list <- lapply(na.omit(unique(res_fraction$subtype)), function(st) {
  result_list <- lapply(na.omit(unique(res_fraction$fraction)), function(fr) {
    message("Computing dNdS for ", fr, " and ", st)
    qres_subset = res_fraction %>% filter(fraction == fr, subtype == st) %>% select(case_barcode,chrom,pos,ref,mut) %>% distinct()
    dnds_subset = dndscv(qres_subset, refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = NULL)
    message(".. Ran dNdSCV for ", n_distinct(qres_subset$case_barcode) - length(dnds_subset$exclsamples))
    globaldnds <- cbind(fraction = fr, subtype = st, dnds_subset$globaldnds)
    sel_cv <- cbind(fraction = fr, subtype = st, dnds_subset$sel_cv)
    gene_ci <- cbind(fraction = fr, subtype = st, geneci(dnds_subset, gene_list = known_cancergenes))
    list(globaldnds, sel_cv, gene_ci)
  })
  res1 <- data.table::rbindlist(lapply(result_list,'[[',1))
  res2 <- data.table::rbindlist(lapply(result_list,'[[',2), fill = TRUE)
  res3 <- data.table::rbindlist(lapply(result_list,'[[',3))
  list(res1,res2,res3)
})
dnds_fraction_global <- data.table::rbindlist(lapply(result_list,'[[',1))
dnds_fraction_sel_cv <- data.table::rbindlist(lapply(result_list,'[[',2))
dnds_fraction_gene_ci <- data.table::rbindlist(lapply(result_list,'[[',3))
print(dnds_fraction_global)
print(dnds_fraction_sel_cv)

###################
## Run dNdS seperately for all primaries and recurrences and IDH status
###################

result_list <- lapply(na.omit(unique(res_sample$subtype)), function(st) {
  result_list <- lapply(na.omit(unique(res_sample$sample_type)), function(sampt) {
    message("Computing dNdS for ", sampt, " and ", st)
    qres_subset = res_sample %>% filter(sample_type == sampt, subtype == st) %>% select(aliquot_barcode,chrom,pos,ref,mut) %>% distinct()
    dnds_subset = dndscv(qres_subset, refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = NULL)
    globaldnds <- cbind(sample_type = sampt, subtype = st, dnds_subset$globaldnds)
    sel_cv <- cbind(sample_type = sampt, subtype = st, dnds_subset$sel_cv)
    gene_ci <- cbind(sample_type = sampt, subtype = st, geneci(dnds_subset, gene_list = known_cancergenes))
    list(globaldnds, sel_cv, gene_ci)
  })
  res1 <- data.table::rbindlist(lapply(result_list,'[[',1))
  res2 <- data.table::rbindlist(lapply(result_list,'[[',2))
  res3 <- data.table::rbindlist(lapply(result_list,'[[',3))
  list(res1,res2,res3)
})
dnds_sample_global <- data.table::rbindlist(lapply(result_list,'[[',1))
dnds_sample_sel_cv <- data.table::rbindlist(lapply(result_list,'[[',2))
dnds_sample_gene_ci <- data.table::rbindlist(lapply(result_list,'[[',3))
print(dnds_sample_global)
print(dnds_sample_sel_cv)

###################
## Plot dNdS for all data
###################

p <- dnds_all_global %>%
  ggplot(aes(x = name, y = mle, ymin = cilow, ymax = cihigh, color = name)) +
  geom_pointrange(position=position_dodge(width = 0.5)) +
  geom_hline(yintercept = 1) +
  theme_bw(base_size = 18) +
  labs(y = "Global dN/dS", color = "Variant Type")
plot(p)

###################
## Plot all genes
###################

df = dnds_all_sel_cv[1:25,]
df$gene_name = factor(df$gene_name, levels = unique(df$gene_name))
p<-ggplot(df, aes(x=gene_name, y=-log10(qglobal_cv)), width=0.75) +
  geom_bar(aes(y=-log10(pglobal_cv)), fill = "gray90", stat="identity", position = position_dodge()) + 
  geom_bar(stat="identity", position = position_dodge()) + 
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = 2) +
  labs(y="-log10(FDR)", x="") +
  guides(fill = F) +
  coord_flip(ylim = c(0,10)) +
  scale_fill_manual(values = c("P" = "#2fb3ca","S" = "#CA932F", "R" = "#CA2FB4"), drop=F) +
  theme_bw(base_size = 16) +
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.ticks = element_line(size = 0))
plot(p)

###################
## Plot dNdS-CV by subtype - by fraction
###################

p<-dnds_fraction_global %>%
  filter(name == 'wall') %>%
  mutate(subtype = factor(subtype, levels = c("IDHmut-codel","IDHmut-noncodel","IDHwt")),
         fraction = factor(fraction, levels = c("P", "S", "R")),
         fraction = fct_recode(fraction, "Primary only" = "P", "Shared" = "S", "Recurrence only" = "R")) %>%
  ggplot(aes(x = subtype, y = mle, ymin = cilow, ymax = cihigh, color = fraction)) +
  geom_pointrange(position=position_dodge(width = 0.5)) +
  scale_color_manual(values = c("Primary only" = "#CA2FB4", "Shared" = "#CA932F", "Recurrence only" = "#2fb3ca"), drop=F) +
  #facet_wrap( ~ subtype) +
  geom_hline(yintercept = 1) +
  theme_bw(base_size = 10) +
  coord_cartesian(ylim=c(0.75,2)) +
  labs(x = "Subtype", y = "Global dN/dS", color = "Variant Fraction")

pdf('~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/EDF7/dnds-fraction-hm.pdf', height = 4, width = 6)
plot(p)
dev.off()

###################
## Plot dNdS-CV by subtype - by sample
###################

p<-dnds_sample_global %>%
  filter(name == 'wall') %>%
  mutate(subtype = factor(subtype, levels = c("IDHmut-codel","IDHmut-noncodel","IDHwt")),
         sample_type = factor(sample_type, levels = c("P", "R")),
         sample_type = fct_recode(sample_type, "Primary" = "P", "Recurrence" = "R")) %>%
  ggplot(aes(x = subtype, y = mle, ymin = cilow, ymax = cihigh, color = sample_type)) +
  geom_pointrange(position=position_dodge(width = 0.5)) +
  scale_color_manual(values = c("Primary" = "#a6611a", "Recurrence" = "#018571"), drop=F) +
  #facet_wrap( ~ subtype) +
  geom_hline(yintercept = 1) +
  theme_bw(base_size = 10) +
  coord_cartesian(ylim=c(0.75,2)) +
  labs(x = "Subtype", y = "Global dN/dS", color = "Sample Type")

pdf('~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/EDF7/dnds-sample-hm.pdf', height = 4, width = 6)
plot(p)
dev.off()

###################
## Plot dNdS-CV genes by subtype - fraction
###################

dat <- dnds_fraction_sel_cv %>%
  mutate(fraction = factor(fraction, levels = c("P", "S", "R"))) %>%
  complete(gene_name,fraction,subtype, fill = list(qglobal_cv=1)) %>%
  group_by(fraction,subtype) %>%
  arrange(qglobal_cv, pglobal_cv) %>%
  ungroup() %>%
  arrange(fraction)

myplots <- lapply(split(dat, paste(dat$subtype, dat$fraction))[c(1,3,2,4,6,5,7,9,8)], function(df){
  df$gene_name = factor(df$gene_name, levels = unique(df$gene_name))
  p <- ggplot(df[1:7,], aes(x=gene_name, y=-log10(qglobal_cv), fill=fraction), width=0.75) +
    geom_bar(aes(y=-log10(pglobal_cv)), fill = "gray90", stat="identity", position = position_dodge()) + 
    geom_bar(stat="identity", position = position_dodge()) + 
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = 2) +
    labs(y="-log10(FDR)", x="") +
    guides(fill = F) +
    facet_wrap(~ subtype, scales = "free_y") + 
    coord_flip(ylim = c(0,10)) +
    scale_fill_manual(values = c("P" = "#2fb3ca","S" = "#CA932F", "R" = "#CA2FB4"), drop=F) +
    theme_bw(base_size = 16) +
    theme( panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           axis.ticks = element_line(size = 0))
  
  return(p)
})

pdf('~/The Jackson Laboratory/GLASS - Documents/Resubmission/Figures/EDF7/dnds-genes-hm.pdf', height = 8, width = 12)
do.call(grid.arrange,(c(myplots, ncol=3)))
dev.off()

## END ##
