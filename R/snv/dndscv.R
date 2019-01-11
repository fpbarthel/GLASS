library(dndscv)
library(tidyverse)
library(DBI)
library(ggthemes)
library(gridExtra)
library(ggplot2)

###################
## DB connection
###################

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

res_fraction <- dbGetQuery(con, read_file("sql/dndscv_input_by_fraction.sql"))
res_sample <- dbGetQuery(con, read_file("sql/dndscv_input_by_sample.sql"))

###################
## Run dNdS for all samples in cohort
###################

dnds_all <- res_fraction %>%
  select(-subtype, -fraction, -evolution) %>% 
  distinct() %>%
  dndscv(refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = 500)

###################
## Examine results
###################
print(dnds_all$globaldnds)
print(dnds_all$nbreg$theta)

sel_cv = dnds_all$sel_cv
print(head(sel_cv, sum(sel_cv$qglobal_cv < 0.05)), digits = 3)

###################
## Compute dNdS for hotspot sites
###################
sitedn_all <- sitednds(dnds_all)

recursites = sitedn_all$recursites
print(head(recursites, 20), digits = 3)
print(recursites[recursites$qval < 0.05, c("gene","aachange","impact","freq","qval")])

###################
## Run dNdS seperately for neutralitytestr groups - by fraction
###################

result_list <- lapply(unique(res_fraction$subtype), function(subt) {
  result_list <- lapply(unique(res_fraction$evolution), function(ev) {
    result_list <- lapply(unique(res_fraction$fraction), function(fr) {
      message("Computing dNdS for ", subt, ", ", fr, " and ", ev)
      
      qres_subset = res_fraction %>% 
        filter(subtype == subt, fraction == fr, evolution == ev,
               complete.cases(evolution)) %>% 
        select(-evolution, -subtype, -fraction) %>% distinct()
      
      message("Found ", nrow(qres_subset), " variants in ", n_distinct(qres_subset$case_barcode), " patients.")
      
      tryCatch({
        dnds_subset <- dndscv(qres_subset, refdb = "hg19", outmats = FALSE, max_coding_muts_per_sample = 500, outp = 1)
      }, error = function(e) message(e))
      
      if(exists('dnds_subset')) {
        globaldnds <- cbind(subtype = subt, fraction = fr, evolution = ev, dnds_subset$globaldnds)
        #sel_cv <- cbind(subtype = st, fraction = fr, evolution = ev, dnds_subset$sel_cv[1:50,])
        list(globaldnds)#, sel_cv)
      }
      
    })
    res1 <- data.table::rbindlist(lapply(result_list,'[[',1))
    #res2 <- data.table::rbindlist(lapply(result_list,'[[',2))
    list(res1)#,res2)
  })
  res1 <- data.table::rbindlist(lapply(result_list,'[[',1))
  #res2 <- data.table::rbindlist(lapply(result_list,'[[',2))
  list(res1)#,res2)
})
dnds_neutralitytestr_fraction_global <- data.table::rbindlist(lapply(result_list,'[[',1))
#dnds_neutralitytestr_fraction_sel_cv <- data.table::rbindlist(lapply(result_list,'[[',2))
print(dnds_neutralitytestr_fraction_global)
#print(dnds_neutralitytestr_fraction_sel_cv)

###################
## Run dNdS seperately for neutralitytestr groups - by aliquot
###################

result_list <- lapply(unique(res_sample$subtype), function(subt) {
  result_list <- lapply(unique(res_sample$evolution), function(ev) {
    result_list <- lapply(unique(res_sample$sample_type), function(st) {
      message("Computing dNdS for ", subt, ", ", st, " and ", ev)
      
      qres_subset = res_sample %>% 
        filter(subtype == subt, sample_type == st, evolution == ev,
               complete.cases(evolution)) %>% 
        select(-evolution, -sample_type, -subtype) %>% 
        distinct()
      
      tryCatch({
        dnds_subset <- dndscv(qres_subset, refdb = "hg19", outmats = FALSE, max_coding_muts_per_sample = 500, outp = 1)
      }, error = function(e) message(e))
      
      if(exists('dnds_subset')) {
        globaldnds <- cbind(subtype = subt, sample_type = st, evolution = ev, dnds_subset$globaldnds)
        #sel_cv <- cbind(subtype = st, sample_type = st, evolution = ev, dnds_subset$sel_cv[1:50,])
        list(globaldnds)#, sel_cv)
      }
      
    })
    res1 <- data.table::rbindlist(lapply(result_list,'[[',1))
    #res2 <- data.table::rbindlist(lapply(result_list,'[[',2))
    list(res1)#,res2)
  })
  res1 <- data.table::rbindlist(lapply(result_list,'[[',1))
  #res2 <- data.table::rbindlist(lapply(result_list,'[[',2))
  list(res1)#,res2)
})
dnds_neutralitytestr_sample_global <- data.table::rbindlist(lapply(result_list,'[[',1))
#dnds_neutralitytestr_sample_sel_cv <- data.table::rbindlist(lapply(result_list,'[[',2))
print(dnds_neutralitytestr_sample_global)
#print(dnds_neutralitytestr_sample_sel_cv)

###################
## Run dNdS seperately for private/shared variants and IDH status
###################

result_list <- lapply(unique(res_fraction$subtype), function(st) {
  result_list <- lapply(unique(res_fraction$fraction), function(fr) {
    message("Computing dNdS for ", fr, " and ", st)
    qres_subset = res_fraction %>% filter(fraction == fr, subtype == st) %>% select(-subtype, -fraction, -evolution) %>% distinct()
    dnds_subset = dndscv(qres_subset, refdb = "hg19", outmats = FALSE, max_coding_muts_per_sample = 500)
    globaldnds <- cbind(fraction = fr, subtype = st, dnds_subset$globaldnds)
    sel_cv <- cbind(fraction = fr, subtype = st, dnds_subset$sel_cv[1:50,])
    list(globaldnds, sel_cv)
  })
  res1 <- data.table::rbindlist(lapply(result_list,'[[',1))
  res2 <- data.table::rbindlist(lapply(result_list,'[[',2))
  list(res1,res2)
})
dnds_fraction_global <- data.table::rbindlist(lapply(result_list,'[[',1))
dnds_fraction_sel_cv <- data.table::rbindlist(lapply(result_list,'[[',2))
print(dnds_fraction_global)
print(dnds_fraction_sel_cv)

###################
## Run dNdS seperately for all primaries and recurrences and IDH status
###################

result_list <- lapply(unique(res_sample$subtype), function(st) {
  result_list <- lapply(unique(res_sample$sample_type), function(sampt) {
    message("Computing dNdS for ", sampt, " and ", st)
    qres_subset = res_sample %>% filter(sample_type == sampt, subtype == st) %>% select(-subtype, -sample_type, -evolution) %>% distinct()
    dnds_subset = dndscv(qres_subset, refdb = "hg19", outmats = FALSE, max_coding_muts_per_sample = 500)
    globaldnds <- cbind(sample_type = sampt, subtype = st, dnds_subset$globaldnds)
    sel_cv <- cbind(sample_type = sampt, subtype = st, dnds_subset$sel_cv[1:50,])
    list(globaldnds, sel_cv)
  })
  res1 <- data.table::rbindlist(lapply(result_list,'[[',1))
  res2 <- data.table::rbindlist(lapply(result_list,'[[',2))
  list(res1,res2)
})
dnds_sample_global <- data.table::rbindlist(lapply(result_list,'[[',1))
dnds_sample_sel_cv <- data.table::rbindlist(lapply(result_list,'[[',2))
print(dnds_sample_global)
print(dnds_sample_sel_cv)

###################
## Cleanup
###################

rm(RefCDS, qres, qres_all, substmodel, sel_cv, result_list, gr_genes, covs)
rm(con)

save.image('results/dndscv/dndscv.RData')   
load('results/dndscv/dndscv.RData')

###################
## Plot dNdS-CV by neutralitytestr group and fraction
###################

dnds_neutralitytestr_fraction_global %>%
  filter(name %in% c('wmis','wtru','wall'), complete.cases(evolution)) %>%
  mutate(name = fct_recode(name, "Missense" = "wmis", "Truncating" = "wtru", "All" = "wall"),
         fraction = factor(fraction, levels = c("P", "S", "R"))) %>%
  ggplot(aes(x = fraction, y = mle, ymin = cilow, ymax = cihigh)) +
  geom_pointrange(aes(color = name), position=position_dodge(width = 0.5)) +
  facet_wrap( ~ evolution + subtype) +
  geom_hline(yintercept = 1) +
  coord_cartesian(ylim = c(-1,3)) +
  theme_bw() +
  labs(x = "Variant Fraction", y = "Global dN/dS", color = "Variant Classification")

###################
## Plot dNdS-CV by neutralitytestr group and sample
###################

dnds_neutralitytestr_sample_global %>%
  filter(name %in% c('wmis','wtru','wall')) %>%
  mutate(name = fct_recode(name, "Missense" = "wmis", "Truncating" = "wtru", "All" = "wall"),
         sample_type = factor(sample_type, levels = c("P", "R"))) %>%
  ggplot(aes(x = sample_type, y = mle, ymin = cilow, ymax = cihigh)) +
  geom_pointrange(aes(color = name), position=position_dodge(width = 0.5)) +
  facet_wrap( ~ evolution + subtype) +
  geom_hline(yintercept = 1) +
  theme_bw() +
  labs(x = "Variant Fraction", y = "Global dN/dS", color = "Variant Classification")


###################
## Plot dNdS-CV by subtype - by fraction
###################

dnds_fraction_global %>%
  filter(name %in% c('wmis','wtru','wall')) %>%
  mutate(subtype = fct_recode(subtype, "IDHmut-codel" = "IDHmut_codel", "IDHmut-noncodel" = "IDHmut_noncodel", "IDHwt" = "IDHwt_noncodel"),
         name = fct_recode(name, "Missense" = "wmis", "Truncating" = "wtru", "All" = "wall"),
         fraction = factor(fraction, levels = c("P", "S", "R"))) %>%
  ggplot(aes(x = name, y = mle, ymin = cilow, ymax = cihigh, color = fraction)) +
  geom_pointrange(position=position_dodge(width = 0.5)) +
  scale_color_manual(values = c("P" = "#2fb3ca","S" = "#CA932F", "R" = "#CA2FB4"), drop=F) +
  facet_wrap( ~ subtype) +
  geom_hline(yintercept = 1) +
  theme_bw(base_size = 18) +
  labs(x = "Variant Classification", y = "Global dN/dS", color = "Variant Fraction")

###################
## Plot dNdS-CV by subtype - by sample
###################

dnds_sample_global %>%
  filter(name %in% c('wmis','wtru','wall')) %>%
  mutate(subtype = fct_recode(subtype, "IDHmut-codel" = "IDHmut_codel", "IDHmut-noncodel" = "IDHmut_noncodel", "IDHwt" = "IDHwt_noncodel"),
         name = fct_recode(name, "Missense" = "wmis", "Truncating" = "wtru", "All" = "wall"),
         sample_type = factor(sample_type, levels = c("P", "R"))) %>%
  ggplot(aes(x = name, y = mle, ymin = cilow, ymax = cihigh, color = sample_type)) +
  geom_pointrange(position=position_dodge(width = 0.5)) +
  scale_color_manual(values = c("P" = "#2fb3ca", "R" = "#CA2FB4"), drop=F) +
  facet_wrap( ~ subtype) +
  geom_hline(yintercept = 1) +
  theme_bw(base_size = 18) +
  labs(x = "Variant Classification", y = "Global dN/dS", color = "Sample Type")

###################
## Plot dNdS-CV genes by subtype - fraction
###################

dat <- dnds_fraction_sel_cv %>%
  mutate(subtype = fct_recode(subtype, "IDHmut-codel" = "IDHmut_codel", "IDHmut-noncodel" = "IDHmut_noncodel", "IDHwt" = "IDHwt_noncodel"),
         fraction = factor(fraction, levels = c("P", "S", "R"))) %>%
  complete(gene_name,fraction,subtype, fill = list(qglobal_cv=1)) %>%
  group_by(fraction,subtype) %>%
  arrange(qglobal_cv) %>%
  ungroup() %>%
  arrange(fraction)

myplots <- lapply(split(dat, paste(dat$subtype, dat$fraction))[c(1,3,2,4,6,5,7,9,8)], function(df){
  df$gene_name = factor(df$gene_name, levels = unique(df$gene_name))
  p <- ggplot(df[1:7,], aes(x=gene_name, y=-log10(qglobal_cv), fill=fraction), width=0.75) +
    geom_bar(stat="identity", position = position_dodge()) + 
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = 2) +
    labs(x="-log10(FDR)", y="") +
    guides(fill = F) +
    facet_wrap(~ subtype, scales = "free_y") + 
    coord_flip(ylim = c(0,10)) +
    scale_fill_manual(values = c("P" = "#2fb3ca","S" = "#CA932F", "R" = "#CA2FB4"), drop=F) +
    theme_bw() 
  
  return(p)
})

do.call(grid.arrange,(c(myplots, ncol=3)))

###################
## Plot dNdS-CV genes by subtype - samples
###################

dat <- dnds_sample_sel_cv %>%
  mutate(subtype = fct_recode(subtype, "IDHmut-codel" = "IDHmut_codel", "IDHmut-noncodel" = "IDHmut_noncodel", "IDHwt" = "IDHwt_noncodel"),
         sample_type = factor(sample_type, levels = c("P", "R"))) %>%
  complete(gene_name,sample_type,subtype, fill = list(qglobal_cv=1)) %>%
  group_by(sample_type,subtype) %>%
  arrange(qglobal_cv) %>%
  ungroup()

myplots <- lapply(split(dat, paste(dat$subtype, dat$sample_type)), function(df){
  df$gene_name = factor(df$gene_name, levels = unique(df$gene_name))
  p <- ggplot(df[1:9,], aes(x=gene_name, y=-log10(qglobal_cv), fill=sample_type), width=0.75) +
    geom_bar(stat="identity", position = position_dodge()) + 
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = 2) +
    labs(x="-log10(FDR)", y="") +
    guides(fill = F) +
    facet_wrap(~ subtype, scales = "free_y") + 
    coord_flip(ylim = c(0,10)) +
    scale_fill_manual(values = c("P" = "#2fb3ca", "R" = "#CA2FB4"), drop=F) +
    theme_bw()
  
  return(p)
})

do.call(grid.arrange,(c(myplots, ncol=2)))

## END ##
