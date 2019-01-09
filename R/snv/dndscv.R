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
qres <- dbGetQuery(con, read_file("sql/dndscv_input.sql"))

###################
## Run dNdS for all samples in cohort
###################

qres_all <- qres %>% select(-idh_codel_subtype, -status) %>% distinct()
dnds_all <- dndscv(qres_all, refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = 500)
dnds_all_known_genes <- dndscv(qres_all, refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = 500, gene_list = known_cancergenes, maxcovs=5)

###################
## Examine results
###################
print(dnds_all$globaldnds)
print(dnds_all$nbreg$theta)

print(dnds_all_known_genes$globaldnds)
print(dnds_all_known_genes$nbreg$theta)

sel_cv = dnds_all$sel_cv
print(head(sel_cv, sum(sel_cv$qglobal_cv < 0.05)), digits = 3)
genes1 = sel_cv$gene_name[sel_cv$qglobal_cv < 0.05]

sel_cv = dnds_all_known_genes$sel_cv
print(head(sel_cv, sum(sel_cv$qglobal_cv < 0.05)), digits = 3)
genes2 = sel_cv$gene_name[sel_cv$qglobal_cv < 0.05]

###################
## Compute dNdS for hotspot sites
###################
sitedn_all <- sitednds(dnds_all)

recursites = sitedn_all$recursites
print(head(recursites, 20), digits = 3)
print(recursites[recursites$qval < 0.05, c("gene","aachange","impact","freq","qval")])

###################
## Run dNdS seperately for private/shared variants
###################

result_list <- lapply(unique(qres$status), function(qres_status) {
  message("Computing dNdS for ", qres_status)
  qres_subset = qres %>% filter(status == qres_status) %>% select(-idh_codel_subtype, -status) %>% distinct()
  dnds_subset <- dndscv(qres_subset, refdb = "hg19", outmats = FALSE, max_coding_muts_per_sample = 500)
  globaldnds <- cbind(status = qres_status, dnds_subset$globaldnds)
  sel_cv <- cbind(status = qres_status, dnds_subset$sel_cv[1:50, ])
  list(globaldnds, sel_cv)
})
dnds_private_shared_global <- data.table::rbindlist(lapply(result_list,'[[',1))
dnds_private_shared_sel_cv <- data.table::rbindlist(lapply(result_list,'[[',2))
print(dnds_private_shared_global)
print(dnds_private_shared_sel_cv)

###################
## Run dNdS seperately for each IDH-codel subtype
###################

result_list <- lapply(unique(qres$idh_codel_subtype), function(qres_subtype) {
  message("Computing dNdS for ", qres_subtype)
  qres_subset <- qres %>% filter(idh_codel_subtype == qres_subtype) %>% select(-idh_codel_subtype, -status) %>% distinct()
  dnds_subset <- dndscv(qres_subset, refdb = "hg19", outmats = FALSE, max_coding_muts_per_sample = 500)
  globaldnds <- cbind(subtype = qres_subtype, dnds_subset$globaldnds)
  sel_cv <- cbind(subtype = qres_subtype, dnds_subset$sel_cv[1:50, ])
  list(globaldnds, sel_cv)
})
dnds_idh_codel_global <- data.table::rbindlist(lapply(result_list,'[[',1))
dnds_idh_codel_sel_cv <- data.table::rbindlist(lapply(result_list,'[[',2))
print(dnds_idh_codel_global)
print(dnds_idh_codel_sel_cv)

###################
## Run dNdS seperately for private/shared variants and IDH status
###################

result_list <- lapply(unique(qres$idh_codel_subtype), function(qres_subtype) {
  result_list <- lapply(unique(qres$status), function(qres_status) {
    message("Computing dNdS for ", qres_status, " and ", qres_subtype)
    qres_subset = qres %>% filter(status == qres_status, idh_codel_subtype == qres_subtype) %>% select(-idh_codel_subtype, -status) %>% distinct()
    dnds_subset = dndscv(qres_subset, refdb = "hg19", outmats = FALSE, max_coding_muts_per_sample = 500)
    globaldnds <- cbind(status = qres_status, subtype = qres_subtype, dnds_subset$globaldnds)
    sel_cv <- cbind(status = qres_status, subtype = qres_subtype, dnds_subset$sel_cv[1:50,])
    list(globaldnds, sel_cv)
  })
  res1 <- data.table::rbindlist(lapply(result_list,'[[',1))
  res2 <- data.table::rbindlist(lapply(result_list,'[[',2))
  list(res1,res2)
})
dnds_private_shared_idh_codel_global <- data.table::rbindlist(lapply(result_list,'[[',1))
dnds_private_shared_idh_codel_sel_cv <- data.table::rbindlist(lapply(result_list,'[[',2))
print(dnds_private_shared_idh_codel_global)
print(dnds_private_shared_idh_codel_sel_cv)

rm(RefCDS, qres, qres_all, substmodel, sel_cv, result_list, gr_genes, covs)
rm(con)

save.image('results/dndscv/dndscv.RData')   
load('results/dndscv/dndscv.RData')

###################
## Plot dNdS-CV
###################

dnds_private_shared_idh_codel_global %>%
  filter(name %in% c('wmis','wtru','wall')) %>%
  mutate(name = fct_recode(name, "Missense" = "wmis", "Truncating" = "wtru", "All" = "wall"),
         status = factor(status, levels = c("P", "S", "R"))) %>%
  ggplot(aes(x = status, y = mle, ymin = cilow, ymax = cihigh)) +
  geom_pointrange(aes(color = name), position=position_dodge(width = 0.5)) +
  facet_wrap( ~ subtype) +
  geom_hline(yintercept = 1) +
  theme_bw() +
  labs(x = "Variant Fraction", y = "Global dN/dS", color = "Variant Classification")


###################
## Plot dNdS-CV genes by subtype
###################

dat <- dnds_private_shared_idh_codel_sel_cv %>%
  complete(gene_name,status,subtype, fill = list(qglobal_cv=1))

myplots <- lapply(split(dat, dat$subtype), function(df){
  genes <- dnds_idh_codel_sel_cv$gene_name[dnds_idh_codel_sel_cv$subtype == unique(df$subtype)][1:10]
  df <- df[df$gene_name %in% genes,]
  df$gene_name = factor(df$gene_name, levels = rev(genes))
  
  p <- ggplot(df, aes(x=gene_name, y=-log10(qglobal_cv), fill=status), width=0.75) +
    geom_bar(stat="identity", position = position_dodge()) + 
    geom_hline(yintercept = -log10(0.01), color = "red", linetype = 2) +
    facet_wrap(~ subtype, scales = "free_y") + 
    coord_flip() +
    theme_bw()
  
  return(p)
})

do.call(grid.arrange,(c(myplots, ncol=3)))


  ggplot(aes(x=gene_name, y=-log10(qglobal_cv), fill=status)) +
  geom_bar(stat="identity", position = position_dodge()) + 
  geom_hline(yintercept = -log10(0.01), color = "red", linetype = 2) +
  facet_wrap(~ subtype, scales = "free_y") + 
  coord_flip() +
  theme_bw()
           

###################################################################################################################################################################################################


#dbWriteTable(con, Id(schema="analysis",table="dndscv_global_by_subtype"), as.data.frame(res))

## Push to database
#dbWriteTable(con, Id(schema="analysis",table="dndscv_global"), dnds_all$globaldnds)
#dbWriteTable(con, Id(schema="analysis",table="dndscv_gene"), dnds_all$sel_cv)
#dbWriteTable(con, Id(schema="analysis",table="dndscv_sites"), sitedn_all$recursites)

## END ##

#dnds_wxs = dnds_all

# dnds_wgs = dnds_all
# 
# print(head(dnds_wxs$sel_cv, 20), digits = 3)
# print(head(dnds_wgs$sel_cv, 20), digits = 3)
# 
# print(dnds_wxs$globaldnds)
# print(dnds_wgs$globaldnds)
# 
# print(dnds_all$globaldnds)
# print(dnds_all$nbreg$theta)

## END ##