library(dndscv)
library(tidyverse)
library(DBI)
library(ggthemes)
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
  sel_cv <- cbind(status = qres_status, dnds_subset$sel_cv[dnds_subset$sel_cv$qglobal_cv<0.05, c("gene_name","qglobal_cv")])
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
  sel_cv <- cbind(subtype = qres_subtype, dnds_subset$sel_cv[dnds_subset$sel_cv$qglobal_cv<0.05, c("gene_name","qglobal_cv")])
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
    dnds_subset = dndscv(qres_subset, refdb = "hg19", outmats = FALSE, max_coding_muts_per_sample = 500, outp = 1)
    cbind(status = qres_status, subtype = qres_subtype, dnds_subset$globaldnds)
  })
  res <- data.table::rbindlist(result_list)
})
dnds_private_shared_idh_codel <- data.table::rbindlist(result_list)
dnds_private_shared_idh_codel

dnds_private_shared_idh_codel %>%
  filter(name %in% c('wmis','wtru','wall')) %>%
  mutate(name = fct_recode(name, "Missense" = "wmis", "Truncating" = "wtru", "All" = "wall")) %>%
  ggplot(aes(x = status, y = mle, ymin = cilow, ymax = cihigh)) +
  geom_pointrange(aes(color = name), position=position_dodge(width = 0.5)) +
  facet_wrap( ~ subtype) +
  geom_hline(yintercept = 1) +
  theme_bw() +
  labs(x = "Variant Fraction", y = "Maximum likelihood estimation", color = "Variant Classification")
  
save.image('results/dndscv/dndscv_ascii.RData', ascii = TRUE)   

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