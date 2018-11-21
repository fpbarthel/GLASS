#######################################################
# Compare the copy number calls using different parameters.
# Date: 2018.11.20
# Author: Kevin J.
#######################################################

# Different copy number call sets:
mybasedir = '/Volumes/verhaak-lab/GLASS-analysis/results'
cnresultsdir = '/cnv/plotcr'

#######################################################
# Necessary packages:
library(tidyverse)
library(GenomicRanges)
library(openxlsx)
library(DBI)
library(EnvStats)

#######################################################
# Establish connection with the database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

# Retrieve the case_sources and biospecimen_aliquots from the Database.
case_sources = dbReadTable(con,  Id(schema="clinical",table="case_sources"))
aliquots = dbReadTable(con,  Id(schema="biospecimen",table="aliquots"))

# The directory may change for each iteration of copy number calls.
files = list.files("/Volumes/verhaak-lab/GLASS-analysis/results/cnv/plotcr", full.names = T, pattern = 'standardizedMAD.txt', recursive = T)

## Each .txt file only contains the standardized MAD.
mad_dat = mclapply(files, function(f){
  dat = tryCatch(read.delim(f,as.is=T, header=F), error=function(e) e)
  if(inherits(dat,'error')) {
    message(f, '\n', dat, '\n')
    return()
  }
  # Truncate the file name to just the aliquot_id.
  dat = dat %>%
    mutate(aliquot_id = gsub(".standardizedMAD.txt", "", basename(f)))
  
  return(dat)
  
}, mc.cores=20)

## Combine all the samples from the GLASS cohort.
glass_mad = data.table::rbindlist(mad_dat)

# Revise any cohort groups by combining tissue source sites from TCGA.
colnames(glass_mad)[1] = "plotcr_standardized_mad"
glass_mad_meta = glass_mad %>% 
  mutate(project_id= substr(aliquot_id, 1, 4),
         subject_id = substring(aliquot_id, 1, 12), 
         source_id = substring(aliquot_id, 6, 7),
         seq_type = substring(aliquot_id, 21, 23),
         group_id = paste0(project_id, "-", source_id)) %>% 
         mutate(group_id = recode(group_id, "TCGA-06"="TCGA-GBM", "TCGA-14" = "TCGA-GBM", "TCGA-19" = "TCGA-GBM" ,"TCGA-DH" = "TCGA-LGG", "TCGA-DU"= "TCGA-LGG", "TCGA-FG" = "TCGA-LGG",
                                  "TCGA-TM" = "TCGA-LGG", "TCGA-TQ" = "TCGA-LGG")) %>% 
         mutate(group_id = paste0(group_id, "-", seq_type))

# Identify whether there were any missing samples in the aliquots table. There should be the same number as aliquots.
anti_join(aliquots, glass_mad_meta, by=c("aliquot_barcode"="aliquot_id"))

# Annotate the copy number data with cohort identifiers.
glass_mad_annot = glass_mad_meta %>% 
  inner_join(aliquots, by=c("aliquot_id"="aliquot_barcode")) %>% 
  filter(aliquot_id!="GLSS-MD-LP07-R1-01D-WXS-D9JSYR") %>% 
  mutate(sample_type = substr(aliquot_id, 14, 15),
         tumor_normal = ifelse(sample_type=="NB" | sample_type=="NM", "normal", "tumor"))

# Separate the data into WXS and WGS to adjust the scales.
glass_mad_annot_wxs = glass_mad_annot %>% 
  filter(grepl("-WXS-", aliquot_id))
glass_mad_annot_wgs = glass_mad_annot %>% 
  filter(grepl("-WGS-", aliquot_id))

# Group classification for PON as defined by the databases for WXS data.
ggplot(glass_mad_annot_wxs, aes(x=as.factor(aliquot_batch), y=plotcr_standardized_mad)) + geom_boxplot() + geom_jitter(aes(color=aliquot_batch), width = 0.2) + theme_bw() +
  ggtitle("glass-wxs-pon") + ylab("Standardized MAD") + xlab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(color="aliquot_batch") + 
   facet_grid(~tumor_normal) + stat_n_text()

# Group classification for PON as defined by the database for WGS data.
ggplot(glass_mad_annot_wgs, aes(x=as.factor(aliquot_batch), y=plotcr_standardized_mad)) + geom_boxplot() + geom_jitter(aes(color=aliquot_batch), width = 0.2) + theme_bw() +
  ggtitle("glass-wgs-pon") + ylab("Standardized MAD") + xlab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(color="aliquot_batch") + 
  facet_grid(~tumor_normal) + stat_n_text()

##########################
####  Delta MAD    ###### 
##########################
# The directory may change for each iteration of copy number calls.
files = list.files("/Volumes/verhaak-lab/GLASS-analysis/results/cnv/plotcr", full.names = T, pattern = 'deltaMAD.txt', recursive = T)

## Each .txt file only contains the deltaMAD (standardized MAD - denoised MAD).
mad_dat = mclapply(files, function(f){
  dat = tryCatch(read.delim(f,as.is=T, header=F), error=function(e) e)
  if(inherits(dat,'error')) {
    message(f, '\n', dat, '\n')
    return()
  }
  # Truncate the file name to just the aliquot_id.
  dat = dat %>%
    mutate(aliquot_id = gsub(".deltaMAD.txt", "", basename(f)))
  
  return(dat)
  
}, mc.cores=20)

## Combine all the samples from the GLASS cohort.
glass_mad = data.table::rbindlist(mad_dat)

# Revise any cohort groups by combining tissue source sites from TCGA.
colnames(glass_mad)[1] = "plotcr_delta_mad"
glass_mad_meta = glass_mad %>% 
  mutate(project_id= substr(aliquot_id, 1, 4),
         subject_id = substring(aliquot_id, 1, 12), 
         source_id = substring(aliquot_id, 6, 7),
         seq_type = substring(aliquot_id, 21, 23),
         group_id = paste0(project_id, "-", source_id)) %>% 
  mutate(group_id = recode(group_id, "TCGA-06"="TCGA-GBM", "TCGA-14" = "TCGA-GBM", "TCGA-19" = "TCGA-GBM" ,"TCGA-DH" = "TCGA-LGG", "TCGA-DU"= "TCGA-LGG", "TCGA-FG" = "TCGA-LGG",
                           "TCGA-TM" = "TCGA-LGG", "TCGA-TQ" = "TCGA-LGG")) %>% 
  mutate(group_id = paste0(group_id, "-", seq_type))

# Identify whether there were any missing samples in the aliquots table. There should be the same number as aliquots.
anti_join(aliquots, glass_mad_meta, by=c("aliquot_barcode"="aliquot_id"))

# Annotate the copy number data with Hoon's cohort identifiers.
glass_mad_annot = glass_mad_meta %>% 
  inner_join(aliquots, by=c("aliquot_id"="aliquot_barcode")) %>% 
  filter(aliquot_id!="GLSS-MD-LP07-R1-01D-WXS-D9JSYR") %>% 
  mutate(sample_type = substr(aliquot_id, 14, 15),
         tumor_normal = ifelse(sample_type=="NB" | sample_type=="NM", "normal", "tumor"))

# Separate the data into WXS and WGS to adjust the scales.
glass_mad_annot_wxs = glass_mad_annot %>% 
  filter(grepl("-WXS-", aliquot_id))
glass_mad_annot_wgs = glass_mad_annot %>% 
  filter(grepl("-WGS-", aliquot_id))

# Group classification for PON as defined by the databases.
ggplot(glass_mad_annot_wxs, aes(x=as.factor(aliquot_batch), y=plotcr_delta_mad)) + geom_boxplot() + geom_jitter(aes(color=aliquot_batch), width = 0.2) + theme_bw() +
  ggtitle("glass-wxs-pon") + ylab("Delta MAD") + xlab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(color="aliquot_batch") + 
  facet_grid(~tumor_normal) + stat_n_text()

# Group classification for PON as defined by the database.
ggplot(glass_mad_annot_wgs, aes(x=as.factor(aliquot_batch), y=plotcr_delta_mad)) + geom_boxplot() + geom_jitter(aes(color=aliquot_batch), width = 0.2) + theme_bw() +
  ggtitle("glass-wgs-pon") + ylab("Delta MAD") + xlab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(color="aliquot_batch") + 
  facet_grid(~tumor_normal) + stat_n_text()

##########################
##  Scaled Delta MAD    ## 
##########################
# The directory may change for each iteration of copy number calls.
files = list.files("/Volumes/verhaak-lab/GLASS-analysis/results/cnv/plotcr", full.names = T, pattern = 'scaledDeltaMAD.txt', recursive = T)

## Each .txt file only contains the scaledDeltaMAD (standardized MAD - denoised MAD)/standardized MAD.
mad_dat = mclapply(files, function(f){
  dat = tryCatch(read.delim(f,as.is=T, header=F), error=function(e) e)
  if(inherits(dat,'error')) {
    message(f, '\n', dat, '\n')
    return()
  }
  # Truncate the file name to just the aliquot_id.
  dat = dat %>%
    mutate(aliquot_id = gsub(".scaledDeltaMAD.txt", "", basename(f)))
  
  return(dat)
  
}, mc.cores=20)

## Combine all the samples from the GLASS cohort.
glass_mad = data.table::rbindlist(mad_dat)

# Revise any cohort groups by combining tissue source sites from TCGA.
colnames(glass_mad)[1] = "plotcr_scaled_delta_mad"
glass_mad_meta = glass_mad %>% 
  mutate(project_id= substr(aliquot_id, 1, 4),
         subject_id = substring(aliquot_id, 1, 12), 
         source_id = substring(aliquot_id, 6, 7),
         seq_type = substring(aliquot_id, 21, 23),
         group_id = paste0(project_id, "-", source_id)) %>% 
  mutate(group_id = recode(group_id, "TCGA-06"="TCGA-GBM", "TCGA-14" = "TCGA-GBM", "TCGA-19" = "TCGA-GBM" ,"TCGA-DH" = "TCGA-LGG", "TCGA-DU"= "TCGA-LGG", "TCGA-FG" = "TCGA-LGG",
                           "TCGA-TM" = "TCGA-LGG", "TCGA-TQ" = "TCGA-LGG")) %>% 
  mutate(group_id = paste0(group_id, "-", seq_type))

# Identify whether there were any missing samples in the aliquots table. There should be the same number as aliquots.
failed_samples = anti_join(aliquots, glass_mad_meta, by=c("aliquot_barcode"="aliquot_id"))

# Annotate the copy number data with Hoon's cohort identifiers.
glass_mad_annot = glass_mad_meta %>% 
  inner_join(aliquots, by=c("aliquot_id"="aliquot_barcode")) %>% 
  filter(aliquot_id!="GLSS-MD-LP07-R1-01D-WXS-D9JSYR") %>% 
  mutate(sample_type = substr(aliquot_id, 14, 15),
         tumor_normal = ifelse(sample_type=="NB" | sample_type=="NM", "normal", "tumor"))

# Separate the data into WXS and WGS to adjust the scales.
glass_mad_annot_wxs = glass_mad_annot %>% 
  filter(grepl("-WXS-", aliquot_id))
glass_mad_annot_wgs = glass_mad_annot %>% 
  filter(grepl("-WGS-", aliquot_id))

# Group classification for PON as defined by the databases.
ggplot(glass_mad_annot_wxs, aes(x=as.factor(aliquot_batch), y=plotcr_scaled_delta_mad)) + geom_boxplot() + geom_jitter(aes(color=aliquot_batch), width = 0.2) + theme_bw() +
  ggtitle("glass-wxs-pon") + ylab("Scaled Delta MAD") + xlab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(color="aliquot_batch") + 
  facet_grid(~tumor_normal) + stat_n_text()

# Group classification for PON as defined by the database.
ggplot(glass_mad_annot_wgs, aes(x=as.factor(aliquot_batch), y=plotcr_scaled_delta_mad)) + geom_boxplot() + geom_jitter(aes(color=aliquot_batch), width = 0.2) + theme_bw() +
  ggtitle("glass-wgs-pon") + ylab("Scaled Delta MAD") + xlab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(color="aliquot_batch") + 
  facet_grid(~tumor_normal) + stat_n_text()

# Load in WGS metrics

