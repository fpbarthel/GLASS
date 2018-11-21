#######################################################
# Identify the St. Jude pediatric brain tumors that have 
# availabe germline, diagnosis, and recurrent tumors.
# Date: 2018.05.14
# Author: Kevin J
#######################################################

# project directory.
setwd("/Users/johnsk/Documents/Life-History/GLASS-WG/")
StJude_dataset_path = "data/st-jude-data/StJude.20180511.xlsx"

#######################################################

library(tidyverse)
library(openxlsx)

#######################################################

# Roel provided the life history working group a file containing the St. Jude data set to which we were
# granted access.
StJude_avail_data = readWorkbook(StJude_dataset_path, sheet = 1, startRow = 1, colNames = TRUE)

################################
# The goal is to identify paired primary-recurrent samples that have WGS data.
################################
# What datasets are available?
table(StJude_avail_data$sj_diseases) 

# HGG = High Grade Glioma.
StJude_HGG = StJude_avail_data %>% 
  filter(sequencing_type=="WGS" & sj_diseases=="HGG" & file_type=="BAM")
StJude_HGG_bams = filter(StJude_HGG, !grepl('bai', file_path))
table(StJude_HGG_bams$subject_name, StJude_HGG_bams$sample_type)
# Looks like 7-pairs: 3 with sample at autopsy; 4 with samples at recurrence.

# LGG = Low Grade Glioma.
StJude_LGG = StJude_avail_data %>% 
  filter(sequencing_type=="WGS" & sj_diseases=="LGG" & file_type=="BAM")
StJude_LGG_bams = filter(StJude_LGG, !grepl('bai', file_path))
table(StJude_LGG_bams$subject_name, StJude_LGG_bams$sample_type)
# 3 LGGs with primary and relapse.

# Might there be patients that started off with LGG and progressed to HGG?
StJude_glioma = StJude_avail_data %>% 
  filter(sequencing_type=="WGS" & sj_diseases%in% c("LGG","HGG") & file_type=="BAM")
StJude_glioma_bams = filter(StJude_glioma, !grepl('bai', file_path))
StJude_glioma_table = table(StJude_glioma_bams$subject_name, StJude_glioma_bams$sample_type)

# Gather names on possible trios.
StJude_glioma_trios = StJude_glioma_bams %>% 
  group_by(subject_name) %>% summarise(Trio = n_distinct(sample_type)) %>% filter(Trio>2)
StJude_glioma_trio_names = StJude_glioma_trios$subject_name

# Check available data for all trio names. 
StJude_glioma_trio_data = StJude_avail_data[StJude_avail_data$subject_name%in%StJude_glioma_trio_names, ]
StJude_glioma_trio_data_wgs = StJude_glioma_trio_data %>% filter(sequencing_type=="WGS" & file_type=="BAM")

# Write out files to be downloaded from St. Jude Cloud.
write.csv(StJude_glioma_trio_data_wgs, "data/st-jude-data/st-jude_glioma_trio_data_wgs.csv")

