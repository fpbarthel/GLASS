#######################################################
# Examine Titan's output for potential ploidy biases in WXS data.
# Date: 2018.11.26 
# Author: Kevin J.
#######################################################

# Floris noted that TITAN was providing unusual ploidy outputs for GLASS:
# https://github.com/gavinha/TitanCNA/issues/46

# Directory for GLASS analysis and titan files.
mybasedir = 'Volumes/verhaak-lab/GLASS-analysis/'
datadir   = 'results/cnv/titanfinal/params'
pattern   = '.params.txt$'

#######################################################

# Necessary packages:
library(tidyverse)
library(DBI)
library(parallel)
library(openxlsx)
library(EnvStats)

#######################################################
# Establish connection with the database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

# Retrieve the biospecimen_aliquots from the Database.
pairs = dbReadTable(con,  Id(schema="analysis",table="pairs"))
aliquots = dbReadTable(con,  Id(schema="biospecimen",table="aliquots"))
surgeries = dbReadTable(con,  Id(schema="clinical",table="surgeries"))

## Read in the "*.params.txt" files.
files = list.files(datadir, full.names = T, pattern = pattern, recursive=T)

## Each params file has issues.
titan_cn_dat = mclapply(files, function(f){
  dat = tryCatch(read.delim(f, as.is=T, header= F), error=function(e) e)
  if(inherits(dat,'error')) {
    message(f, '\n', dat, '\n')
    return()
  }
  # Truncate the file name to just the sample_id.
  dat = dat %>%
    mutate(pair_barcode = gsub(".params.txt", "", basename(f))) 
  return(dat)
  
}, mc.cores=20)

## Combine all the samples from the GLASS cohort.
titan_glass_cn = data.table::rbindlist(titan_cn_dat)
colnames(titan_glass_cn)[1] <- "variable"
colnames(titan_glass_cn)[2] <- "value"

# Analyze puirty and ploidy estimates by seq_type and cohort:
phi = titan_glass_cn %>% 
  filter(variable=="Average tumour ploidy estimate:") %>%  
  mutate(phi = round(as.numeric(value))) %>% 
  inner_join(pairs, by="pair_barcode") %>%
  inner_join(aliquots, by=c("tumor_barcode"="aliquot_barcode"))

# Frequency of Phi across cohorts: 2 (144); 3 (70); 4 (372)
table(phi$phi)

ploidy = titan_glass_cn %>% 
  filter(variable=="Average tumour ploidy estimate:") %>%   
  inner_join(pairs, by="pair_barcode") %>%
  inner_join(aliquots, by=c("tumor_barcode"="aliquot_barcode"))

purity = titan_glass_cn %>% 
  filter(variable=="Normal contamination estimate:")  %>%   
  inner_join(pairs, by="pair_barcode") %>%
  inner_join(aliquots, by=c("tumor_barcode"="aliquot_barcode"))

# Group classification for PON as defined by the databases for WXS data.
ggplot(ploidy, aes(x=as.factor(aliquot_batch), y=as.numeric(value))) + geom_boxplot() + geom_jitter(aes(color=aliquot_batch), width = 0.2) + theme_bw() +
  ggtitle("") + ylab("Titan ploidy") + xlab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(color="aliquot_batch") + 
  facet_grid(~aliquot_analysis_type) + stat_n_text()

ggplot(purity, aes(x=as.factor(aliquot_batch), y=as.numeric(value))) + geom_boxplot() + geom_jitter(aes(color=aliquot_batch), width = 0.2) + theme_bw() +
  ggtitle("") + ylab("Titan purity") + xlab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(color="aliquot_batch") + 
  facet_grid(~aliquot_analysis_type) + stat_n_text()

# Identify those samples that have both WGS and WXS.
# Subset to having both TITAN calls for WGS and WXS.
shared_wgs_wxs = ploidy %>% 
  mutate(sample_seq_barcode = paste(substr(tumor_barcode, 1, 15), aliquot_analysis_type, sep="-")) %>% 
  select(sample_seq_barcode, sample_barcode) %>% 
  distinct(sample_seq_barcode) %>% 
  mutate(sample_barcode = substr(sample_seq_barcode, 1, 15)) %>% 
  count(sample_barcode) %>% 
  filter(n > 1)
# Quick hack to get only those TCGA samples that have both WGS and WXS.
shared_wgs_wxs_tcga = shared_wgs_wxs %>% 
  filter(grepl("TCGA", sample_barcode)) %>% 
  mutate(subject_id = substr(sample_barcode, 1, 12)) 
unique(shared_wgs_wxs_tcga$subject_id)

# Ladder-plot for samples with both WGS and WXS.
polidy_wgs_wxs = ploidy %>% 
  filter(grepl('GLSS-HF-3016|GLSS-HK-|TCGA-06-0125|TCGA-06-0152|TCGA-06-0171|TCGA-06-0190|TCGA-06-0210|TCGA-06-0211|TCGA-06-0221|TCGA-14-1034|TCGA-14-1402|TCGA-19-1389|TCGA-DH-A669|TCGA-DU-5870|TCGA-DU-5872|TCGA-DU-6397|TCGA-DU-6404|TCGA-DU-6407|TCGA-DU-7304|TCGA-FG-5965|TCGA-FG-A4MT|TCGA-TM-A7CF|"TCGA-TQ-A7RK|TCGA-TQ-A7RV|TCGA-TQ-A8XE', tumor_barcode)) %>% 
  mutate(cohort = substr(aliquot_batch, 1, 7))
purity_wgs_wxs = purity %>% 
  filter(grepl('GLSS-HF-3016|GLSS-HK-|TCGA-06-0125|TCGA-06-0152|TCGA-06-0171|TCGA-06-0190|TCGA-06-0210|TCGA-06-0211|TCGA-06-0221|TCGA-14-1034|TCGA-14-1402|TCGA-19-1389|TCGA-DH-A669|TCGA-DU-5870|TCGA-DU-5872|TCGA-DU-6397|TCGA-DU-6404|TCGA-DU-6407|TCGA-DU-7304|TCGA-FG-5965|TCGA-FG-A4MT|TCGA-TM-A7CF|"TCGA-TQ-A7RK|TCGA-TQ-A7RV|TCGA-TQ-A8XE', tumor_barcode)) %>% 
  mutate(cohort = substr(aliquot_batch, 1, 7))

# Estimated ploidy when TITAN uses default selectionSolution.
ggplot(polidy_wgs_wxs, aes(x = aliquot_analysis_type, y = as.numeric(value), group=sample_barcode, color= cohort)) +
  geom_line(linetype="solid", size=1) + 
  geom_point(color="black", size=2) +  theme_bw() + ylab("TITAN estimated Ploidy") + xlab("") 

# Estimated purity when TITAN uses default selectionSolution.
ggplot(purity_wgs_wxs, aes(x = aliquot_analysis_type, y = as.numeric(value), group=sample_barcode, color= cohort)) +
  geom_line(linetype="solid", size=1) + ylim(0, 1) +
  geom_point(color="black", size=2) +  theme_bw() + ylab("TITAN estimated Purity") + xlab("")


##############################
# Retrieve all .params files:
##############################
# Read in parameters files for all ploidy 2 and ploidy 3 values.
## Read in the "*.params.txt" files.
param2_files = list.files("/Volumes/fastscratch/johnsk/titan/ploidy2", full.names = T, pattern = pattern, recursive=T)
param3_files = list.files("/Volumes/fastscratch/johnsk/titan/ploidy3", full.names = T, pattern = pattern, recursive=T)
# param4_files = list.files("/Volumes/fastscratch/johnsk/titan/ploidy4", full.names = T, pattern = pattern, recursive=T)


# Isolate just ploidy2.
titan_ploidy2 = mclapply(param2_files, function(f){
  dat = tryCatch(read.delim(f, as.is=T, header= F), error=function(e) e)
  if(inherits(dat,'error')) {
    message(f, '\n', dat, '\n')
    return()
  }
  # Truncate the file name to just the sample_id.
  dat = dat %>%
    mutate(pair_barcode = gsub(".params.txt", "", basename(f))) 
  return(dat)
  
}, mc.cores=20)

## Combine all the samples from the GLASS cohort.
titan_glass_ploidy2 = data.table::rbindlist(titan_ploidy2)
colnames(titan_glass_ploidy2)[1] <- "variable"
colnames(titan_glass_ploidy2)[2] <- "ploidy2"

ploidy2_results = titan_glass_ploidy2 %>% 
  filter(variable=="S_Dbw validity index (Both):") %>%   
  separate(pair_barcode, c("sample_barcode", "cluster"), sep="_", remove = F) %>% 
  group_by(sample_barcode) %>% 
  slice(which.min(ploidy2))

##### Isolate just ploidy3.
titan_ploidy3 = mclapply(param3_files, function(f){
  dat = tryCatch(read.delim(f, as.is=T, header= F), error=function(e) e)
  if(inherits(dat,'error')) {
    message(f, '\n', dat, '\n')
    return()
  }
  # Truncate the file name to just the sample_id.
  dat = dat %>%
    mutate(pair_barcode = gsub(".params.txt", "", basename(f))) 
  return(dat)
  
}, mc.cores=20)

## Combine all the samples from the GLASS cohort.
titan_glass_ploidy3 = data.table::rbindlist(titan_ploidy3)
colnames(titan_glass_ploidy3)[1] <- "variable"
colnames(titan_glass_ploidy3)[2] <- "ploidy3"

ploidy3_results = titan_glass_ploidy3 %>% 
  filter(variable=="S_Dbw validity index (Both):") %>%   
  separate(pair_barcode, c("sample_barcode", "cluster"), sep="_", remove = F) %>% 
  group_by(sample_barcode) %>% 
  slice(which.min(ploidy3))

#####
# Isolate just ploidy4.
#####
titan_ploidy4 = mclapply(param4_files, function(f){
  dat = tryCatch(read.delim(f, as.is=T, header= F), error=function(e) e)
  if(inherits(dat,'error')) {
    message(f, '\n', dat, '\n')
    return()
  }
  # Truncate the file name to just the sample_id.
  dat = dat %>%
    mutate(pair_barcode = gsub(".params.txt", "", basename(f))) 
  return(dat)
  
}, mc.cores=20)

## Combine all the samples from the GLASS cohort.
titan_glass_ploidy4 = data.table::rbindlist(titan_ploidy4)
colnames(titan_glass_ploidy4)[1] <- "variable"
colnames(titan_glass_ploidy4)[2] <- "ploidy4"

ploidy4_results = titan_glass_ploidy4 %>% 
  filter(variable=="S_Dbw validity index (Both):") %>%   
  separate(pair_barcode, c("sample_barcode", "cluster"), sep="_", remove = F) %>% 
  group_by(sample_barcode) %>% 
  slice(which.min(ploidy4))

# 1. Combined the ploidy2 results with the ploidy3 and ploidy 4. 
# 2. Find the minimum S_Dbw validity index across all ploidy.
# 3. What happens when you subset to just Ploidy2 and Ploidy 3?

#### Join data for all ploidy levels:
all_ploidy = ploidy2_results %>%
  inner_join(ploidy3_results, by="sample_barcode") %>% 
  inner_join(ploidy4_results, by="sample_barcode") %>% 
  select(sample_barcode, ploidy2, ploidy3, ploidy4) %>% 
  gather(ploidy, S_Dbw, ploidy2:ploidy4) %>%
  mutate_at(vars(starts_with("S_Dbw")), funs(as.numeric)) %>% 
  group_by(sample_barcode) %>% 
  slice(which.min(S_Dbw)) 

# Frequency of ploidy solution selected.
table(all_ploidy$ploidy)  # 2 (364), 3 (107); 4 (118)

####
colnames(ploidy2_results)[5] <- "ploidy2_cluster"
colnames(ploidy3_results)[5] <- "ploidy3_cluster"
ploidy2_ploidy3 = ploidy2_results %>%
  inner_join(ploidy3_results, by="sample_barcode") %>% 
  select(sample_barcode, ploidy2, ploidy3, ploidy2_cluster, ploidy3_cluster) %>% 
  gather(ploidy, S_Dbw, ploidy2:ploidy3) %>%
  mutate_at(vars(starts_with("S_Dbw")), funs(as.numeric)) %>% 
  group_by(sample_barcode) %>% 
  slice(which.min(S_Dbw)) 

# Frequency of ploidy solution selected.
table(ploidy2_ploidy3$ploidy)  # 2 (422), 3 (167)

tmp = ploidy2_ploidy3 %>% 
inner_join(pairs, by=c("sample_barcode"="pair_barcode")) %>%
  inner_join(aliquots, by=c("tumor_barcode"="aliquot_barcode"))


#### Combine with latest TITAN selectSolution.R for ploidy2/3 only:
## Read in the "*.params.txt" files.
files = list.files("/Volumes/verhaak-lab/GLASS-analysis/results/cnv/titanfinal_ploidy23_thres0.05_defaultalphaK/params/", full.names = T, pattern = pattern, recursive=T)

## Each params file has issues.
titan_cn_dat = mclapply(files, function(f){
  dat = tryCatch(read.delim(f, as.is=T, header= F), error=function(e) e)
  if(inherits(dat,'error')) {
    message(f, '\n', dat, '\n')
    return()
  }
  # Truncate the file name to just the sample_id.
  dat = dat %>%
    mutate(pair_barcode = gsub(".params.txt", "", basename(f))) 
  return(dat)
  
}, mc.cores=20)

## Combine all the samples from the GLASS cohort.
titan_thres05_selectSoln = data.table::rbindlist(titan_cn_dat)
colnames(titan_thres05_selectSoln)[1] <- "select_variable"
colnames(titan_thres05_selectSoln)[2] <- "select_value"

# Analyze puirty and ploidy estimates by seq_type and cohort:
SDbw_soln = ploidy2_ploidy3 %>% 
  mutate(SDbw_phi = gsub("ploidy", "", ploidy)) %>% 
  select(pair_barcode = sample_barcode, SDbw_ploidy2_clust = ploidy2_cluster, SDbw_ploidy3_clust = ploidy3_cluster,  SDbw_global_min = S_Dbw, SDbw_min_phi = SDbw_phi)

titan_phi_consensus = titan_thres05_selectSoln %>% 
  filter(select_variable=="Average tumour ploidy estimate:") %>%  
  mutate(phi = round(as.numeric(select_value))) %>% 
  inner_join(pairs, by="pair_barcode") %>%
  inner_join(aliquots, by=c("tumor_barcode"="aliquot_barcode")) %>% 
  select(pair_barcode, tumor_barcode, selectSoln_ploidy = select_value, selectSoln_phi = phi, aliquot_batch) %>% 
  inner_join(SDbw_soln, by="pair_barcode")
titan_phi_consensus$is_consistent <- titan_phi_consensus$selectSoln_phi==titan_phi_consensus$SDbw_min_phi
table(titan_phi_consensus$is_consistent)

#### Examine for threshold 0.2
files = list.files("/Volumes/verhaak-lab/GLASS-analysis/results/cnv/titanfinal/params/", full.names = T, pattern = pattern, recursive=T)

## Each params file has issues.
titan_cn_dat = mclapply(files, function(f){
  dat = tryCatch(read.delim(f, as.is=T, header= F), error=function(e) e)
  if(inherits(dat,'error')) {
    message(f, '\n', dat, '\n')
    return()
  }
  # Truncate the file name to just the sample_id.
  dat = dat %>%
    mutate(pair_barcode = gsub(".params.txt", "", basename(f))) 
  return(dat)
  
}, mc.cores=20)

## Combine all the samples from the GLASS cohort.
titan_thres1_selectSoln = data.table::rbindlist(titan_cn_dat)
colnames(titan_thres1_selectSoln)[1] <- "select_variable"
colnames(titan_thres1_selectSoln)[2] <- "select_value"

# Analyze puirty and ploidy estimates by seq_type and cohort:
titan_phi_thres1_consensus = titan_thres1_selectSoln %>% 
  filter(select_variable=="Average tumour ploidy estimate:") %>%  
  mutate(phi = round(as.numeric(select_value))) %>% 
  inner_join(pairs, by="pair_barcode") %>%
  inner_join(aliquots, by=c("tumor_barcode"="aliquot_barcode")) %>% 
  select(pair_barcode, tumor_barcode, selectSoln_ploidy = select_value, selectSoln_phi = phi, aliquot_batch) %>% 
  inner_join(SDbw_soln, by="pair_barcode")
titan_phi_thres1_consensus$is_consistent <- titan_phi_thres1_consensus$selectSoln_phi==titan_phi_thres1_consensus$SDbw_min_phi
table(titan_phi_thres1_consensus$is_consistent)
table(titan_phi_thres1_consensus$selectSoln_phi, titan_phi_thres1_consensus$aliquot_batch)
table(titan_phi_thres1_consensus$is_consistent, titan_phi_thres1_consensus$aliquot_batch)
table(titan_phi_thres1_consensus$SDbw_min_phi)
table(titan_phi_thres1_consensus$selectSoln_phi)

write.table(titan_phi_consensus, file="/Users/johnsk/Documents/titan_phi_consensus_thres0_05.txt", sep="\t", row.names = F, col.names = T, quote = F)
write.table(titan_phi_thres1_consensus, file="/Users/johnsk/Documents/titan_phi_consensus_thres0_1.txt", sep="\t", row.names = F, col.names = T, quote = F)

# Determine the minimum between `selectSoln` and `SDbw_min` for each sample.
  tmp = titan_phi_thres1_consensus %>%
  mutate_at(vars(starts_with("SDbw_min_phi")),funs(as.numeric)) %>% 
  select(selectSoln_phi, SDbw_min_phi) %>% 
  rowwise() %>% 
  mutate(min_overall_phi = min(selectSoln_phi, SDbw_min_phi))

### We can make a ggplot2 for paired samples:
polidy_wgs_wxs_compare = all_ploidy %>% 
  mutate(titan_ploidy = as.numeric(gsub("ploidy", "", ploidy))) %>% 
  filter(grepl('GLSS-HF-3016|GLSS-HK-|TCGA-06-0125|TCGA-06-0152|TCGA-06-0171|TCGA-06-0190|TCGA-06-0210|TCGA-06-0211|TCGA-06-0221|TCGA-14-1034|TCGA-14-1402|TCGA-19-1389|TCGA-DH-A669|TCGA-DU-5870|TCGA-DU-5872|TCGA-DU-6397|TCGA-DU-6404|TCGA-DU-6407|TCGA-DU-7304|TCGA-FG-5965|TCGA-FG-A4MT|TCGA-TM-A7CF|"TCGA-TQ-A7RK|TCGA-TQ-A7RV|TCGA-TQ-A8XE', sample_barcode)) %>% 
  mutate(sample_id = substr(sample_barcode, 1, 15),
         seq_type = substr(sample_barcode, 27, 29),
         cohort = substr(sample_barcode, 1, 7))

ggplot(polidy_wgs_wxs_compare, aes(x = seq_type, y = titan_ploidy, group=sample_id, alpha= 0.05)) +
  geom_line(linetype="solid", size=1, position=position_jitter(w=0.025, h=0)) + 
  geom_point(color="black", size=2) +  theme_bw() + ylab("TITAN Ploidy") + xlab("") 

## Only for ploidy2 and 3:
polidy2_3_wgs_wxs_compare = ploidy2_ploidy3 %>% 
  mutate(titan_ploidy = as.numeric(gsub("ploidy", "", ploidy))) %>% 
  filter(grepl('GLSS-HF-3016|GLSS-HK-|TCGA-06-0125|TCGA-06-0152|TCGA-06-0171|TCGA-06-0190|TCGA-06-0210|TCGA-06-0211|TCGA-06-0221|TCGA-14-1034|TCGA-14-1402|TCGA-19-1389|TCGA-DH-A669|TCGA-DU-5870|TCGA-DU-5872|TCGA-DU-6397|TCGA-DU-6404|TCGA-DU-6407|TCGA-DU-7304|TCGA-FG-5965|TCGA-FG-A4MT|TCGA-TM-A7CF|"TCGA-TQ-A7RK|TCGA-TQ-A7RV|TCGA-TQ-A8XE', sample_barcode)) %>% 
  mutate(sample_id = substr(sample_barcode, 1, 15),
         seq_type = substr(sample_barcode, 27, 29),
         cohort = substr(sample_barcode, 1, 7))

ggplot(polidy2_3_wgs_wxs_compare, aes(x = seq_type, y = titan_ploidy, group=sample_id, alpha= 0.05)) +
  geom_line(linetype="solid", size=1, position=position_jitter(w=0.025, h=0)) + 
  geom_point(color="black", size=2) +  theme_bw() + ylab("TITAN Ploidy") + xlab("") 


#########################
# Examine the log-likelihood for ploidy2 and ploidy3
#########################
ploidy2_results = titan_glass_ploidy2 %>% 
  filter(variable=="Log likelihood:") %>%   
  separate(pair_barcode, c("sample_barcode", "cluster"), sep="_", remove = F) %>% 
  group_by(sample_barcode) %>% 
  slice(which.min(ploidy2))


ploidy3_results = titan_glass_ploidy3 %>% 
  filter(variable=="Log likelihood:") %>%   
  separate(pair_barcode, c("sample_barcode", "cluster"), sep="_", remove = F) %>% 
  group_by(sample_barcode) %>% 
  slice(which.min(ploidy3))

ploidy2_ploidy3 = ploidy2_results %>%
  inner_join(ploidy3_results, by="sample_barcode") %>% 
  select(sample_barcode, ploidy2, ploidy3) %>% 
  gather(ploidy, log_like, ploidy2:ploidy3) %>%
  mutate_at(vars(starts_with("log_like")), funs(as.numeric)) %>% 
  group_by(sample_barcode) %>% 
  slice(which.min(log_like)) 
table(ploidy2_ploidy3$ploidy)

all_ploidy = ploidy2_results %>%
  inner_join(ploidy3_results, by="sample_barcode") %>% 
  inner_join(ploidy4_results, by="sample_barcode") %>% 
  select(sample_barcode, ploidy2, ploidy3, ploidy4) %>% 
  gather(ploidy, log_like, ploidy2:ploidy4) %>%
  mutate_at(vars(starts_with("log_like")), funs(as.numeric)) %>% 
  group_by(sample_barcode) %>% 
  slice(which.min(log_like)) 
table(all_ploidy$ploidy)

#############################
# Check WGD events from Meyerson paper
#############################
taylor_data = readWorkbook("/Users/johnsk/Documents/taylor-cancer-cell-2018.xlsx", sheet = 1, rowNames = F, colNames = TRUE)

glioma_dat = taylor_data %>% 
  filter(Type%in%c("GBM", "LGG")) %>% 
  mutate(sample_type_num = substr(Sample, 14, 15), 
         sample_type = recode(sample_type_num, "01" = "TP", "02" = "R1"),
         sample_barcode = paste(substr(Sample, 1, 12), sample_type, sep ="-"))
  
# From the initial ploidy variable defined above:
public_data_merge = ploidy %>%
  inner_join(glioma_dat, by="sample_barcode") %>% 
  select(tumor_barcode, titan_ploidy = value,  Sample, Type, `AneuploidyScore(AS)`, Genome_doublings, Purity) %>% 
  inner_join(purity, by="tumor_barcode")

public_data_merge %>% 
  group_by(Genome_doublings) %>% 
  summarize(Mean = mean(as.numeric(titan_ploidy), na.rm=TRUE))
  
#### Merge Meyerson data with S_Dbw_validity_index selected samples:
tmp = ploidy2_ploidy3 %>% 
  inner_join(pairs, by=c("sample_barcode"= "pair_barcode")) %>%
  inner_join(aliquots, by=c("tumor_barcode"="aliquot_barcode")) %>% 
  select(pair_barocde = sample_barcode.x, ploidy, S_Dbw, sample_barcode = sample_barcode.y) %>% 
  inner_join(glioma_dat, by="sample_barcode") %>% 
  select(pair_barocde, sample_barcode, Type, `AneuploidyScore(AS)`, titan_ploidy = ploidy, S_Dbw, Genome_doublings, Purity)

