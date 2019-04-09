##################################################
# Determine whether purity is associated with mutational frequency.
# Updated: 2019.03.28
# Author: Kevin J.
##################################################

# Working directory for this analysis in the GLASS-analysis project. 
mybasedir = "/Volumes/verhaak-lab/GLASS-analysis/"
setwd(mybasedir)

##################################################
# Necessary packages:
library(tidyverse)
library(DBI)

##################################################
# Establish connection with database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB2")

# Load in the necessary information regarding purity
mutation_freq = dbGetQuery(con, "SELECT * FROM analysis.mut_freq")
ssm2_count = dbGetQuery(con, "SELECT * FROM variants.ssm2_count")
n_distinct(ssm2_count$aliquot_barcode)
coverage = dbGetQuery(con, "SELECT * FROM analysis.coverage")
coverage_tumor = coverage %>% 
  filter(!grepl("-NB-", aliquot_barcode))
n_distinct(coverage_tumor$aliquot_barcode)


pairs = dbReadTable(con,  Id(schema="analysis", table="pairs"))
subtypes = dbReadTable(con,  Id(schema="clinical", table="subtypes"))
blocklist = dbReadTable(con,  Id(schema="analysis", table="blocklist"))
seqz_params = dbReadTable(con,  Id(schema="variants", table="seqz_params"))
titan_params = dbReadTable(con,  Id(schema="variants", table="titan_params"))

# Results from anti-join are dependent upon the order of the objects.
anti_join(seqz_params, titan_params, by="pair_barcode")
anti_join(titan_params, seqz_params, by="pair_barcode")

# Merge all samples together and remove any blocklisted samples.
seqz_purity = seqz_params %>% 
  inner_join(pairs, by= "pair_barcode") %>% 
  inner_join(blocklist, by=c("tumor_barcode"="aliquot_barcode")) %>% 
  inner_join(mutation_freq, by=c("tumor_barcode"="aliquot_barcode")) %>% 
  mutate(case_barcode = substr(pair_barcode, 1, 12),
         sample_barcode = substr(pair_barcode, 1, 15),
         seq_type = substr(pair_barcode, 27, 29)) %>% 
  left_join(subtypes, by="case_barcode")

# Correlation between tumor purity and mutational frequency in general.
cor.test(seqz_purity$cellularity, seqz_purity$coverage_adj_mut_freq, method="spearman")
# Visualize the relationship between tumor purity and mutational frequency.
ggplot(seqz_purity, aes(x=cellularity, y=log10(coverage_adj_mut_freq), color=seq_type)) + geom_point() + theme_bw() + 
  ylab("log10(mutation freq.)") + xlab("Sequenza - cellularity") + ggtitle("GLASS samples (n=665)") +
  annotate("text", x = 0.75, y=-1, label = "rho=0.19, P=2.8E-07") +  geom_smooth(method='lm')

# Determine the relationship between tumor purity and mutational freq. adjusted 
seqz_purity_primary = seqz_purity %>% 
  filter(grepl("-TP-", pair_barcode))
seqz_purity_primary_codel = seqz_purity_primary %>% 
  filter(idh_codel_subtype == "IDHmut-codel")
seqz_purity_primary_noncodel = seqz_purity_primary %>% 
  filter(idh_codel_subtype == "IDHmut-noncodel")
seqz_purity_primary_wt = seqz_purity_primary %>% 
  filter(idh_codel_subtype == "IDHwt")

cor.test(seqz_purity_primary_codel$cellularity, seqz_purity_primary_codel$coverage_adj_mut_freq, method="spearman")
cor.test(seqz_purity_primary_noncodel$cellularity, seqz_purity_primary_noncodel$coverage_adj_mut_freq, method="spearman")
cor.test(seqz_purity_primary_wt$cellularity, seqz_purity_primary_wt$coverage_adj_mut_freq, method="spearman")

seqz_purity_recurrences = seqz_purity %>% 
  filter(!grepl("-TP-", pair_barcode))
seqz_purity_recurrences_codel = seqz_purity_recurrences %>% 
  filter(idh_codel_subtype == "IDHmut-codel")
seqz_purity_recurrences_noncodel = seqz_purity_recurrences %>% 
  filter(idh_codel_subtype == "IDHmut-noncodel")
seqz_purity_recurrences_wt = seqz_purity_recurrences %>% 
  filter(idh_codel_subtype == "IDHwt")

cor.test(seqz_purity_recurrences_codel$cellularity, seqz_purity_recurrences_codel$coverage_adj_mut_freq, method="spearman")
cor.test(seqz_purity_recurrences_noncodel$cellularity, seqz_purity_recurrences_noncodel$coverage_adj_mut_freq, method="spearman")
cor.test(seqz_purity_recurrences_wt$cellularity, seqz_purity_recurrences_wt$coverage_adj_mut_freq, method="spearman")


########################
# TITAN purity
########################
titan_purity = titan_params %>% 
  inner_join(pairs, by= "pair_barcode") %>% 
  inner_join(blocklist, by=c("tumor_barcode"="aliquot_barcode")) %>% 
  inner_join(mutation_freq, by=c("tumor_barcode"="aliquot_barcode")) %>% 
  mutate(case_barcode = substr(pair_barcode, 1, 12),
         sample_barcode = substr(pair_barcode, 1, 15),
         seq_type = substr(pair_barcode, 27, 29)) %>% 
  left_join(subtypes, by="case_barcode")

# Correlation between tumor purity and mutational frequency in general.
cor.test(titan_purity$purity, titan_purity$coverage_adj_mut_freq, method="spearman")
# Visualize the relationship between tumor purity and mutational frequency.
ggplot(titan_purity, aes(x=purity, y=log10(coverage_adj_mut_freq), color=seq_type)) + geom_point() + theme_bw() + 
  ylab("log10(mutation freq.)") + xlab("TITAN - purity") + ggtitle("GLASS samples (n=656)") +
  annotate("text", x = 0.75, y=-1, label = "rho=-0.05, P=0.23") +  geom_smooth(method='lm')

# Determine the relationship between tumor purity and mutational freq. adjusted 
titan_purity_primary = titan_purity %>% 
  filter(grepl("-TP-", pair_barcode))
titan_purity_primary_codel = titan_purity_primary %>% 
  filter(idh_codel_subtype == "IDHmut-codel")
titan_purity_primary_noncodel = titan_purity_primary %>% 
  filter(idh_codel_subtype == "IDHmut-noncodel")
titan_purity_primary_wt = titan_purity_primary %>% 
  filter(idh_codel_subtype == "IDHwt")

cor.test(titan_purity_primary_codel$purity, titan_purity_primary_codel$coverage_adj_mut_freq, method="spearman")
cor.test(titan_purity_primary_noncodel$purity, titan_purity_primary_noncodel$coverage_adj_mut_freq, method="spearman")
cor.test(titan_purity_primary_wt$purity, titan_purity_primary_wt$coverage_adj_mut_freq, method="spearman")

titan_purity_recurrences = titan_purity %>% 
  filter(!grepl("-TP-", pair_barcode))
titan_purity_recurrences_codel = titan_purity_recurrences %>% 
  filter(idh_codel_subtype == "IDHmut-codel")
titan_purity_recurrences_noncodel = titan_purity_recurrences %>% 
  filter(idh_codel_subtype == "IDHmut-noncodel")
titan_purity_recurrences_wt = titan_purity_recurrences %>% 
  filter(idh_codel_subtype == "IDHwt")

cor.test(titan_purity_recurrences_codel$purity, titan_purity_recurrences_codel$coverage_adj_mut_freq, method="spearman")
cor.test(titan_purity_recurrences_noncodel$purity, titan_purity_recurrences_noncodel$coverage_adj_mut_freq, method="spearman")
cor.test(titan_purity_recurrences_wt$purity, titan_purity_recurrences_wt$coverage_adj_mut_freq, method="spearman")




# What's the breakdown per subtype
low_seqz_purity %>% 
  group_by(idh_codel_subtype) %>% 
  summarise(counts = n())

plot(low_seqz_purity$cellularity, low_seqz_purity$coverage_adj_mut_freq)


titan_low_purity = low_seqz_purity %>% 
  inner_join(titan_purity, by = "pair_barcode") %>%  
  select(pair_barcode, cellularity, purity, idh_codel_subtype, coverage_adj_mut_freq = coverage_adj_mut_freq.y)

plot(titan_low_purity$cellularity, titan_low_purity$purity)
plot(titan_low_purity$cellularity, titan_low_purity$coverage_adj_mut_freq)
plot(titan_low_purity$purity, titan_low_purity$coverage_adj_mut_freq)
hist(titan_low_purity$purity)

titan_purity = titan_params %>% 
  inner_join(pairs, by= "pair_barcode") %>% 
  inner_join(blocklist, by=c("tumor_barcode"="aliquot_barcode")) %>% 
  inner_join(mutation_freq, by=c("tumor_barcode"="aliquot_barcode")) %>% 
  mutate(case_barcode = substr(pair_barcode, 1, 12),
         sample_barcode = substr(pair_barcode, 1, 15)) 

cor.test(titan_purity$purity, log10(titan_purity$coverage_adj_mut_freq), method="pearson")
cor.test(titan_purity$purity, titan_purity$coverage_adj_mut_freq, method="spearman")

plot(titan_purity$purity, titan_purity$coverage_adj_mut_freq)

hist(titan_purity$purity)

#############################################
taylor_data = readWorkbook("/Users/johnsk/Documents/Life-History/titan-analyses/data/taylor-cancer-cell-2018.xlsx", sheet = 1, rowNames = F, colNames = TRUE)

# Extract only the iniformation for GBM and LGG. Recode some of their metadata to fit our own.
glioma_dat = taylor_data %>% 
  filter(Type%in%c("GBM", "LGG")) %>% 
  mutate(sample_type_num = substr(Sample, 14, 15), 
         sample_type = recode(sample_type_num, "01" = "TP", "02" = "R1"),
         sample_barcode = paste(substr(Sample, 1, 12), sample_type, sep ="-")) %>% 
  select(Sample, sample_type, sample_barcode, ABSOLUTE_purity = Purity)

# Examine how the AneuploidyScore data determined by ABSOLUTE match with our estimates of aneuploidy.
glass_taylor_seqz_purity = seqz_purity %>%
  inner_join(glioma_dat, by="sample_barcode")
cor.test(glass_taylor_seqz_purity$cellularity, glass_taylor_seqz_purity$ABSOLUTE_purity, method="s")
ggplot(glass_taylor_seqz_purity, aes(ABSOLUTE_purity, cellularity)) + geom_point() + theme_bw() + xlab("ABSOLUTE purity") +
  ylab("Sequenza purity") +  ggtitle("TCGA glioma samples (n=98)") + geom_abline(slope = 1) + ylim(0,1) + xlim(0, 1)


glass_taylor_titan_purity = titan_purity %>%
  inner_join(glioma_dat, by="sample_barcode")
cor.test(glass_taylor_titan_purity$purity, glass_taylor_titan_purity$ABSOLUTE_purity, method="s")
ggplot(glass_taylor_titan_purity, aes(ABSOLUTE_purity, purity)) + geom_point() + theme_bw() + xlab("ABSOLUTE purity") +
  ylab("TITAN purity") + ggtitle("TCGA glioma samples (n=99)") + geom_abline(slope = 1) + ylim(0,1) + xlim(0, 1)




