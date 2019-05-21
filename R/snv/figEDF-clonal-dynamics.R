##############################################
# Clonal dynamics for selected samples.
# Updated: 2019.05.21
# Author: Kevin J.
##################################################

# Working directory for this analysis in the GLASS-analysis project. 
mybasedir = "/Volumes/verhaak-lab/GLASS-analysis/"
setwd(mybasedir)

#######################################################

# Necessary packages:
library(tidyverse)
library(DBI)
library(parallel)
library(clonevol)

#######################################################
# Establish connection with the GLASS database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB2")

# Case examples where we have multi-sector data as well as multiple samples.
tumor_pairs = dbReadTable(con,  Id(schema="analysis", table="tumor_pairs"))
surgeries = dbReadTable(con,  Id(schema="clinical", table="surgeries"))
gold_set = dbReadTable(con,  Id(schema="analysis", table="gold_set"))
# If interested in investigating the clonal residence of specific mutations.
# passanno = dbGetQuery(con, "SELECT * FROM variants.passanno")
variant_classifications = dbReadTable(con,  Id(schema="variants", table="variant_classifications"))

# Query to investigate the mutational overlap between multi-sector and multi-longitudinal samples.
tumor_mut_comparison_anno = dbGetQuery(con, "SELECT
tmc.tumor_pair_barcode,
tmc.case_barcode,
tmc.tumor_barcode_a,
tmc.tumor_barcode_b,
idh_codel_subtype,
received_alk,
hypermutator_status,
0 AS time_birth,
ca.case_age_diagnosis_years AS time_initial,
ROUND(ca.case_age_diagnosis_years + (tmc.surgical_interval_mo / 12.0),2) AS time_recurrence,
0 AS mf_birth,
mf1.coverage_adj_mut_freq AS mf_initial,
mf2.coverage_adj_mut_freq AS mf_recurrence,
tmc.count_a,
tmc.count_b,
tmc.union_ab,
tmc.intersection_ab,
tmc.setdiff_a,
tmc.setdiff_b,
mf1.cumulative_coverage AS cov_a,
mf2.cumulative_coverage AS cov_b,
LEAST(mf1.cumulative_coverage, mf2.cumulative_coverage) AS min_cov,
ROUND(setdiff_a::decimal / mf1.cumulative_coverage * 1e6, 4) AS mf_private_a,
ROUND(setdiff_b::decimal / mf2.cumulative_coverage * 1e6, 4) AS mf_private_b,
ROUND(intersection_ab::decimal / LEAST(mf1.cumulative_coverage, mf2.cumulative_coverage) * 1e6, 4) AS mf_shared
FROM analysis.tumor_mut_comparison tmc
INNER JOIN analysis.tumor_pairs stp ON tmc.tumor_pair_barcode = stp.tumor_pair_barcode
LEFT JOIN analysis.tumor_clinical_comparison ctp ON ctp.tumor_pair_barcode = stp.tumor_pair_barcode
LEFT JOIN analysis.mut_freq mf1 ON mf1.aliquot_barcode = tmc.tumor_barcode_a 
LEFT JOIN analysis.mut_freq mf2 ON mf2.aliquot_barcode = tmc.tumor_barcode_b 
LEFT JOIN clinical.subtypes su ON su.case_barcode = stp.case_barcode
LEFT JOIN clinical.cases ca ON ca.case_barcode = stp.case_barcode")

# Get the new clinical_tumor_pairs table.
clinical_tumor_pairs_query = read_file("sql/clinical/clinical-tumor-pairs-db2.sql")
clinical_tumor_pairs <- dbGetQuery(con, clinical_tumor_pairs_query)

#############################
##### Example 1: GLSS-MD-0027
#############################
# Load the loci associated with -TP-, -R1-, and -R2- for this sample.
loci <- read_tsv("results/pyclone/run/GLSS-MD-0027/tables/loci.tsv") 

# Offset clusters so that they are 1-based.
locidf <- loci %>% 
  mutate(sample_id = gsub("-", "_", substr(sample_id, 14, 18)),
         cluster = cluster_id + 1,
         variant_allele_frequency = variant_allele_frequency*100) %>%
  select(-cluster_id) %>%
  # Avoid alphabetical ordering.
  mutate(sample_id = gsub("TP", "S1", sample_id)) %>%
  mutate(sample_id = gsub("R1", "S2", sample_id)) %>%
  mutate(sample_id = gsub("R2", "S3", sample_id)) %>%
  gather(variable, value, -mutation_id, -sample_id, -cluster) %>%
  unite(temp, sample_id, variable) %>%
  spread(temp, value) %>%
  # Remove clusters with low number of mutations (9, 10, 11, 12).
  filter(!cluster%in%c(9, 10, 11, 12)) %>% 
  arrange(cluster) 

# Shorten vaf column names as they will be projected otherwise.
vaf.col.names <- grep('_variant_allele_frequency', colnames(locidf), value=T)
ccf.col.names <- grep('_cellular_prevalence$', colnames(locidf), value=T)
# Rename the samples.
sample.names <- gsub('_variant_allele_frequency', '', vaf.col.names)
locidf[, sample.names] <- locidf[, ccf.col.names]
locidf[, sprintf("%s_VAF", sample.names)] <- locidf[, vaf.col.names]
ccf.col.names <- sample.names
vaf.col.names <- sprintf("%s_VAF", sample.names)

# Retrieve the average CCF for a particular cluster:
sample_cluster_ccf = locidf %>% 
  select(-S1_01_variant_allele_frequency, -S2_01_variant_allele_frequency, -S3_01_variant_allele_frequency, -S1_01_cellular_prevalence_std,
         -S2_01_cellular_prevalence_std, -S3_01_cellular_prevalence_std) %>% 
  gather(sample, cellular_prevalence, c(S1_01_cellular_prevalence, S2_01_cellular_prevalence, S3_01_cellular_prevalence)) %>% 
  mutate(sample = gsub("_01_cellular_prevalence", "", sample),
         sample = recode(sample, "S1" = 0, "S2" = 55, "S3" = 76)) %>% 
  group_by(cluster, sample) %>% 
  summarise(median_ccf = median(cellular_prevalence)) %>% 
  ungroup()

# Count the number of mutations per cluster.
cluster_num = locidf %>% 
  group_by(cluster) %>% 
  summarise(mut_count = n())

# Join the number of mutations per cluster with the cluster name.
cluster_ccf_annot = sample_cluster_ccf %>% 
  left_join(cluster_num, by="cluster") %>% 
  select(cluster, mut_count) %>% 
  distinct()

# Create a visualization that displays the time with points indicating a surgical sample. Create gray areas for treatment.
pdf('/Users/johnsk/Documents/pyclone-example-MD-0027.pdf', width = 9, height = 6, useDingbats = FALSE, title='')
ggplot(data=sample_cluster_ccf, aes(x=sample, y=median_ccf, group=as.factor(cluster), color=as.factor(cluster))) + 
  geom_line(size=1.5) + geom_point(size=3) + theme_bw() + xlab("Time (months)") + ylab("Median Cluster CCF") +labs(color="Mutational clusters\n(mutation number)") +
  coord_cartesian(xlim=c(-5, 80), ylim = c(0, 1)) + annotate("text", x = c(0, 55, 76), y = c(0.775, 0.775, 0.775), label = c("Initial", "R1", "R2"), col = "red", size=8) +
  scale_color_discrete(labels= paste(cluster_ccf_annot$cluster, " (n=", cluster_ccf_annot$mut_count, ")", sep = "")) +
  annotate("rect", xmin = 57, xmax = 69, ymin = -5, ymax = 75, alpha = .2) + annotate("text", x = 63, y = 0.95, label = "TMZ\n12 cycles", col = "black") + ggtitle("GLSS-MD-0027")
dev.off()


#############################
##### Example 2: GLSS-AT-GP01
#############################

# A multi-sector IDHwt tumor with multiple samples from biopsy and recurrence.
loci_at <- read_tsv("results/pyclone/run/GLSS-AT-GP01/tables/loci.tsv") 

# Again, offset the mutational cluster labelling to be 1-based.
# There are no low-mutation clusters.
loci_at_df <- loci_at %>% 
  mutate(sample_id = gsub("-", "_", substr(sample_id, 14, 18)),
         cluster = cluster_id + 1,
         variant_allele_frequency = variant_allele_frequency*100) %>%
  select(-cluster_id) %>%
  mutate(sample_id = gsub("TP", "S1", sample_id)) %>%
  mutate(sample_id = gsub("R1", "S2", sample_id)) %>%
  gather(variable, value, -mutation_id, -sample_id, -cluster) %>%
  unite(temp, sample_id, variable) %>%
  spread(temp, value) %>%
  arrange(cluster) 

# Restructure some of the variable names.
vaf.col.names <- grep('_variant_allele_frequency', colnames(loci_at_df), value=T)
ccf.col.names <- grep('_cellular_prevalence$', colnames(loci_at_df), value=T)
sample.names <- gsub('_variant_allele_frequency', '', vaf.col.names)
loci_at_df[, sample.names] <- loci_at_df[, ccf.col.names]
loci_at_df[, sprintf("%s_VAF", sample.names)] <- loci_at_df[, vaf.col.names]
ccf.col.names <- sample.names
vaf.col.names <- sprintf("%s_VAF", sample.names)

# Retrieve the average CCF for a particular cluster:
sample_cluster_ccf = loci_at_df %>% 
  gather(sample, cellular_prevalence, c(S1_01_cellular_prevalence, S1_02_cellular_prevalence, S1_03_cellular_prevalence,
                                        S2_01_cellular_prevalence, S2_02_cellular_prevalence, S2_03_cellular_prevalence, 
                                        S2_04_cellular_prevalence)) %>% 
  mutate(sample = gsub("_cellular_prevalence", "", sample),
         sample = recode(sample, "S1_01" = -0.25, "S1_02" = 0, "S1_03" = 0.25, "S2_01" = 4, "S2_02" = 4.25, "S2_03" = 4.5, "S2_04" = 4.75)) %>% 
  group_by(cluster, sample) %>% 
  summarise(median_ccf = median(cellular_prevalence)) %>% 
  ungroup()

# Tabulate the number of mutations per PyClone cluster.
cluster_num = loci_at_df %>% 
  group_by(cluster) %>% 
  summarise(mut_count = n())
# Annotate them with the appropriate name.
cluster_ccf_annot = sample_cluster_ccf %>% 
  left_join(cluster_num, by="cluster") %>% 
  select(cluster, mut_count) %>% 
  distinct()

# It's challenging to create a time-based PyClone cluster for multisector samples.
# I tried to off-set them each by a few days to display the heterogeneity, but it's not helpful.
ggplot(data=sample_cluster_ccf, aes(x=sample, y=median_ccf, group=as.factor(cluster), color=as.factor(cluster))) + 
  geom_line(size=1.5) + geom_point(size=3) + theme_bw() + xlab("Time (months)") + ylab("Median Cluster CCF") +labs(color="Mutational clusters\n(mutation number)") +
scale_color_discrete(labels= paste(cluster_ccf_annot$cluster, " (n=", cluster_ccf_annot$mut_count, ")", sep = "")) 

# ALTERNATIVE: 
sample_cluster_ccf = loci_at_df %>% 
  gather(sample, cellular_prevalence, c(S1_01_cellular_prevalence, S1_02_cellular_prevalence, S1_03_cellular_prevalence,
                                        S2_01_cellular_prevalence, S2_02_cellular_prevalence, S2_03_cellular_prevalence, 
                                        S2_04_cellular_prevalence)) %>% 
  mutate(sample = gsub("_cellular_prevalence", "", sample)) %>% 
  group_by(cluster, sample) %>% 
  summarise(median_ccf = median(cellular_prevalence)) %>% 
  ungroup()
  
sample_cluster_ccf$time_point = ifelse(grepl("S1", sample_cluster_ccf$sample), "Initial", "Recurrence")

cluster_num = loci_at_df %>% 
  group_by(cluster) %>% 
  summarise(mut_count = n())


cluster_ccf_annot = sample_cluster_ccf %>% 
  left_join(cluster_num, by="cluster") %>% 
  select(cluster, mut_count) %>% 
  distinct()


# Get average CCF per cluster, and plot line thickness by mutation count within cluster. Join cluster 5 with passanno.
pdf('/Users/johnsk/Documents/pyclone-example-AT-GP01-multisector.pdf', width = 9, height = 6, useDingbats = FALSE, title='')
ggplot(data=sample_cluster_ccf, aes(x=sample, y=median_ccf, group=as.factor(cluster), color=as.factor(cluster))) + 
  geom_line(size=1.5) + geom_point(size=3) + theme_bw() + xlab("Multi-sector samples") + ylab("Median Cluster CCF") +labs(color="Mutational clusters\n(mutation number)") +
  scale_color_discrete(labels= paste(cluster_ccf_annot$cluster, " (n=", cluster_ccf_annot$mut_count, ")", sep = "")) + facet_wrap(~time_point, scales="free") + ggtitle("GLSS-AT-GP01")
dev.off()

# Now, collapse each of the multi-sectors into a representative time point sample.
sample_cluster_ccf_avg = sample_cluster_ccf %>% 
  group_by(time_point, cluster) %>% 
  summarise(avg_ccf = mean(median_ccf))

# This is illustrative of wiping out the heterogeneity in that exists in the multi-sector samples.
pdf('/Users/johnsk/Documents/pyclone-example-AT-GP01-timepoint.pdf', width = 9, height = 6, useDingbats = FALSE, title='')
ggplot(data=sample_cluster_ccf_avg, aes(x=time_point, y=avg_ccf, group=as.factor(cluster), color=as.factor(cluster))) + 
  geom_line(size=1.5) + geom_point(size=3) + theme_bw() + xlab("") + ylab("Multisector Average\n Cluster CCF") +labs(color="Mutational clusters\n(mutation number)") +
  scale_color_discrete(labels= paste(cluster_ccf_annot$cluster, " (n=", cluster_ccf_annot$mut_count, ")", sep = "")) + ggtitle("GLSS-AT-GP01")
dev.off()


#############################
##### Example 3: GLSS-SF-0004
#############################

# This is a IDH-mut noncodel, non-hypermutator.
loci_sf <- read_tsv("results/pyclone/run/GLSS-SF-0004/tables/loci.tsv") 

# Remove some of the small mutational clusters.
locidf_sf <- loci_sf %>% 
  mutate(sample_id = gsub("-", "_", substr(sample_id, 14, 18)),
         cluster = cluster_id + 1,
         variant_allele_frequency = variant_allele_frequency*100) %>%
  select(-cluster_id) %>%
  mutate(sample_id = gsub("TP", "S1", sample_id)) %>%
  mutate(sample_id = gsub("R1", "S2", sample_id)) %>%
  mutate(sample_id = gsub("R2", "S3", sample_id)) %>%
  mutate(sample_id = gsub("R3", "S4", sample_id)) %>%
  gather(variable, value, -mutation_id, -sample_id, -cluster) %>%
  unite(temp, sample_id, variable) %>%
  spread(temp, value) %>%
  filter(!cluster%in%c(10, 11, 12)) %>% 
  arrange(cluster) 


# ASIDE: If interested in determing in which cluster the IDH mutation resides.
idh1_passanno = passanno %>% 
  filter(gene_symbol == "IDH1") %>% 
  mutate_if(bit64::is.integer64, as.double)
# Combine the annotation with the loci for each sample.
sf_0004_loci = loci_sf %>% 
  inner_join(idh1_passanno, by=c("mutation_id"="variant_id"))


# Some data cleaning.
vaf.col.names <- grep('_variant_allele_frequency', colnames(locidf_sf), value=T)
ccf.col.names <- grep('_cellular_prevalence$', colnames(locidf_sf), value=T)
sample.names <- gsub('_variant_allele_frequency', '', vaf.col.names)
locidf_sf[, sample.names] <- locidf_sf[, ccf.col.names]
locidf_sf[, sprintf("%s_VAF", sample.names)] <- locidf_sf[, vaf.col.names]
ccf.col.names <- sample.names
vaf.col.names <- sprintf("%s_VAF", sample.names)


# Retrieve the average CCF for a particular cluster:
sample_cluster_ccf = locidf_sf %>% 
  select(-S1_01_variant_allele_frequency, -S2_01_variant_allele_frequency, -S3_01_variant_allele_frequency, -S4_01_variant_allele_frequency,
         -S1_01_cellular_prevalence_std, -S2_01_cellular_prevalence_std, -S3_01_cellular_prevalence_std, -S4_01_cellular_prevalence_std) %>% 
  gather(sample, cellular_prevalence, c(S1_01_cellular_prevalence, S2_01_cellular_prevalence, S3_01_cellular_prevalence, S4_01_cellular_prevalence)) %>% 
  mutate(sample = gsub("_01_cellular_prevalence", "", sample),
         sample = recode(sample, "S1" = 0, "S2" = 15, "S3" = 35, "S4" = 44)) %>% 
  group_by(cluster, sample) %>% 
  summarise(median_ccf = median(cellular_prevalence)) %>% 
  ungroup()

# Define the number of mutations per clusters.
cluster_num = locidf_sf %>% 
  group_by(cluster) %>% 
  summarise(mut_count = n())

# Annotate the clusters.
cluster_ccf_annot = sample_cluster_ccf %>% 
  left_join(cluster_num, by="cluster") %>% 
  select(cluster, mut_count) %>% 
  distinct()

# Provide a clinical and mutational cluster representation as example cases.
pdf('/Users/johnsk/Documents/pyclone-example-SF-0004.pdf', width = 9, height = 6, useDingbats = FALSE, title='')
ggplot(data=sample_cluster_ccf, aes(x=sample, y=median_ccf, group=as.factor(cluster), color=as.factor(cluster))) + 
  geom_line(size=1.5) + geom_point(size=3) + theme_bw() + xlab("Time (months)") + ylab("Median Cluster CCF") +labs(color="Mutational clusters\n(mutation number)") +
  coord_cartesian(xlim=c(-5, 50), ylim = c(0, 1)) + annotate("text", x = c(0, 15, 35, 44), y = c(0.775, 0.775, 0.775, 0.775), label = c("Initial", "R1", "R2", "R3"), col = "red", size=8) +
  scale_color_discrete(labels= paste(cluster_ccf_annot$cluster, " (n=", cluster_ccf_annot$mut_count, ")", sep = "")) + ggtitle("GLSS-SF-0004") +
  annotate("rect", xmin = 17, xmax = 24, ymin = -5, ymax = 75, alpha = .2) + annotate("rect", xmin = 36, xmax = 42, ymin = -5, ymax = 75, alpha = .2) + annotate("text", x = c(21, 39), y = 0.95, label = c("TMZ\n7 cycles", "TMZ\n6 cycles"), col = "black") 
dev.off()

