#######################################################
# Inspect Mutect2 and Varsca2 results, compare with prior TCGA results.
# Date: 2018.08.15
# Author: Kevin J.
#######################################################
# Local directory for github repo.
mybasedir = "~/mnt/scratchhelix/johnsk/GLASS_WG_floris/results/"
setwd(mybasedir) 

# We are comparing filtered Mutect2 and Varscan2 calls.
Mutect2dir =  "~/mnt/scratchhelix/johnsk/GLASS_WG_floris/results/mutect2/m2filter/"
VarScan2dir = "~/mnt/scratchhelix/johnsk/GLASS_WG_floris/results/varscan2/final/"

# Completed life-history barcodes.
life_history_barcodes = "/Users/johnsk/Documents/Life-History/GLASS-WG/data/sequencing-information/master_life_history_uniform_naming_complete.txt"

#######################################################

# Load necessary packages.
library(tidyverse)
library(maftools)

#######################################################

# Generate theme for ggplot:
my_theme <- function() {
  p <- theme(
    axis.text=element_text(size=rel(1.5), color='black'),
    axis.title.y=element_text(size=rel(1.5), margin=margin(0, 10, 0, 0)),
    axis.title.x=element_text(size=rel(1.5), margin=margin(10, 0, 0, 0)),
    axis.line = element_line(colour="black", size=1),
    axis.ticks.length = unit(.3, 'cm'),
    axis.ticks.margin = unit(.3, 'cm'),
    legend.position='right',
    legend.text=element_text(size=rel(1.75)),
    legend.title=element_text(size=rel(1.75), face='bold'),
    legend.key=element_rect(fill='transparent'),
    strip.text=element_text(size=rel(1.75)),
    panel.border=element_blank(),
    panel.grid.major=element_line(colour="grey60", size=0.1, linetype='solid'),
    panel.grid.minor=element_line(colour="grey60", size=0.1, linetype='dotted'),
    panel.background=element_rect(fill="transparent", colour = NA),
    plot.background=element_rect(fill="transparent", colour = NA)
  )
  return (p)
}

#######################################################
# Inspect Mutect2 filters applied to both SNVs and small indels.
setwd(Mutect2dir)

# Create list of names with the ".filtered2.vep_filters.txt". Sample "GLSS-MD-LP03" was removed due to poor quality.
filenames <- list.files(path = "results/mutect2/filters/", pattern = "*_filters.txt")
list_names <-substr(filenames, 1, 29)

# Load in all of the files, and set list names. Do once because it's slow.
# These files were made smaller by only presenting the filter column.

mutect2_tbl <- list.files(path = "results/mutect2/filters/", pattern = "*_filters.txt", full.names = TRUE) %>% 
  map(~read_table(., col_names = "filters")) 

mutect2_tbl <- mutect2_tbl %>% 
  set_names(list_names)

#       save(mutect2_tbl, file="/Users/johnsk/Documents/mutect2_tbl_filters.RData")
# Load "mutect2_tbl" file.
#load("/Users/johnsk/Documents/Life-History/mutect2_tbl_filters.RData")

# Each variant has M2 filters associated with its call. Some SNV/indels have many filters justifying why 
# a variant is tossed out. We are interested in the frequency of total filters and solo filters across all cohorts.
# Total filters applied to each sample.
total_filters <- mutect2_tbl %>% map(function(x) strsplit(x$filters, ";")) %>% map(unlist) %>% map(table)

# Create a df with each column being a filter and each row being a sample.
total_filter_df <- as.data.frame(do.call(dplyr::bind_rows, total_filters), stringsAsFactors = FALSE) 
total_filter_tidy <- total_filter_df %>% 
  mutate(sample_names = names(total_filters)) %>% 
  gather("total_filter", "total_count", artifact_in_normal:orientation_bias)
# Create a variable on which total and solo filter data sets can be merged.
total_filter_tidy$sample_filter <- paste(total_filter_tidy$sample_names, total_filter_tidy$total_filter, sep="-")

# We are also interested in the solo filter applied per sample. 
solo_filter_df <- tibble(sample_names = names(mutect2_tbl),  filters = mutect2_tbl) %>% 
  unnest() %>% 
  mutate(num_filter_applied = map(strsplit(filters, ";"), length), snv_id = row_number()) %>% 
  unnest() %>% 
  group_by(sample_names) %>% 
  filter(num_filter_applied == 1) %>%  
  count(solo_filters = filters) %>% 
  spread(solo_filters, n)

# Tidy the "solo_filter_df" so that it can be merged with the "total_filter_df".
solo_filter_tidy <- solo_filter_df %>% 
  gather("solo_filter", "solo_count", artifact_in_normal:t_lod)
solo_filter_tidy$sample_filter <- paste(solo_filter_tidy$sample_names, solo_filter_tidy$solo_filter, sep="-")

# We need to plot total/solo filters by cohort (boxplot).
all_filters_tidy <- inner_join(total_filter_tidy, solo_filter_tidy, by=c("sample_filter" = "sample_filter"))
all_filters_tidy$cohort <- sprintf("%s-%s", substr(all_filters_tidy$sample_names.x, 1,7), substr(all_filters_tidy$sample_names.x, 27,29))
#all_filters_tidy$cohort[grepl("HF", all_filters_tidy$sample_names.x)] <- "HF"
#all_filters_tidy$cohort[grepl("HK", all_filters_tidy$sample_names.x)] <- "HK"
#all_filters_tidy$cohort[grepl("GLSS-MD", all_filters_tidy$sample_names.x)] <- "MD-LP"
#all_filters_tidy$cohort[grepl("TCGA", all_filters_tidy$sample_names.x)] <- "TCGA"

# Create tables to provide metrics about M2 filters.
total_filters_table <- all_filters_tidy %>% 
  group_by(cohort, total_filter) %>% 
  summarize(Samples = n(),
            AvgTotalFilter = mean(total_count),
            TotalFilterStDev = sd(total_count),
            TotalFilterMin = min(total_count), 
            TotalFilterMax = max(total_count))
write.table(total_filters_table, file = sprintf("%s.%s", "/Users/johnsk/Documents/Life-History/GLASS-WG/data/sequencing-information/M2-filters/m2-total-filters", "tsv"), sep="\t", row.names = F, col.names = T, quote = F)
# Table for solo filters.
solo_filters_table <- all_filters_tidy %>% 
  group_by(cohort, solo_filter) %>% 
  summarize(Samples = n(),
            AvgTotalFilter = mean(solo_count),
            TotalFilterStDev = sd(solo_count),
            TotalFilterMin = min(solo_count), 
            TotalFilterMax = max(solo_count))
write.table(solo_filters_table, file = sprintf("%s.%s", "/Users/johnsk/Documents/Life-History/GLASS-WG/data/sequencing-information/M2-filters/m2-solo-filters", "tsv"), sep="\t", row.names = F, col.names = T, quote = F)

all_filters_tidy <- all_filters_tidy %>% group_by(cohort) %>% mutate(Samples = n_distinct(sample_names.x)) %>% ungroup() %>% mutate(cohort = sprintf("%s (n=%s samples)", cohort, Samples))

# Boxplots of total filters by cohort:
ggplot(all_filters_tidy, aes(x=total_filter, y=total_count, color=cohort))  + facet_wrap(~cohort, scales = "free_y") + 
  geom_boxplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("Total Filters Applied") +xlab("Filter Type") +
  ggtitle("Total M2 Filters") + guides(color=FALSE)

ggplot(all_filters_tidy, aes(x=total_filter, y=total_count, color=cohort))  + facet_grid(~cohort) + 
  geom_boxplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1))  + ylim(0,100000) +
    ylab("Total Filters Applied") +xlab("Filter Type") + ggtitle("Total M2 Filters Zoomed")

# Boxplots of solo filters applied:
ggplot(all_filters_tidy, aes(x=solo_filter, y=solo_count, color=cohort))+ facet_grid(~cohort) + 
  geom_boxplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("Solo Filters Applied") +xlab("Filter Type") +
  ggtitle("Solo M2 Filters")
ggplot(all_filters_tidy, aes(x=solo_filter, y=solo_count, color=cohort))+ facet_grid(~cohort) + 
  geom_boxplot(outlier.shape = NA) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0,10000) + 
  ylab("Solo Filters Applied") +xlab("Filter Type") + ggtitle("Solo M2 Filters Zoomed")

# Inspect filter frequencies for small indels too?

###############################
####    TCGA maf files      ###
###############################
# Link "list_names" object with original TCGA names using barcodes.
vcf_names <- substring(list_names, 1, 15)
life_history_barcode_sheet = read.delim(life_history_barcodes, as.is=T)
# Identify those samples that are used in the GLASS-WG project.
glass_tcga_samples <- life_history_barcode_sheet %>% 
  filter(Barcode%in%vcf_names) %>% 
  filter(grepl("TCGA", Original_ID)) %>% 
  select(Original_ID) %>% 
  .[["Original_ID"]]

# Load TCGA LGG-GBM mutation calls from exome data. Downloaded by Floris.
tcga_gbm_snv_calls <- read.table("/Users/johnsk/Documents/Life-History/TCGA_GBM.broad.mit.edu.n308.aggregated.maf.txt", header = T, sep ="\t", fill = TRUE)

# There seems to be relatively few variants per TCGA sample.
variant_totals <- tcga_gbm_snv_calls %>% 
  group_by(Tumor_Sample_Barcode) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n))
 
# Only three samples exist in the TCGA_GBM.broad.mit[...] data set.
glass_tcga_samples[glass_tcga_samples%in%as.character(tcga_snv_calls$Tumor_Sample_Barcode)]
# Only two mutations in one sample ("TCGA-06-0125-01A-01D-1490-08")?
tcga_gbm_snv_calls %>% 
  filter(Tumor_Sample_Barcode%in%glass_tcga_samples)

# What about the LGG data?
tcga_lgg_snv_calls = read.maf(maf = "/Users/johnsk/Documents/Life-History/PR_TCGA_LGG_PAIR_Capture_All_Pairs_QCPASS_v7.aggregated.capture.tcga.uuid.automated.somatic.maf")
# Example comparisons between TCGA-LGG and the data I am working with.
example_compare = tcgaCompare(maf = tcga_lgg_snv_calls, cohortName = "LOW_GRADE")
# These plots can take a long time to render.
oncoplot(maf = tcga_lgg_snv_calls, top = 10, fontSize = 12)
oncostrip(maf = tcga_lgg_snv_calls, genes = c("IDH1","TP53"))

# Load VarScan2 data.
setwd(VarScan2dir)
vs2_filenames <- list.files(pattern = "*.final.vcf$")
vs2_list_names <-substr(vs2_filenames, 1, 22)

# Load in all of the files, and set list names. 
varscan2_tbl <- list.files(pattern = "*.final.vcf$")[1:2] %>% 
map(~read.table(.)) %>% 
set_names(vs2_list_names[1:2])


# Load in the sample-specific barcodes: 



# Breaks down the per sample filter frequency. 
mutect2_HF_split <- map(mutect2_split_tmp[1:2], unlist)
HF_vec <- table(unlist(mutect2_HF_split))
HF_df <- data.frame(as.list(HF_vec))

mutect2_HK_split <- map(mutect2_split_tmp[3:12], unlist)
HK_vec <- table(unlist(mutect2_HK_split))
HK_df <- data.frame(as.list(HK_vec))

mutect2_MD_split <- map(mutect2_split_tmp[13:33], unlist)
MD_vec <- table(unlist(mutect2_MD_split))
MD_df <- data.frame(as.list(MD_vec))

mutect2_TCGA_split <- map(mutect2_split_tmp[34:83], unlist)
TCGA_vec <- table(unlist(mutect2_TCGA_split))
TCGA_df <- data.frame(as.list(TCGA_vec))

# Bring them all together.
filter_type_df <- bind_rows(TCGA_df, MD_df, HK_df, HF_df)
rownames(filter_type_df) <- c("TCGA", "MD", "HK", "HF")
filter_type_df_t <- as.data.frame(t(filter_type_df))
filter_type_df_t$filter_type <- rownames(filter_type_df_t) 
filter_type_df_tidy <- filter_type_df_t %>% gather(Cohort, Observations, TCGA:HF)


filter_type_df_mat <- t(filter_type_df)
filter_type_df_mat[is.na(filter_type_df_mat)] <- 0
filter_proportions <- as.data.frame(prop.table(filter_type_df_mat, 2))
filter_proportions$filter_type <- rownames(filter_proportions); rownames(filter_proportions) <- NULL
filter_type_df_tidy <- filter_proportions %>% gather(Cohort, Frequencies, TCGA:HF)

# Plot the relative frequencies with special theme.
special_theme <- function() {
  p <- theme(
    axis.text=element_text(size=rel(1.5), color='black'),
    axis.title.y=element_text(size=rel(1.5), margin=margin(0, 10, 0, 0)),
    axis.title.x=element_text(size=rel(1), margin=margin(10, 0, 0, 0)),
    axis.line = element_line(colour="black", size=1),
    axis.ticks.length = unit(.3, 'cm'),
    axis.ticks.margin = unit(.3, 'cm'),
    legend.position='right',
    legend.text=element_text(size=rel(1)),
    legend.title=element_text(size=rel(1), face='bold'),
    legend.key=element_rect(fill='transparent'),
    strip.text=element_text(size=rel(1)),
    panel.border=element_blank(),
    panel.grid.major=element_line(colour="grey60", size=0.1, linetype='solid'),
    panel.grid.minor=element_line(colour="grey60", size=0.1, linetype='dotted'),
    panel.background=element_rect(fill="transparent", colour = NA),
    plot.background=element_rect(fill="transparent", colour = NA)
  )
  return (p)
}


# Relevel the factor.
filter_type_df_tidy$Cohort <- factor(filter_type_df_tidy$Cohort, levels = c("MD", "HF", "HK", "TCGA"))
#  There are too many filters in the ggplot2 window.
ggplot(filter_type_df_tidy, aes(x=as.factor(Cohort), y=Frequencies, fill=filter_type)) + 
  geom_bar(stat = "identity")  + xlab("") + ylab("Relative Proportions") + special_theme() + guides(fill=guide_legend(ncol=4)) +
  theme(legend.position="bottom") + guides(fill=guide_legend(title="Filter Type"))

#######################################
#######     Mutect2         ###########
#######################################
# Total Mutect2 SNVs
total_snv_mutect2 <- read.table("mutect2_snv_totals.txt", sep="\t", header=F, stringsAsFactors = F)
total_snv_samplenames <-  read.table("total_snv_samplenames.txt", sep="\t", header=F, stringsAsFactors = F)
colnames(total_snv_mutect2) <- total_snv_samplenames$V1

# Total Mutect2 PASSED SNVs
passed_snv_mutect2 <- read.table("mutect2_snv_PASSED.txt", sep="\t", header=F, stringsAsFactors = F)
passed_snv_samplenames <-  read.table("PASSED_snv_samplenames.txt", sep="\t", header=F, stringsAsFactors = F)
colnames(passed_snv_mutect2) <- passed_snv_samplenames$V1

# Total Mutect2 PASSED indels
indel_mutect2 <- read.table("mutect2_indel_PASSED.txt", sep="\t", header=F, stringsAsFactors = F)
indel_samplenames <-  read.table("PASSED_indel_samplenames.txt", sep="\t", header=F, stringsAsFactors = F)
colnames(indel_mutect2) <- indel_samplenames$V1

# Look at the histograms for each of these values.
mutect2_all <- as.data.frame(t(bind_rows(total_snv_mutect2, passed_snv_mutect2, indel_mutect2)))
colnames(mutect2_all) <- c("total_snv", "passed_snv", "passed_indel")
mutect2_all$sample_names <- as.character(rownames(mutect2_all))
mutect2_all$Cohort <- paste(sapply(strsplit(mutect2_all$sample_names, "-"), "[[", 1), sapply(strsplit(mutect2_all$sample_names, "-"), "[[", 2), sep = "-")
mutect2_all$Cohort <- as.factor(gsub("TCGA-..", "TCGA", mutect2_all$Cohort))
mutect2_all$sample_type <- sapply(strsplit(mutect2_all$sample_names, "-"), "[[", 4)
mutect2_all$sample_type <- factor(mutect2_all$sample_type, levels = c("TP", "R1", "R2"))

# Plot total SNVs, passed SNVs, and indels by cohort.
ggplot(mutect2_all, aes(x = Cohort, y = total_snv, colour=Cohort)) + geom_boxplot() + ylab("Total SNVs") + 
  xlab("") + my_theme() + ggtitle("Mutect2 total SNVs") + ylim(0, 4E+06)
ggplot(mutect2_all, aes(x = Cohort, y = passed_snv, colour=Cohort)) + geom_boxplot() + ylab("Total SNVs") + 
  xlab("") + my_theme() + ggtitle("Mutect2 passed filter SNVs")
ggplot(mutect2_all, aes(x = Cohort, y = passed_indel, colour=Cohort)) + geom_boxplot() + ylab("Total SNVs") + 
  xlab("") + my_theme() + ggtitle("Mutect2 passed indels")

# Plot total SNVs, passed SNVs, and indels by primary-recurrent status by cohort. 
ggplot(mutect2_all, aes(x = sample_type, y = total_snv, colour=Cohort)) + geom_boxplot() + ylab("Total SNVs") + 
  xlab("") + my_theme() + ggtitle("Mutect2 total SNVs") + facet_wrap(~Cohort)
ggplot(mutect2_all, aes(x = sample_type, y = passed_snv, colour=Cohort)) + geom_boxplot() + ylab("Total SNVs") + 
  xlab("") + my_theme() + ggtitle("Mutect2 passed filter SNVs") + facet_wrap(~Cohort)  + ylim(0, 2E+04)
ggplot(mutect2_all, aes(x = sample_type, y = passed_indel, colour=Cohort)) + geom_boxplot() + ylab("Total SNVs") + 
  xlab("") + my_theme() + ggtitle("Mutect2 passed indels") + facet_wrap(~Cohort)

#############################
####### VarScan2 ############
#############################
# Compare filtered call in: results/varscan2/final/ with those prefiltered results in: results/varscan2/vcf/
# Set Varscan2 directory.
VarScan2dir = "~/mnt/scratchhelix/johnsk/GLASS_WG_floris/results/varscan2/final/"
setwd(VarScan2dir)

# Total varscan2 SNVs
total_snv_varscan2 <- read.table("varscan2_snv_totals.txt", sep="\t", header=F, stringsAsFactors = F)
total_snv_samplenames <-  read.table("total_snv_samplenames.txt", sep="\t", header=F, stringsAsFactors = F)
colnames(total_snv_varscan2) <- total_snv_samplenames$V1

# Total varscan2 PASSED SNVs
passed_snv_varscan2 <- read.table("varscan2_snv_PASSED.txt", sep="\t", header=F, stringsAsFactors = F)
passed_snv_samplenames <-  read.table("PASSED_snv_samplenames.txt", sep="\t", header=F, stringsAsFactors = F)
colnames(passed_snv_varscan2) <- passed_snv_samplenames$V1

# Total varscan2 PASSED indels
indel_varscan2 <- read.table("varscan2_indel_PASSED.txt", sep="\t", header=F, stringsAsFactors = F)
indel_samplenames <-  read.table("PASSED_indel_samplenames.txt", sep="\t", header=F, stringsAsFactors = F)
colnames(indel_varscan2) <- indel_samplenames$V1

# Look at the histograms for each of these values.
varscan2_all <- as.data.frame(t(bind_rows(total_snv_varscan2, passed_snv_varscan2, indel_varscan2)))
colnames(varscan2_all) <- c("total_snv", "passed_snv", "passed_indel")
varscan2_all$sample_names <- as.character(rownames(varscan2_all))
varscan2_all$Cohort <- paste(sapply(strsplit(varscan2_all$sample_names, "-"), "[[", 1), sapply(strsplit(varscan2_all$sample_names, "-"), "[[", 2), sep = "-")
varscan2_all$Cohort <- as.factor(gsub("TCGA-..", "TCGA", varscan2_all$Cohort))
varscan2_all$sample_type <- sapply(strsplit(varscan2_all$sample_names, "-"), "[[", 4)
varscan2_all$sample_type <- factor(varscan2_all$sample_type, levels = c("TP", "R1", "R2"))

# Plot total SNVs, passed SNVs, and indels by cohort.
ggplot(varscan2_all, aes(x = Cohort, y = total_snv, colour=Cohort)) + geom_boxplot() + ylab("Total SNVs") + 
  xlab("") + my_theme() + ggtitle("VarScan2 total SNVs")
ggplot(varscan2_all, aes(x = Cohort, y = passed_snv, colour=Cohort)) + geom_boxplot() + ylab("Total SNVs") + 
  xlab("") + my_theme() + ggtitle("VarScan2 passed filter SNVs")
ggplot(varscan2_all, aes(x = Cohort, y = passed_indel, colour=Cohort)) + geom_boxplot() + ylab("Total SNVs") + 
  xlab("") + my_theme() + ggtitle("VarScan2 passed indels") + ylim(0, 200000)

# Plot total SNVs, passed SNVs, and indels by primary-recurrent status by cohort. 
ggplot(varscan2_all, aes(x = sample_type, y = total_snv, colour=Cohort)) + geom_boxplot() + ylab("Total SNVs") + 
  xlab("") + my_theme() + ggtitle("VarScan2 total SNVs") + facet_wrap(~Cohort)
ggplot(varscan2_all, aes(x = sample_type, y = passed_snv, colour=Cohort)) + geom_boxplot() + ylab("Total SNVs") + 
  xlab("") + my_theme() + ggtitle("VarScan2 passed filter SNVs") + facet_wrap(~Cohort)  + ylim(0, 2E+04)
ggplot(varscan2_all, aes(x = sample_type, y = passed_indel, colour=Cohort)) + geom_boxplot() + ylab("Total SNVs") + 
  xlab("") + my_theme() + ggtitle("VarScan2 passed indels") + facet_wrap(~Cohort)


# Combine the calls from VarScan2 and Mutect2.
varscan2_all$caller <-"VarScan2"
mutect2_all$caller  <- "Mutect2"
all_calls <- bind_rows(varscan2_all, mutect2_all)

#############################
####### ladder plot #########
#############################
# Generate ladder plot for values.
ggplot(all_calls, aes(x = caller, y = passed_snv, group=sample_names, color= Cohort)) +
  geom_line(linetype="solid", size=1)+ 
  geom_point(color="black", size=2) + my_theme() + ylab("PASSED SNVs") + xlab("")

# Zoomed axis.
ggplot(all_calls, aes(x = caller, y = passed_snv, group=sample_names, color= Cohort)) +
  geom_line(linetype="solid", size=1) + 
  geom_point(color="black", size=2) + my_theme() + ylab("PASSED SNVs") + xlab("") + ylim(0, 15000)
