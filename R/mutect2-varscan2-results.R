#######################################################
# Inspect Mutect2 and Varsca2 results, compare with prior results
# Date: 2018.08.07
# Author: Kevin J.
#######################################################
# Local directory for github repo.
# We are comparing filtered Mutect2 and Varscan2 calls (not yet maf files).
mybasedir = "~/mnt/scratchhelix/johnsk/GLASS_WG_floris/results/mutect2/m2filter/"
setwd(mybasedir) 

# Location of files.
m2filter_total_snv_file <- paste(mybasedir, "/mutect2_snv_totals.txt", sep="") 
m2filter_samplename_file <- paste(mybasedir, "/samplenames.txt", sep="") 

#######################################################
# Load necessary packages.
library(tidyverse)
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
# Create list of names with the ".filtered2.vep_filters.txt".
filenames <- list.files(pattern = "*_filters.txt")
list_names <-substr(filenames, 1, 22)

# Load in all of the files, and set list names. 
mutect2_tbl <- list.files(pattern = "*_filters.txt") %>% 
  map(~read_table(., col_names = "filters")) %>% 
  set_names(list_names)
# save(mutect2_tbl, file="/Users/johnsk/Documents/mutect2_tbl_filters.RData")

# For each variant per subject, we separate out the different filters.
mutect2_split_tmp <- mutect2_tbl %>% map(function(x) strsplit(x$filters, ";"))

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
