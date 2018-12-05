#######################################################
# Determination of samples that were mismatch per GATK fingerprinting tool
# Date: 2018.12.05
# Author: Kevin J.
#######################################################

# Current location of files:
workingdir = 'Volumes/verhaak-lab/GLASS-analysis/results/fingerprinting/'

#######################################################
# Necessary packages:
library(tidyverse)
library(DBI)
#######################################################
# Establish connection with the database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

# Retrieve the biospecimen_aliquots from the Database.
mut_freq = dbReadTable(con,  Id(schema="analysis",table="mutation_freq"))
aliquots = dbReadTable(con,  Id(schema="biospecimen",table="aliquots"))

# Read in the large file from fingerprinting results:
cross_check_results = read_tsv("GLASS.crosscheck_metrics", skip = 6)

# Identify those samples, within a case, that do not match each other:
crosscheck_sample_mismatch = cross_check_results %>% 
  mutate(caseA = substr(LEFT_GROUP_VALUE, 1, 12),
         caseB = substr(RIGHT_GROUP_VALUE, 1, 12),
         cases_match = caseA==caseB) %>% 
  filter(RESULT=="EXPECTED_MISMATCH") %>% 
  filter(cases_match=="TRUE")

# Are any of these cases hypermutators? ANSWER: No, not really. Also, nearly all the hypermutators did not display a mismatch.
crosscheck_mut_freq = mut_freq %>% 
  inner_join(crosscheck_sample_mismatch, by=c("aliquot_barcode"="LEFT_GROUP_VALUE"))

# Identify samples, from different cases, that match each other.
crosscheck_match = cross_check_results %>% 
  mutate(caseA = substr(LEFT_GROUP_VALUE, 1, 12),
         caseB = substr(RIGHT_GROUP_VALUE, 1, 12),
         cases_match = caseA==caseB) %>% 
  filter(RESULT=="UNEXPECTED_MATCH") %>% 
  filter(cases_match=="FALSE")

# Produce a list of the samples that either mismatch their case's other samples AND/OR match another subject in the cohort.
# 1. Samples that DO NOT match the other samples in their own case.
list_sample_mismatch = crosscheck_sample_mismatch %>% 
  group_by(RIGHT_GROUP_VALUE) %>% 
  summarise(freq = n()) %>% 
  filter(freq > 1) %>% 
  mutate(case_with_mismatch = substr(RIGHT_GROUP_VALUE, 1, 12), 
         sample_barcode = substr(RIGHT_GROUP_VALUE, 1, 15)) %>% 
  select(case_with_mismatch, mismatched_aliquot_barcode = RIGHT_GROUP_VALUE, mismatched_sample_barcode = sample_barcode)

# 2. Samples that DO match samples from ANOTHER case.
unexpected_matches = crosscheck_match %>%
  mutate(sampleA = substr(LEFT_GROUP_VALUE, 1, 15),
         sampleB = substr(RIGHT_GROUP_VALUE, 1, 15),
         sample_pair = paste(sampleA, sampleB, sep="_")) %>% 
  group_by(sampleA) %>% 
  summarise(freq = n()) %>% 
  filter(freq > 1)

# Extract the information for the mismatched_sample and where it is supposed to map.
list_unexpected_match = crosscheck_match %>%
  mutate(sampleA = substr(LEFT_GROUP_VALUE, 1, 15),
         sampleB = substr(RIGHT_GROUP_VALUE, 1, 15)) %>% 
  filter(sampleA%in%unexpected_matches$sampleA) %>% 
  select(sample_mismatched = sampleA, matching_case = caseB) %>% 
  distinct()

# Samples to be dropped (do not match any other sample in this cohort):
unmatched_samples = as.data.frame(list_sample_mismatch$mismatched_sample_barcode[!list_sample_mismatch$mismatched_sample_barcode%in%list_unexpected_match$sample_mismatched])
colnames(unmatched_samples) = "sample_barcode"

# Create a list of the unmatched aliquots.
unmatched_aliquots = unmatched_samples %>% 
  inner_join(aliquots, by="sample_barcode") %>% 
  filter(grepl("-WGS", aliquot_batch)) %>%  # Hack to only retain the WGS, the GLSS-MD-LP03-R1*-WXS sample was not a mismatch.
  select(aliquot_id:aliquot_batch, unmatched_sample_barcode = sample_barcode)

# Samples that can be renamed (they do match another case):
mismatched_aliquots = list_unexpected_match %>% 
  inner_join(aliquots, by=c("sample_mismatched"="sample_barcode")) %>% 
  select(aliquot_id:aliquot_batch, sample_mismatched, matching_case)

# Write out files to be able to post to github:
write.table(unmatched_aliquots, file = "/Users/johnsk/Documents/Life-History/unmatched_aliquots.txt", sep="\t", row.names = F, col.names = T, quote = F)
write.table(mismatched_aliquots, file = "/Users/johnsk/Documents/Life-History/mismatched_aliquots.txt", sep="\t", row.names = F, col.names = T, quote = F)

