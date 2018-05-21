#######################################################
# Generate a metadata json file for Life-History sequencing samples.
# Date: 2018.05.21
# Authors: Kevin, Samir.
#######################################################

library(tidyverse)
library(jsonlite)
library(purrr)
library(listviewer)

#######################################################
# An example of tabular data containing read groups, file locations, and basic subject covariate information.
metadf <- read_tsv("/Users/johnsk/Documents/Life-History/GLASS-WG/data/table_to_json_test.tsv")

# Some columns did not need to be included in JSON.
# This is a fabricated test that included a few BAM files 
# from these libraries to test parsing ability.
metajson <- metadf %>%
  select(-one_of("Mate_ID")) %>%
  group_by(Patient_ID, Cohort, Sex, Age) %>%
  rename(fileNames = Fastq_Filenames) %>%
  group_by(FlowCell_ID, Lane_ID, add = TRUE) %>%
  mutate(files = list(list(fileType, fileNames))) %>%    
  ungroup() %>%
  select(-one_of("fileType", "fileNames")) %>%
  filter(!duplicated(files)) %>%
  nest(-Patient_ID, -Cohort, -Sex, -Age, .key = PatientLevel) %>%
  mutate(PatientLevel = purrr::map(PatientLevel, ~ .x %>%
                                     group_by(Sample_Type) %>%
                                     nest(.key = SampleLevel)))
# Interactively assess list tree structure.
listviewer::jsonedit(metajson)

# If opening in SublimeText use app (shift + cmd + p) to get in JSON pretty format.
write_json(metajson, "/Users/johnsk/Documents/Life-History/GLASS-WG/data/metadata-json.json")


