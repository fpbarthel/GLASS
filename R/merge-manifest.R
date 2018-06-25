## Merge manifest
## Author: Floris Barthel
## Date: Jun 24 2018

setwd("~/projects/GLASS-WG/")

library(tidyverse)

manifest_dir      = "data/manifest"
cases_prefix      = "cases"
samples_prefix    = "samples"
aliquots_prefix   = "aliquots"
readgroups_prefix = "readgroups"
files_prefix      = "files"
pairs_prefix      = "pairs"

cases = list.files(manifest_dir, pattern = sprintf("%s.tsv", cases_prefix), recursive = T, full.names = T) %>%
  map(read.delim, as.is=T) %>%
  reduce(bind_rows) %>% ## map_dfr maybe bettre
  distinct() ## NEED TO REMOVE DUPLICATE ROWS IN SOURCE FILES

samples = list.files(manifest_dir, pattern = sprintf("%s.tsv", samples_prefix), recursive = T, full.names = T) %>%
  map(read.delim, as.is=T) %>%
  reduce(bind_rows)

aliquots = list.files(manifest_dir, pattern = sprintf("%s.tsv", aliquots_prefix), recursive = T, full.names = T) %>%
  map(read.delim, as.is=T) %>%
  reduce(bind_rows)

readgroups = list.files(manifest_dir, pattern = sprintf("%s.tsv", readgroups_prefix), recursive = T, full.names = T) %>%
  map(read.delim, as.is=T) %>%
  reduce(bind_rows)

files = list.files(manifest_dir, pattern = sprintf("%s.tsv", files_prefix), recursive = T, full.names = T) %>%
  map(read.delim, as.is=T) %>%
  reduce(bind_rows)

pairs = list.files(manifest_dir, pattern = sprintf("%s.tsv", pairs_prefix), recursive = T, full.names = T) %>%
  map(read.delim, as.is=T) %>%
  reduce(bind_rows)

print(sprintf("Exporting manifest as json files for snakemake use."))
write(jsonlite::toJSON(aliquots, pretty = T), file = sprintf("%s/%s.json", manifest_dir, aliquots_prefix))
write(jsonlite::toJSON(files, pretty = T), file = sprintf("%s/%s.json", manifest_dir, files_prefix))
write(jsonlite::toJSON(cases, pretty = T), file = sprintf("%s/%s.json", manifest_dir, cases_prefix))
write(jsonlite::toJSON(pairs, pretty = T), file = sprintf("%s/%s.json", manifest_dir, pairs_prefix))
write(jsonlite::toJSON(readgroups, pretty = T), file = sprintf("%s/%s.json", manifest_dir, readgroups_prefix))
write(jsonlite::toJSON(samples, pretty = T), file = sprintf("%s/%s.json", manifest_dir, samples_prefix))

## END ##