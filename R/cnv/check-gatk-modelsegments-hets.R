files = list.files("results/cnv/modelsegments/", pattern = "hets.tsv|hets.normal.tsv", recursive = T, full.names = T)
df = data.frame(fn = files, case = substr(basename(files),1,12), analysis = substr(basename(files),21,23), size = file.size(files))

library(tidyverse)

df2 <- df %>% group_by(case, analysis) %>% mutate(var = var(size)) %>% ungroup()


thets <- read_tsv("results/cnv/modelsegments/TCGA-DU-7304-R1-12D-WGS-TNHHDG/TCGA-DU-7304-R1-12D-WGS-TNHHDG.hets.tsv", comment = "@", col_types = "ciiicc")
nhets <- read_tsv("results/cnv/modelsegments/TCGA-DU-7304-R1-12D-WGS-TNHHDG/TCGA-DU-7304-R1-12D-WGS-TNHHDG.hets.normal.tsv", comment = "@", col_types = "ciiicc")

nhets <- nhets %>%
  mutate(ct = REF_COUNT + ALT_COUNT, vaf = ALT_COUNT / ct)

thets <- thets %>%
  mutate(ct = REF_COUNT + ALT_COUNT, vaf = ALT_COUNT / ct)

plot(density(nhets$vaf))
plot(density(thets$vaf))
