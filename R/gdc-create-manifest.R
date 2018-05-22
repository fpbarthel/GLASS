## Create a manifest of TCGA whole genome (WGS) files from LGG and GBM cohorts
## Limit to primary-recurrent triplets and 2nd recurrences

library(GenomicDataCommons)
library(listviewer)
library(tidyverse)

setwd("/Volumes/Helix-Projects/GLASS-WG")

unpaired_json = "data/ref/TCGA_WGS_GDC_legacy_UUIDs.json"

## Make sure to include case ids
# "cases.project.project_id" = project (eg. TCGA-LGG)
# "cases.samples.sample_type" = sample type (eg. Primary Tumor)
# "cases.samples.portions.analytes.aliquots.submitter_id" = TCGA barcode

add_fields = c("cases.project.project_id",
               "cases.samples.sample_type", 
               "cases.samples.portions.analytes.aliquots.submitter_id")

## Inspect data types
files(legacy = TRUE) %>% facet(c('data_type')) %>% aggregations()

## Get a list of all WXS/RNA-Seq aligned BAM files from primary tumors
fq = files(legacy = TRUE) %>% 
  GenomicDataCommons::filter( ~ cases.project.project_id %in% c("TCGA-LGG", "TCGA-GBM") & 
                                experimental_strategy == "WGS" & 
                                data_type == "Aligned reads") %>% 
  GenomicDataCommons::select(c(default_fields(files()), add_fields))

message(sprintf("Found %s hits", GenomicDataCommons::count(fq)))

## Extract results
fres = results_all(fq)

## Inspect list using listviewer::jsonedit
jsonedit(fres)

## Flatten nested variables
fres$project = map(fres$cases, "project") %>% map_chr("project_id")
fres$sample_type = map(fres$cases, "samples") %>% map(unlist) %>% map_chr("sample_type")
fres$aliquot_id = map(fres$cases, "samples") %>% map(unlist) %>% map_chr("portions.analytes.aliquots.submitter_id")

## Convert to dataframe (finally!)
df = as.data.frame(fres[-which(names(fres) %in% c("cases", "acl", "analysis"))], stringsAsFactors = F) %>%
  select(id, aliquot_id, project, sample_type, experimental_strategy, file_size, md5sum, file_name, created_datetime, updated_datetime) %>%
  filter(grepl("TCGA", project)) %>%
  mutate(sample_id = substr(aliquot_id,1,16),
         case_id = substr(aliquot_id,1,12),
         format = "BAM")

## Pair primary-recurrent-2ndrecurrence samples
filtered_files = df %>% 
  group_by(sample_id) %>%
  mutate(p = order(file_size, decreasing = T)) %>%
  ungroup() %>%
  group_by(case_id) %>%
  mutate(hasRec = any(sample_type == "Recurrent Tumor")) %>%
  ungroup() %>%
  filter(hasRec, p == 1) %>%
  select(-hasRec, -p, -created_datetime, -updated_datetime, -experimental_strategy) %>%
  mutate(file_size_readable = gdata::humanReadable(file_size, standard="Unix"),
         sample_type_code = recode_factor(substr(aliquot_id, 14,16), "01A" = "TP", "01B" = "TP", "02A" = "R1", "02B" = "R2", "10A" = "NB", "10B" = "NB", "10D" = "NB"))

nested_filtered_files = filtered_files %>% ### filter(case_id == "TCGA-DU-6404") %>% ## CONSIDER INCLUDING FOR EASIER TEST
  nest(-case_id, -project, .key = samples) %>%
  mutate(samples = map(samples, nest, -sample_type_code, -sample_type, -sample_id, .key=files))

jsonlite::toJSON(nested_filtered_files, pretty = T)
write(jsonlite::toJSON(nested_filtered_files, pretty = T), file = unpaired_json)

# paired_files = filtered_files %>%
#   mutate(sample_type_numeric = recode_factor(substr(aliquot_id, 14,16), "01A" = "P", "01B" = "P", "02A" = "R1", "02B" = "R2", "10A" = "N", "10B" = "N", "10D" = "N")) %>%
#   select(case_id, project, sample_type_numeric, id) %>%
#   spread(sample_type_numeric, id)


  

