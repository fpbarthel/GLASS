#######################################################
# Analyse gene drivers over time.
# Date: 2019.01.16 
# Author: Kevin J.
#######################################################

# Necessary packages:
library(tidyverse)
library(DBI)
library(gridExtra)
library(nlme)
library(EnvStats)

#######################################################

# Establish connection with Floris' database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

# Load additional tables from the database.
pairs = dbReadTable(con,  Id(schema="analysis",table="pairs"))
tumor_pairs = dbReadTable(con,  Id(schema="analysis", table="tumor_pairs"))
aliquots = dbReadTable(con,  Id(schema="biospecimen", table="aliquots"))
surgeries = dbReadTable(con,  Id(schema="clinical", table="surgeries"))
cases = dbReadTable(con,  Id(schema="clinical", table="cases"))
titan_param = dbReadTable(con,  Id(schema="analysis",table="titan_params"))
mut_freq = dbReadTable(con,  Id(schema="analysis",table="mutation_freq"))
neutrality_aliquots = dbReadTable(con,  Id(schema="analysis",table="neutrality_aliquots"))
neutrality_tumor_pairs = dbReadTable(con,  Id(schema="analysis",table="neutrality_tumor_pairs"))

# These tables **MAY** change, especially the driver table.
clinal_tumor_pairs = dbGetQuery(con,"SELECT * FROM analysis.clinical_by_tumor_pair")
drivers = dbGetQuery(con, "SELECT * FROM analysis.driver_status")
silver_set = dbGetQuery(con, "SELECT * FROM analysis.silver_set")

# Combine with pairs file to access "tumor_barcode".
titan_info = titan_param %>% 
  left_join(pairs, by="pair_barcode")

# Retrieve the correct subtype variable.
case_level_subtype = surgeries %>% 
  select(case_barcode, idh_codel_subtype) %>% 
  distinct() %>% 
  filter(!is.na(idh_codel_subtype)) 

# Be able to categorize a tumor_pair by its surgery number.
surgeries_sub = surgeries %>% 
  select(sample_barcode, surgery_number)

# Test whether there is an increase in the number of drivers per tumor type.
drivers_silver = drivers %>% 
  select(-tumor_barcode_a, -tumor_barcode_b, -case_barcode) %>% 
  inner_join(silver_set, by="tumor_pair_barcode") %>% 
  mutate(sample_barcode_a = substr(tumor_barcode_a, 1, 15),
         sample_barcode_b = substr(tumor_barcode_b, 1, 15)) %>% 
  left_join(surgeries_sub, by=c("sample_barcode_a" = "sample_barcode")) %>% 
  left_join(surgeries_sub, by=c("sample_barcode_b" = "sample_barcode")) %>% 
  mutate(surgery_pair = paste(surgery_number.x, surgery_number.y, sep="-")) %>% 
  left_join(case_level_subtype, by="case_barcode") %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, "IDHmut_codel" = "IDHmut codel", "IDHmut_noncodel" = "IDHmut noncodel", "IDHwt_noncodel" = "IDHwt"))


# There is a significant difference in the number of drivers in the tumor pair. This is probably driven by IDH being counted as driver (every codel/noncodel has +1)
kruskal.test(as.numeric(drivers_silver$driver_count), as.factor(drivers_silver$idh_codel_subtype))
ggplot(drivers_silver, aes(x=as.factor(idh_codel_subtype), y=as.numeric(driver_count))) + geom_boxplot() + theme_bw() +
  ylab("Total driver count") + xlab("")
  
# Examine whether there are any particular gains or losses in the non-shared IDHwt tumors.
no_shared_drivers = drivers_silver %>% 
  filter(is.na(driver_shared)) 
