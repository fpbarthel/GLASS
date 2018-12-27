library(DBI)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

dat <- dbGetQuery(con, read_file("sql/mutation_freq_private_shared.sql"))
clindata <- dbGetQuery(con, "SELECT DISTINCT case_barcode, idh_codel_subtype FROM clinical.surgeries WHERE idh_codel_subtype IS NOT NULL")
dat <- dat %>% left_join(clindata) %>% filter(mf_a < 10, mf_b < 10)

ggplot(dat, aes(x = surgical_interval_mo, y = mf_shared)) + 
  geom_point() + 
  geom_smooth(method="lm") +
  facet_wrap(~idh_codel_subtype, scales = "free")

ggplot(dat, aes(x = surgical_interval_mo, y = mf_private_b)) + 
  geom_point() + 
  geom_smooth(method="lm") +
  facet_wrap(~idh_codel_subtype, scales = "free")

ggplot(dat, aes(x = surgical_interval_mo, y = mf_private_a)) + 
  geom_point() + 
  geom_smooth(method="lm") +
  facet_wrap(~idh_codel_subtype, scales = "free")