
library(tidyverse)
library(DBI)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")
dat <- dbGetQuery(con, read_file("sql/mf_longitudinal_analysis.sql"))

dat <- dat %>% 
  gather(v, value, time_birth:mf_recurrence) %>%
  separate(v, c("var", "descriptor")) %>%
  spread(var, value)

p <- ggplot(data = dat, aes(x = time, y = mf, group = tumor_pair_barcode, color = descriptor)) +
  #geom_line() +
  stat_smooth(aes(group = 1), method = "lm") +
  #stat_summary(aes(group = 1), fun.y = mean, geom = "point",
  #             shape = 17, size = 3) +
  facet_wrap(~hypermutator_status ) + 
  coord_cartesian(ylim = c(0,10))

p
