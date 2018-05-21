library(tidyverse)
library(stringi)

df = data.frame(client  = c(rep("A",4), rep("B",4)),
                case    = rep(c("P","R"),4),
                file    = stri_rand_strings(8, 6, '[A-Z]'))

df %>%
  split(.$patient) %>%
  map(split, .$sample)
