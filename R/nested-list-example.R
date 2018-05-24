library(tidyverse)
library(gapminder)

# json = gapminder %>% 
#   filter(continent == "Oceania") %>% ## Limit data to Oceania to get a smaller table
#   nest(-continent, .key = countries) %>%
#   mutate(countries = map(countries, nest, -country, .key=years))
# 
# jsonlite::toJSON(json, pretty = T)

library(tidyverse)
library(stringi)

n_patient = 2
n_samples = 3
n_readgroup = 4
n_mate = 2

df = data.frame(patient   = rep(rep(LETTERS[1:n_patient], n_samples),2),
                sample    = rep(rep(seq(1:n_samples), each = n_patient),2),
                readgroup = rep(stri_rand_strings(n_patient * n_samples * n_readgroup, 6, '[A-Z]'),2),
                mate      = rep(1:n_mate, each = n_patient * n_samples * n_readgroup)) %>%
  mutate(file = sprintf("%s.%s.%s_%s", patient, sample, readgroup, mate)) %>%
  arrange(file)

json = df %>% 
  nest(-patient, .key = samples) %>%
  mutate(samples = map(samples, nest, -sample, .key=readgroups))

json3 <- df %>% nest(-(1:3),.key=mate) %>% nest(-(1:2),.key=readgroups) %>% nest(-1,.key=samples)

jsonlite::toJSON(json3,pretty=T)

vars <- names(df)[-1] # or whatever variables you want to nest, order matters!


nest_by <- function(df, ..., reverse = T) {
  
  var_pairs <- map((length(vars)-1):1,~vars[.x:(.x+1)])
  json4 <- reduce(var_pairs,~{nm<-.y[1];nest(.x,.y,.key=!!enquo(nm))},.init=df)
}

test <- function(...) {
  enquo(...)
}

jsonlite::toJSON(json4,pretty=T)

json = df %>%
  group_by(patient) %>%
  group_by(sample, add = T) %>%
  nest()

jsonlite::toJSON(json, pretty = T)

# [
#   {
#     "patient" : "A",
#     "samples" : [
#       {
#         "sample" : "P",
#         "files" : [
#           {
#             "file" : "ZZEVYQ"
#           },
#           {
#             "file" : "XRYBUE"
#           }
#         ] 
#       },
#       {
#         "sample" : "R",
#         "files" : [
#           {
#             "file" : "KIUXRU"
#           },
#           {
#             "file" : "ZCHBKN"
#           }
#         ]
#       }
#     ]
#   },
#   {
#     "patient" : "B",
#     "samples" : [
#       {
#         "sample" : "P",
#         "files" : [
#           {
#             "file" : "WZYAPM"
#           },
#           {
#             "file" : "CYEJCK"
#           }
#           ] 
#       },
#       {
#         "sample" : "R",
#         "files" : [
#           {
#             "file" : "EKDFYT"
#           },
#           {
#             "file" : "XFAYXX"
#           }
#         ]
#       }
#     ]
#   }
# ]