library(tidyverse)
library(stringi)

df = data.frame(patient  = c(rep("A",4), rep("B",4)),
                sample    = rep(c("P","R"),4),
                file      = stri_rand_strings(8, 6, '[A-Z]'))

dfout <- df %>% group_by(patient, sample) %>% 
  summarize(files=list(map(file, ~list(file=.x)))) %>% 
  summarize(samples=list(map2(sample, files, ~list(samples=.x, files=.y))))

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