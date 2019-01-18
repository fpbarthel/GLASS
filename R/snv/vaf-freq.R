library(tidyverse)
library(DBI)
library(ggthemes)
library(ggplot2)
library(RColorBrewer)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")
vaf_res <- dbGetQuery(con, read_file('sql/vaf_compare.sql'))

ggplot(vaf_res, aes(vaf_a, vaf_b)) +
  geom_point(aes(color=variant_classification)) + 
  geom_abline(slope=1, alpha=0.2, linetype=2) +
  labs(x="Primary VAF", y="Optimal recurrence VAF", color = "Variant Classification") +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  theme_bw(base_size = 18) +
  theme(axis.text=element_text(size=10)) +
  facet_wrap(~gene_symbol) + 
  scale_color_manual(values=c("5'Flank" = brewer.pal(9, "Paired")[9],
                              "Frame_Shift_Del" = brewer.pal(7, "Paired")[1],
                              "Frame_Shift_Ins" = brewer.pal(7, "Paired")[2],
                              "In_Frame_Del" = brewer.pal(7, "Paired")[3],
                              "In_Frame_Ins" = brewer.pal(7, "Paired")[4],
                              "Missense_Mutation" = brewer.pal(7, "Paired")[5],
                              "Nonsense_Mutation" = brewer.pal(7, "Paired")[6],
                              "Splice_Site" = brewer.pal(9, "Paired")[7],
                              "Translation_Start_Site" = brewer.pal(9, "Paired")[8]))
