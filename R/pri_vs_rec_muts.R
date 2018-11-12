library(tidyverse)
library(odbc)
library(DBI)
library(ggthemes)
library(ggplot2)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

# q = "SELECT gt.aliquot_barcode, ts.case_barcode, idh_status, codel_status, gene_symbol, ts.sample_type, sample_type_count, COUNT(*)
# FROM analysis.called_genotypes gt
# INNER JOIN biospecimen.aliquots ta ON ta.aliquot_barcode = gt.aliquot_barcode
# INNER JOIN biospecimen.samples ts ON ts.sample_barcode = ta.sample_barcode
# INNER JOIN (SELECT sample_type,COUNT(*) AS sample_type_count FROM biospecimen.samples GROUP BY sample_type) stc ON stc.sample_type = ts.sample_type
# LEFT JOIN clinical.surgeries s ON s.sample_barcode = ts.sample_barcode
# WHERE variant_classification IN ('In_Frame_Ins','In_Frame_Del','Frame_Shift_Del','Frame_Shift_Ins','Missense_Mutation','Nonstop_Mutation','Nonsense_Mutation')
# AND gene_symbol IN ('TP53','ATRX','RB1','EGFR','PTEN','TERT','NF1','CIC','FUBP1','PDGFRA','PIK3CA','PIK3R1')
# GROUP BY gt.aliquot_barcode, ts.case_barcode, idh_status, codel_status, gene_symbol, ts.sample_type, sample_type_count
# ORDER BY count desc"

q = "SELECT ts.case_barcode, s.idh_status, s.codel_status, gene_symbol, ts.sample_type, sample_type_count, COUNT(*)
FROM analysis.called_genotypes gt
INNER JOIN biospecimen.aliquots ta ON ta.aliquot_barcode = gt.aliquot_barcode
INNER JOIN biospecimen.samples ts ON ts.sample_barcode = ta.sample_barcode
LEFT JOIN clinical.surgeries s ON s.sample_barcode = ts.sample_barcode
INNER JOIN (SELECT sample_type, idh_status, codel_status, COUNT(*) AS sample_type_count
FROM biospecimen.samples sa
LEFT JOIN clinical.surgeries su ON sa.sample_barcode = su.sample_barcode
GROUP BY sample_type, idh_status, codel_status) stc ON stc.sample_type = ts.sample_type AND stc.idh_status = s.idh_status AND stc.codel_status = s.codel_status
WHERE variant_classification IN ('In_Frame_Ins','In_Frame_Del','Frame_Shift_Del','Frame_Shift_Ins','Missense_Mutation','Nonstop_Mutation','Nonsense_Mutation')
AND gene_symbol IN ('TP53','ATRX','RB1','EGFR','PTEN','NF1','CIC','FUBP1','PDGFRA','PIK3CA','PIK3R1')
GROUP BY ts.case_barcode, s.idh_status, s.codel_status, gene_symbol, ts.sample_type, sample_type_count
ORDER BY count desc"

qres <- dbGetQuery(con, q)

df <- qres %>%
  group_by(case_barcode) %>%
  mutate(idh_codel_grp = ifelse(any(idh_status == "IDH.mt") && any(codel_status == "codel"), "IDHmut-codel", 
                                ifelse(any(idh_status == "IDH.mt"), "IDHmut",
                                       ifelse(any(idh_status == "IDH.wt"), "IDHwt", NA))),
         n = n()) %>%
  ungroup() %>%
  filter(complete.cases(idh_codel_grp)) %>%
  mutate(sample_type = factor(sample_type, c("TP",sprintf("R%s",1:4),"M1"))) %>%
  group_by(gene_symbol, idh_codel_grp, sample_type, sample_type_count) %>%
  summarize(freq = n()/sample_type_count[1]) %>%
  ungroup()# %>% 
  #group_by(gene_symbol) %>%
  #mutate(any_high_freq = any(freq[sample_type=="TP"] > 0.1)) %>%
  #ungroup() %>%
  #filter(any_high_freq)

ggplot(df, aes(x=gene_symbol,y=freq)) + 
  geom_bar(stat="identity", position = "dodge", aes(fill=sample_type)) + 
  facet_wrap(~idh_codel_grp, ncol=1)
