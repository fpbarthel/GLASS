##################################################
# Perform a correlation test between the PyClone mutational cluster rankings
# Updated: 2019.07.14
# Author: Kevin J.
##################################################

# Working directory for this analysis in the GLASS-analysis project. 
mybasedir = "/Volumes/verhaak-lab/GLASS-analysis/"
setwd(mybasedir)

#######################################################

# Necessary packages:
library(tidyverse)
library(DBI)

#######################################################
# Establish connection with the GLASS database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB2")


# Query to generate PyClone ranks based on mean CCF of the cluster.
pyclone_rank = dbGetQuery(con, "WITH
selected_tumor_pairs AS
                                       (
                                       SELECT tumor_pair_barcode, tumor_barcode_a, tumor_barcode_b, ss.case_barcode, idh_codel_subtype
                                       FROM analysis.gold_set ss
                                       INNER JOIN clinical.subtypes st ON st.case_barcode = ss.case_barcode
                                       ),
                                       pyclone_clusters AS
                                       (
                                       SELECT stp.case_barcode, stp.idh_codel_subtype, pc1.cluster_id, pc1.size AS size, pc1.mean AS ccf_a, pc2.mean AS ccf_b,
                                       (RANK() OVER (PARTITION BY stp.case_barcode ORDER BY pc1.mean DESC))::integer AS rank_a,
                                       (RANK() OVER (PARTITION BY stp.case_barcode ORDER BY pc2.mean DESC))::integer AS rank_b
                                       FROM selected_tumor_pairs stp
                                       INNER JOIN variants.pyclone_cluster pc1 ON pc1.aliquot_barcode = stp.tumor_barcode_a
                                       INNER JOIN variants.pyclone_cluster pc2 ON pc2.aliquot_barcode = stp.tumor_barcode_b AND pc2.cluster_id = pc1.cluster_id
                                       WHERE pc1.size > 1 AND (pc1.mean > 0.1 OR pc2.mean > 0.1)
                                       )
                                       SELECT * FROM pyclone_clusters")

# Test whether the ccf_rankings are strongly associated with one another for all gold_set samples.
res = cor.test(pyclone_rank$rank_a, pyclone_rank$rank_b, method="kendall")

# Separate into subtypes.
pyclone_rank_IDHwt = pyclone_rank %>% 
  filter(idh_codel_subtype=="IDHwt")
pyclone_rank_IDHnoncodel = pyclone_rank %>% 
  filter(idh_codel_subtype=="IDHmut-noncodel")
pyclone_rank_IDHcodel = pyclone_rank %>% 
  filter(idh_codel_subtype=="IDHmut-codel")

# Apply cor.test per subtype.
cor.test(pyclone_rank_IDHwt$rank_a, pyclone_rank_IDHwt$rank_b, method="kendall")
cor.test(pyclone_rank_IDHnoncodel$rank_a, pyclone_rank_IDHnoncodel$rank_b, method="kendall")
cor.test(pyclone_rank_IDHcodel$rank_a, pyclone_rank_IDHcodel$rank_b, method="kendall")

# Remove the first, potentially truncal mutational cluster.
pyclone_rank_non_trunk = pyclone_rank %>% 
    filter(rank_a!="1")
# Is there still a significant, albeit weak relationship between initial and recurrence?
cor.test(pyclone_rank_non_trunk$rank_a, pyclone_rank_non_trunk$rank_b, method="kendall")

# Remove the first, truncal mutational cluster and perform ccf analysis by subtypes.
pyclone_rank_noT_IDHwt = pyclone_rank_non_trunk %>% 
  filter(idh_codel_subtype=="IDHwt")
pyclone_rank_noT_IDHnoncodel = pyclone_rank_non_trunk %>% 
  filter(idh_codel_subtype=="IDHmut-noncodel")
pyclone_rank_noT_IDHcodel = pyclone_rank_non_trunk %>% 
  filter(idh_codel_subtype=="IDHmut-codel")

# Apply cor.test per subtype.
cor.test(pyclone_rank_noT_IDHwt$rank_a, pyclone_rank_noT_IDHwt$rank_b, method="kendall")
cor.test(pyclone_rank_noT_IDHnoncodel$rank_a, pyclone_rank_noT_IDHnoncodel$rank_b, method="kendall")
cor.test(pyclone_rank_noT_IDHcodel$rank_a, pyclone_rank_noT_IDHcodel$rank_b, method="kendall")
