# To get mutations in cell cycle genes from Figure 1.
mutation_genes = read_file("sql/heatmap/build_heatmap_data_mutation.sql")
mutations_selected = dbGetQuery(con, mutation_genes)

# To query genes that are found to be altered in the cell cycle.
cell_cycle_cnv_titan = dbGetQuery(con, "SELECT * FROM analysis.cnv_by_gene WHERE gene_symbol IN ('CDK4','CCND2','CDK6','CDKN2A','RB1')")
cell_cycle_cnv_gatk = dbGetQuery(con, "SELECT * FROM analysis.cnv_by_gene_gatk WHERE gene_symbol IN ('CDK4','CCND2','CDK6','CDKN2A','RB1')")
cell_cycle_cnv_titan = dbGetQuery(con, "SELECT * FROM analysis.cnv_by_gene WHERE gene_symbol = 'CDKN2A'")
cell_cycle_cnv_gatk = dbGetQuery(con, "SELECT * FROM analysis.cnv_by_gene_gatk WHERE gene_symbol = 'CDKN2A'")

cnv_titan = cell_cycle_cnv_titan %>% 
  inner_join(pairs, by="pair_barcode") %>% 
  select(tumor_barcode, gene_symbol, copy_number, corrected_cn, titan_call)

cell_cycle_cnv_merged = cell_cycle_cnv_gatk %>% 
  select(aliquot_barcode, corrected_cn_gatk = corrected_cn, cn_call_gatk = cn_call) %>% 
  inner_join(cnv_titan, by=c("aliquot_barcode"="tumor_barcode"))

table(cell_cycle_cnv_merged$cn_call_gatk, cell_cycle_cnv_merged$titan_call)