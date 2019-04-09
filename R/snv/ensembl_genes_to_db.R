library(biomaRt)
library(DBI)

ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
term = getBM(c('ensembl_gene_id','hgnc_symbol'),mart=ensembl)

con <- DBI::dbConnect(odbc::odbc(), "GLASSv2")
dbWriteTable(con, Id(schema='ref',table='ensembl_genes'),term)
