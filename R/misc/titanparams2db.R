### push titan params to db

library(tidyverse)
library(DBI)
library(odbc)
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

paramfiles <- list.files('results/cnv/titanfinal/params', full.names = TRUE)

## taken from Gavin's script
formatParams <- function(params){
  id <- colnames(params)
  barcode <- strsplit(id, "_cluster")[[1]][1]
  cellPrev <- strsplit(params[grepl("Clonal cluster cellular prevalence", 
                                    rownames(params)), 1], " ")[[1]]
  numClust <- length(cellPrev)
  cellPrev <- paste0(format(cellPrev, digits=4), collapse=",")
  norm <- as.numeric(params[grepl("Normal contamination estimate", rownames(params)), 1])
  purity <- 1 - norm
  ploidy <- as.numeric(params[grepl("Average tumour ploidy estimate", rownames(params)), 1])
  loglik <- as.numeric(params[grepl("likelihood", rownames(params)), 1])
  sdbw <- as.numeric(params[grepl("S_Dbw validity index \\(Both\\)", rownames(params)), 1])
  return(data.frame(id=id, barcode=barcode, numClust=numClust, cellPrev=cellPrev, 
              purity=purity, norm=norm, ploidy=ploidy, loglik=loglik, sdbw=sdbw,
              stringsAsFactors = FALSE))
}

datlist = lapply(paramfiles, function(f) {
  phi <- read.delim(f, header=F, row.names=1, stringsAsFactors=F, sep="\t")
  colnames(phi) <- gsub(".params.txt", "", basename(f))	
  return(formatParams(phi))
})

dat = data.table::rbindlist(datlist) %>% 
  as.data.frame() %>%
  select(pair_barcode = barcode,
         num_clones = numClust,
         cellular_prevalence = cellPrev,
         purity,
         normal_contamination = norm,
         ploidy,
         loglik,
         sdbw)

dbWriteTable(con, Id(schema="analysis",table="titan_params"), dat, append = FALSE)
