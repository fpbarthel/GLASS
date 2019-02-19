
library(VariantAnnotation)
library(tidyverse)

## Parse snakemake
if(exists("snakemake")) {
  vcff = snakemake@input[["vcf"]]
  genof = snakemake@output[["geno"]]
  infof = snakemake@output[["info"]]
  case_barcode = snakemake@wildcards[["case_barcode"]]
} else {
  vcff = "results/mutect2/dropgt/GLSS-SF-0001.filtered.normalized.sorted.vcf.gz"
  genof = "GLSS-SF-0001_geno.tsv"
  infof = "GLSS-SF-0001_info.tsv"
  case_barcode = "GLSS-SF-0001"
}

## Read freebayes and consensus input as VRanges
vcf = readVcf(vcff, "hg19")
message("Loaded ", basename(vcff))

## Create output dataframe

vcf_info <- data.frame(case_barcode,
                       chrom = as.character(seqnames(vcf)),
                       pos = sprintf("[%s,%s]", start(vcf), end(vcf)),
                       alt = unstrsplit(CharacterList(alt(vcf)), sep=","),
                       filter = rowRanges(vcf)$FILTER,
                       tlod = unlist(info(vcf)$TLOD),
                       nlod = unlist(info(vcf)$NLOD),
                       saaf_none = as.numeric(map(as.list(info(vcf)$SAAF), 3)),
                       stringsAsFactors = F)

extractSample <- function(aliquot_barcode) {
  df = data.frame(aliquot_barcode,
                  chrom = as.character(seqnames(vcf)),
                  pos = sprintf("[%s,%s]", start(vcf), end(vcf)),
                  alt = unstrsplit(CharacterList(alt(vcf)), sep=","),
                  ref_count = as.integer(map(geno(vcf)$AD[,aliquot_barcode],1)),
                  alt_count = as.integer(map(geno(vcf)$AD[,aliquot_barcode],2)),
                  stringsAsFactors = F)
  return(df)
}

vcf_geno <- data.table::rbindlist(lapply(rownames(colData(vcf)), extractSample))

write.table(vcf_geno, file = genof, quote = F, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(vcf_info, file = infof, quote = F, sep = "\t", row.names = FALSE, col.names = FALSE)

## Write to database
# .libPaths('/home/barthf/R/x86_64-pc-linux-gnu-library/3.3')

#con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")
#DBI::dbWriteTable(con, DBI::Id(schema="analysis",table="snv_genotypes"), df, append=T)

## END ##