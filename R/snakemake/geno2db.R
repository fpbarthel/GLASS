
library(VariantAnnotation)
library(tidyverse)

## Parse snakemake
if(exists("snakemake")) {
  vcff = snakemake@input[["vcf"]]
  ssvcff = snakemake@input[["ssvcf"]]
  genof = snakemake@output[["geno"]]
  infof = snakemake@output[["info"]]
  case_barcode = snakemake@wildcards[["case_barcode"]]
} else {
  vcff = "results/mutect2/dropgt/GLSS-SF-0001.filtered.normalized.sorted.vcf.gz"
  ssvcff = list("results/mutect2/ssdropgt/GLSS-SF-0001-TP-01-NM-01D-WXS.filtered.normalized.sorted.vcf.gz","results/mutect2/ssdropgt/GLSS-SF-0001-R1-01-NM-01D-WXS.filtered.normalized.sorted.vcf.gz")
  genof = "GLSS-SF-0001_geno.tsv"
  infof = "GLSS-SF-0001_info.tsv"
  case_barcode = "GLSS-SF-0001"
}

## Read freebayes and consensus input as VRanges
vcf = readVcf(vcff, "hg19")
message("Loaded ", basename(vcff))

ssvcf = lapply(ssvcff, function(fn) {
  out <- readVcf(fn, "hg19")
  message("Loaded ", fn)
  return(out[which(rowRanges(out)$FILTER == "PASS")])
  #return(out)
})
names(ssvcf) = sapply(ssvcf, function(v) samples(header(v))[2])
message("Loaded ", paste(ssvcff, collapse = ", "))

## Create output dataframe

vcf_info <- data.frame(case_barcode,
                       chrom = as.character(seqnames(vcf)),
                       pos = sprintf("[%s,%s]", start(vcf), end(vcf)),
                       alt = unstrsplit(CharacterList(alt(vcf)), sep=","),
                       filter = rowRanges(vcf)$FILTER,
                       tlod = unlist(info(vcf)$TLOD),
                       nlod = unlist(info(vcf)$NLOD)
                       stringsAsFactors = F)

# vcf_ssinfo <- lapply(ssvcf, function(v) {
#   data.frame(aliquot_barcode = samples(header(v))[2],
#              chrom = as.character(seqnames(v)),
#              pos = sprintf("[%s,%s]", start(v), end(v)),
#              alt = unstrsplit(CharacterList(alt(v)), sep=","),
#              filter = rowRanges(v)$FILTER,
#              tlod = unlist(info(v)$TLOD),
#              nlod = unlist(info(v)$NLOD),
#              saaf_none = as.numeric(map(as.list(info(v)$SAAF), 3)),
#              stringsAsFactors = F)
# }) %>% data.table::rbindlist() %>% as.data.frame()

extractSample <- function(aliquot_barcode) {
  ssm2_call = NA
  if(aliquot_barcode %in% names(ssvcf))
    ssm2_call = overlapsAny(vcf, ssvcf[[aliquot_barcode]])
  df = data.frame(aliquot_barcode,
                  case_barcode = substr(aliquot_barcode,1,12),
                  chrom = as.character(seqnames(vcf)),
                  pos = sprintf("[%s,%s]", start(vcf), end(vcf)),
                  alt = unstrsplit(CharacterList(alt(vcf)), sep=","),
                  ref_count = as.integer(map(geno(vcf)$AD[,aliquot_barcode],1)),
                  alt_count = as.integer(map(geno(vcf)$AD[,aliquot_barcode],2)),
                  af = as.numeric(map(geno(vcf)$AF[,aliquot_barcode],1)),
                  ssm2_call,
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