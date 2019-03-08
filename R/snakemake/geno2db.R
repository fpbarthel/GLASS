
library(VariantAnnotation)
library(tidyverse)

## Parse snakemake
if(exists("snakemake")) {
  vcff = snakemake@input[["vcf"]]
  ssvcff = snakemake@input[["ssvcf"]]
  mfcvf = snakemake@output[["mfcov"]]
  genof = snakemake@output[["geno"]]
  infof = snakemake@output[["info"]]
  case_barcode = snakemake@wildcards[["case_barcode"]]
} else {
  vcff = "results/mutect2/dropgt/GLSS-DF-0014.filtered.normalized.sorted.vcf.gz"
  ssvcff = list("results/mutect2/ssdropgt/GLSS-DF-0014-R1-01-NB-01D-WXS.filtered.normalized.sorted.vcf.gz","results/mutect2/ssdropgt/GLSS-DF-0014-TP-01-NB-01D-WXS.filtered.normalized.sorted.vcf.gz")
  mfcvf = "GLSS-DF-0014_mfcvf.tsv"
  genof = "GLSS-DF-0014_geno.tsv"
  infof = "GLSS-DF-0014_info.tsv"
  case_barcode = "GLSS-DF-0014"
  
  vcff = "results/mutect2/dropgt/GLSS-MD-0135.filtered.normalized.sorted.vcf.gz"
  ssvcff = list("results/mutect2/ssdropgt/GLSS-MD-0135-R1-01-NB-01D-WXS.filtered.normalized.sorted.vcf.gz")
  mfcvf = "GLSS-MD-0135_mfcvf.tsv"
  genof = "GLSS-MD-0135_geno.tsv"
  infof = "GLSS-MD-0135_info.tsv"
  case_barcode = "GLSS-MD-0135"
}

## Read multi-sample M2 into CollapsedVCF object
vcf = readVcf(vcff, "hg19")
vr = as(vcf, "VRanges")
message("Loaded ", basename(vcff))

## Read single-sample M2 into list of CollapsedVCF objects
ssvcf = lapply(ssvcff, function(fn) {
  out <- readVcf(fn, "hg19")
  message("Loaded ", fn)
  return(out[which(rowRanges(out)$FILTER == "PASS")]) # Single-sample calls are for validation only
})
names(ssvcf) = sapply(ssvcf, function(v) samples(header(v))[2])

## remove empty
idx = which(sapply(ssvcf,length) == 0)
ssvcf = ssvcf[-idx]

ssvr = lapply(ssvcf, as, "VRanges")
message("Loaded ", paste(ssvcff, collapse = ", "))

## Perform some sanity checks
if(!exists("snakemake")) {
  names(ssvcf)
  length(vcf)
  
  length(ssvcf[[1]])
  table(overlapsAny(ssvcf[[1]], vcf)) ## R1
  table(overlapsAny(subset(ssvr[[1]],sampleNames=="GLSS-DF-0014-R1-01D-WXS-ZZNMIY"), subset(vr, sampleNames=="GLSS-DF-0014-R1-01D-WXS-ZZNMIY")))
  table(subset(ssvr[[1]],sampleNames=="GLSS-DF-0014-R1-01D-WXS-ZZNMIY") %in% subset(vr, sampleNames=="GLSS-DF-0014-R1-01D-WXS-ZZNMIY"))

  length(ssvcf[[2]])
  table(overlapsAny(ssvcf[[2]], vcf)) ## TP
  table(overlapsAny(subset(ssvr[[2]],sampleNames=="GLSS-DF-0014-TP-01D-WXS-H9KDMM"), subset(vr, sampleNames=="GLSS-DF-0014-TP-01D-WXS-H9KDMM")))
  table(subset(ssvr[[2]],sampleNames=="GLSS-DF-0014-TP-01D-WXS-H9KDMM") %in% subset(vr, sampleNames=="GLSS-DF-0014-TP-01D-WXS-H9KDMM"))
  
  subset(mfcov, dp == 14)
  #length(geno())
  #table(sapply(geno(ssvcf[[1]])$AD[,2],sum) > 14)
  table(sapply(geno(ssvcf[[2]])$AD[,2],sum) > 14)
  table(overlapsAny(ssvcf[[2]], vcf))
  
  table(overlapsAny(vcf[sapply(geno(vcf)$AD[,"GLSS-DF-0014-TP-01D-WXS-H9KDMM"],sum) > 14], ssvcf[["GLSS-DF-0014-TP-01D-WXS-H9KDMM"]]))
  table(subset(ssvr[[2]],sampleNames="GLSS-DF-0014-TP-01D-WXS-H9KDMM" && refDepth(ssvr[[2]]) + altDepth(ssvr[[2]]) > 14) %in% subset(vr, sampleNames=="GLSS-DF-0014-TP-01D-WXS-H9KDMM"))
  table(subset(vr, sampleNames(vr) == "GLSS-DF-0014-TP-01D-WXS-H9KDMM" && refDepth(vr) + altDepth(vr) > 14) %in% subset(ssvr[[2]],sampleNames="GLSS-DF-0014-TP-01D-WXS-H9KDMM" && refDepth(ssvr[[2]]) + altDepth(ssvr[[2]]) > 14))
}

## Create output dataframe
vcf_info <- data.frame(case_barcode,
                       chrom = as.character(seqnames(vcf)),
                       pos = sprintf("[%s,%s]", start(vcf), end(vcf)),
                       alt = unstrsplit(CharacterList(alt(vcf)), sep=","),
                       filter = rowRanges(vcf)$FILTER,
                       tlod = unlist(info(vcf)$TLOD),
                       nlod = unlist(info(vcf)$NLOD),
                       stringsAsFactors = F)

extractSample <- function(aliquot_barcode) {
  
  if(aliquot_barcode %in% names(ssvcf)) {
    
    ## Define single same variant object
    ssv = subset(ssvr[[aliquot_barcode]], sampleNames(ssvr[[aliquot_barcode]]) == aliquot_barcode)
    
    ## Define multi-sample variant object
    msv = subset(vr, sampleNames(vr) == aliquot_barcode)
    
    idx = match(msv,ssv)
    
    df = data.frame(aliquot_barcode,
                    case_barcode = substr(aliquot_barcode,1,12),
                    chrom = as.character(seqnames(msv)),
                    pos = sprintf("[%s,%s]", start(msv), end(msv)),
                    alt = alt(msv),
                    ref_count = refDepth(msv),
                    alt_count = altDepth(msv),
                    af = as.numeric(msv$AF),
                    ssm2_pass_call = !is.na(idx),
                    ssm2_saaf_none = as.numeric(map(as.list(ssv$SAAF),3))[idx],
                    stringsAsFactors = F)
  
  } else {
    
    df = data.frame(aliquot_barcode,
                    case_barcode = substr(aliquot_barcode,1,12),
                    chrom = as.character(seqnames(vcf)),
                    pos = sprintf("[%s,%s]", start(vcf), end(vcf)),
                    alt = unstrsplit(CharacterList(alt(vcf)), sep=","),
                    ref_count = as.integer(map(geno(vcf)$AD[,aliquot_barcode],1)),
                    alt_count = as.integer(map(geno(vcf)$AD[,aliquot_barcode],2)),
                    af = as.numeric(map(geno(vcf)$AF[,aliquot_barcode],1)),
                    ssm2_pass_call = NA,
                    ssm2_saaf_none = NA,
                    stringsAsFactors = F)
  }

  return(df)
}
vcf_geno <- data.table::rbindlist(lapply(rownames(colData(vcf)), extractSample))

mfcov <- lapply(ssvcf, function(v) {
  t_adsum = sapply(geno(v)$AD[,2],sum)
  data.frame(aliquot_barcode = samples(header(v))[2],
             dp = seq(0:250),
             nmut = sapply(seq(0:250), function(i) sum(t_adsum > i)))
}) %>% data.table::rbindlist() %>% as.data.frame()

write.table(vcf_geno, file = genof, quote = F, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(mfcov,    file = mfcvf, quote = F, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(vcf_info, file = infof, quote = F, sep = "\t", row.names = FALSE, col.names = FALSE)

## END ##