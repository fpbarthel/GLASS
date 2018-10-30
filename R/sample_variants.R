
library(VariantAnnotation)

sample = "GLSS-SF-0081-TP-01D-WXS-RFB24P"
pair = "GLSS-SF-0081-TP-01-NB-01D-WXS"

fbf = sprintf("results/mutect2/freebayes/%s.normalized.sorted.vcf.gz", sample)
mtf = sprintf("results/mutect2/final/%s.final.vcf", pair)
csf = "results/mutect2/consensusvcf/consensus.normalized.sorted.vcf.gz"

mt = readVcf(mtf, "hg19")
fb = readVcf(fbf, "hg19")
cs = readVcf(csf, "hg19")

hits = findOverlaps(rowRanges(fb), rowRanges(cs), type = "equal", select = "first")
idx = which(!is.na(map))

rowRanges(fb)$hits = hits
