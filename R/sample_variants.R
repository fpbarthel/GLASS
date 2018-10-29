
sample = "GLSS-MD-0138-R1-01D-WXS-JT1ZPW"
pair = "GLSS-MD-0138-R1-01-NB-01D-WXS"

fbf = sprintf("results/mutect2/freebayes/%s.vcf.gz", sample)
mtf = sprintf("results/mutect2/final/%s.final.vcf", pair)

mt = readVcf(mtf, "hg19")
fb = readVcf(fbf, "hg19")

test = findOverlaps(fb, mt, type = "equal")
