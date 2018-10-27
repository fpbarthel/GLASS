
sample = "GLSS-MD-0138-R1-01D-WXS-JT1ZPW"

fbf = sprintf("results/mutect2/freebayes/%s.vcf.gz", sample)
mtf = sprintf("results/mutect2/final/%s.final.vcf", sample)

mt = readVcf(mtf, "hg19")