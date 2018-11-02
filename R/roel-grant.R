library(maftools)
setwd('/fastscratch/verhaak-lab/GLASS-WG/')
maf = read.maf('results/mutect2/vcf2maf/GLASS.maf')

pdf('test.pdf', width = 12, height = 12)
plotmafSummary(maf)
plot.new()
oncoplot(maf, top=20)
geneCloud(maf)
dev.off()

?read.maf
