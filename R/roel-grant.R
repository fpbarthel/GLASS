library(maftools)
library(tidyverse)
setwd('/fastscratch/verhaak-lab/GLASS-WG/')

genedf = openxlsx::read.xlsx("data/ref/glioma_driver_genes.xlsx")
anno = openxlsx::read.xlsx("data/ref/glass_wg_subtypes.xlsx")
codel = read.delim("data/ref/glass_1p19q_codeletion_status.txt", as.is=T) %>% select(aliquot_id=sample_id,codel_status)
anno = anno %>%
  mutate(sample_type = factor(substr(anno$aliquot_id,14,15), labels = c("TP", sprintf("R%s",1:4), "M1"))) %>%
  filter(sample_type != "NB") %>%
  left_join(codel) %>%
  mutate(aliquot_id = gsub("-", ".", aliquot_id),
         IDHCodel = ifelse(codel_status=='codel', 'IDHmut-codel', ifelse(Sequencing.IDH=='Mutant', 'IDHmut', ifelse(Sequencing.IDH=='WT', 'IDHwt', NA)))) %>%
  select(Tumor_Sample_Barcode=aliquot_id, sample_type, IDHCodel)

mafdf = read_tsv('results/mutect2/vcf2maf/GLASS.maf',skip=1)
mafdf = mafdf %>%
  filter(Tumor_Sample_Barcode %in% anno$aliquot_id,
         Hugo_Symbol %in% genedf$gene) %>%
  mutate(Tumor_Sample_Barcode = gsub("-", ".", Tumor_Sample_Barcode))

cndf = read_tsv('data/ref/maf_tools_input_wgs.txt')
cndf = cndf %>%
  rename(Tumor_Sample_Barcode = sample_id) %>%
  mutate(Tumor_Sample_Barcode = gsub("-", ".", Tumor_Sample_Barcode))

maf = read.maf(mafdf, cnTable = , removeSilent = TRUE, isTCGA = FALSE)
ssmaf = subsetMaf(maf, genes = genedf$gene, mafObj = TRUE)

fabcolors = RColorBrewer::brewer.pal(n = 6,name = 'Spectral')
names(fabcolors) = c(levels(maf@data$sample_type))

pdf('test.pdf', width = 12, height = 12)
plotmafSummary(ssmaf)
plot.new()
oncoplot(ssmaf, top=20, annotation = anno, sortByAnnotation = TRUE) #annotationColor = fabcolors, 
dev.off()


