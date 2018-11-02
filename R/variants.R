#df = read.delim('/fastscratch/verhaak-lab/GLASS-WG/results/mutect2/pon/pon.coding_only.anno.tsv', as.is=T, sep="\t")
#df2 = read.delim('/fastscratch/verhaak-lab/GLASS-C7/results/cnv/igv_convert/GLSS-SM-R095-R1-01-NB-01D-WXS.igv.seg',as.is=T)

#variants = data.table::fread('/fastscratch/verhaak-lab/GLASS-WG/results/mutect2/annoconsensusvcf/consensus.tsv', skip=49)

#newvar = variants %>% separate(Location, c("chrom", "pos"), sep=":", remove = T, convert = T) #mutate(chrom)
#newvar = newvar %>% select(chrom, pos, ref = )

library(VariantAnnotation)
library(ensemblVEP)
library(tidyverse)
library(DBI)

setwd('/fastscratch/verhaak-lab/GLASS-WG/')

vcf = readVcf("results/mutect2/consensusvcf/consensus.normalized.sorted.vcf.gz", "hg19")
maf = read.delim("results/mutect2/annoconsensusvcf/consensus.normalized.sorted.maf", as.is = T, comment.char = '#')


df = data.frame(chrom = seqnames(vcf),
                start = start(vcf),
                end = end(vcf),
                ref = ref(vcf),
                alt = unstrsplit(CharacterList(alt(vcf)), sep=","),
                gene_id = maf$Gene,
                gene_symbol = maf$Hugo_Symbol,
                variant_classification = maf$Variant_Classification,
                variant_type = maf$Variant_Type,
                cdna_position = maf$cDNA_position,
                cds_position = maf$CDS_position,
                protein_position = maf$Protein_position,
                amino_acids = maf$Amino_acids, 
                codons = maf$Codons,
                hgvs_c = maf$HGVSc,
                hgvs_p = maf$HGVSp_Short,
                polyphen = maf$PolyPhen,
                sift = maf$SIFT,
                stringsAsFactors = F)

.libPaths('/home/barthf/R/x86_64-pc-linux-gnu-library/3.3')

idx = duplicated(df)
# > table(idx)
# idx
#   FALSE    TRUE 
# 2687020      27 
df = df[-idx,]

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")
dbWriteTable(con, Id(schema="analysis",table="snvs"), df, append=T)

## END ##