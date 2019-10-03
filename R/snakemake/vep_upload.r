#This script takes the output of the annotate_vep rule in the Snakemake mutect2-post.smk module reformats it for uploading to the db (variants.vep table)
#Additionally generates a .tsv file for backup
#-----------------------------------------------------

library(VariantAnnotation)

library(ensemblVEP)
library(tidyverse)
library(DBI)


setwd('/projects/varnf/GLASS/GLASS/')

## Parse snakemake
maff = "results/mutect2/annoconsensusvcf/consensus.normalized.sorted.vep.maf"
vcff = "results/mutect2/consensusvcf/consensus.normalized.sorted.vcf.gz"
tsvf = "results/mutect2/maf2db/consensus.normalized.sorted.vep.tsv"

vcf = readVcf(vcff, "hg19")
maf = read.delim(maff, as.is = T, comment.char = '#')

message("Read file ", basename(vcff))
message("Read file ", basename(maff))

df = data.frame(chrom = as.character(seqnames(vcf)),
				pos = sprintf("[%s,%s]", start(vcf), end(vcf)),
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

#Change chromosome X to chromosome 23
df[which(df[,"chrom"]=='X'),"chrom"] <- 23
df[,"chrom"] = as.numeric(df[,"chrom"])

#Manual edit to match GLASS variant_classifications table; this is now done in SQL
#df[which(df[,"variant_classification"]=="Splice_Region"),"variant_classification"]  <- "Splice_Site"

write.table(df, file = tsvf, quote = F, sep = "\t", row.names = F, col.names = T)

message("Wrote output ", basename(tsvf))

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")
dbWriteTable(con, Id(schema="variants",table="vep"), df, overwrite=TRUE)

