
library(VariantAnnotation)
#setwd("/fastscratch/verhaak-lab/GLASS-WG")

## Parse snakemake
f1f = snakemake@input[["freebayes1"]] # "results/mutect2/freebayes/batch1/GLSS-MD-LP10-R1-01D-WGS-CYN7CI.normalized.sorted.vcf.gz" #
f2f = snakemake@input[["freebayes2"]] # "results/mutect2/freebayes/batch2/GLSS-MD-LP10-R1-01D-WGS-CYN7CI.normalized.sorted.vcf.gz" #
f3f = snakemake@input[["freebayes3"]] # "results/mutect2/freebayes/batch3/GLSS-MD-LP10-R1-01D-WGS-CYN7CI.normalized.sorted.vcf.gz" #
f4f = snakemake@input[["freebayes4"]] # "results/mutect2/freebayes/batch4/GLSS-MD-LP10-R1-01D-WGS-CYN7CI.normalized.sorted.vcf.gz" #
f5f = snakemake@input[["freebayes5"]] # "results/mutect2/freebayes/batch5/GLSS-MD-LP10-R1-01D-WGS-CYN7CI.normalized.sorted.vcf.gz" #
csf = snakemake@input[["consensus"]] # "results/mutect2/consensusvcf/consensus.normalized.sorted.vcf.gz" #
mtf = snakemake@params[["mutect2"]] # "results/mutect2/m2post/GLSS-MD-LP10-R1-01-NB-01D-WXS.normalized.sorted.vcf.gz" #
tsv = snakemake@output[["tsv"]]
spl = snakemake@wildcards[["aliquot_barcode"]] # "GLSS-MD-LP10-R1-01D-WGS-CYN7CI" #

## Read freebayes and consensus input as VRanges
f1 = readVcfAsVRanges(f1f, "hg19")
message("Loaded ", basename(f1f))

f2 = readVcfAsVRanges(f2f, "hg19")
message("Loaded ", basename(f2f))

f3 = readVcfAsVRanges(f3f, "hg19")
message("Loaded ", basename(f3f))

f4 = readVcfAsVRanges(f4f, "hg19")
message("Loaded ", basename(f4f))

f5 = readVcfAsVRanges(f5f, "hg19")
message("Loaded ", basename(f5f))

cs = readVcfAsVRanges(csf, "hg19", param=ScanVcfParam(fixed = "ALT", info = NA, geno = "AD"))
message("Loaded ", basename(csf))

## If sample is a tumor sample, read mutect calls
paired = FALSE
if(!is.na(mtf) & file.exists(mtf)) {
  paired = TRUE
  mt = readVcfAsVRanges(mtf, "hg19", param=ScanVcfParam(fixed = "ALT", info = NA, geno = "AD"))
  message("Loaded ", basename(mtf))
}

# table(overlapsAny(f2,f3))
# table(f2 %in% f3)
# table(f3 %in% f2)
# subset(f2, f2 %in% f3)
# subset(f3, f3 %in% f2)
# 
# table(overlapsAny(f2,f4))
# table((f2 %in% f4))
# 
# table(overlapsAny(f2,f5))
# table((f2 %in% f5))
# 
# table(overlapsAny(f3,f4))
# table((f3 %in% f4))
# 
# table(overlapsAny(f3,f5))
# table((f3 %in% f5))
# 
# table(overlapsAny(f4,f5))
# table((f4 %in% f5))

## Merge new freebayes calls
fb = c(f2,f3,f4,f5)
fb = fb[!duplicated(fb)]

## Drop calls already present in batch1
fb = fb[!(fb %in% f1)]

## Count overlap between freebayes calls and consensus callset
hitscs = fb %in% cs

## Print some numbers
prop_cs = round(length(fb)/length(cs)*100,1)
prop_fb = round(sum(hitscs)/length(hitscs)*100,1)
message("Found ", length(cs), " variants in consensus callset")
message("Found ", length(fb), " freebayes calls (", prop_cs, "% of callset), amongst which ",
        sum(hitscs), " (", prop_fb, "%) matched calls from the consensus callset.")

## Subset freebayes output by only variants present in db
fb = fb[hitscs]

## Clear some memory
rm(cs,f1)

## If a mutect callset is available, quantify overlap between mutect and freebayes
if(paired) {
  hitsmt = fb %in% mt
  prop_mt = round(sum(hitsmt)/length(mt)*100,1)
  message("Found ", length(mt), " filtered Mutect2 calls, of which ", sum(hitsmt), " (",
          prop_mt, "%) exactly match calls from the consensus callset.")
  
  ## Annotate M2-called variants
  fb$called = hitsmt
} else {
  fb$called = FALSE
}

## Create output dataframe
df = data.frame(aliquot_barcode = spl, chrom = seqnames(fb), start = start(fb), end = end(fb), alt = alt(fb),
                genotype = fb$GT, read_depth = totalDepth(fb), ref_count = refDepth(fb), alt_count = altDepth(fb),
                called = fb$called,
                stringsAsFactors = F)

## Drop variants without read counts
df = df[which(!is.na(df$read_depth)),]

## Clear more memory
rm(fb)

write.table(df, file = tsv, quote = F, sep = "\t", row.names = F, col.names = T)

## Write to database
# .libPaths('/home/barthf/R/x86_64-pc-linux-gnu-library/3.3')

#con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")
#DBI::dbWriteTable(con, DBI::Id(schema="analysis",table="snv_genotypes"), df, append=T)

## Write a trigger with the number of rows added
#cat(nrow(df), file = trf)
#message("Printed number of rows (", nrow(df), ") to file: ", basename(trf))

## END ##