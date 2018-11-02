
library(VariantAnnotation)
#setwd("/fastscratch/verhaak-lab/GLASS-WG")

## Parse snakemake
fbf = snakemake@input[["freebayes"]]
csf = snakemake@input[["consensus"]]
mtf = snakemake@params[["mutect2"]]
tsv = snakemake@output[["tsv"]]
spl = snakemake@wildcards[["aliquot_barcode"]]

## Read freebayes and consensus input as VRanges
fb = readVcfAsVRanges(fbf, "hg19")
cs = readVcfAsVRanges(csf, "hg19", param=ScanVcfParam(fixed = "ALT", info = NA, geno = "AD"))
message("Loaded ", basename(fbf))
message("Loaded ", basename(csf))

## If sample is a tumor sample, read mutect calls
paired = FALSE
if(!is.na(mtf) & file.exists(mtf)) {
  paired = TRUE
  mt = readVcfAsVRanges(mtf, "hg19", param=ScanVcfParam(fixed = "ALT", info = NA, geno = "AD"))
  message("Loaded ", basename(mtf))
}

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
rm(cs)

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