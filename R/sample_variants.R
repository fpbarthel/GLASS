
library(VariantAnnotation)

sample = "GLSS-SF-0081-TP-01D-WXS-RFB24P"
pair = "GLSS-SF-0081-TP-01-NB-01D-WXS"

setwd("/fastscratch/verhaak-lab/GLASS-WG")

fbf = sprintf("results/mutect2/freebayes/%s.normalized.sorted.vcf.gz", sample)
mtf = sprintf("results/mutect2/m2post/%s.normalized.sorted.vcf.gz", pair)
csf = "results/mutect2/consensusvcf/consensus.normalized.sorted.vcf.gz"

fb = readVcfAsVRanges(fbf, "hg19")
cs = readVcfAsVRanges(csf, "hg19", param=ScanVcfParam(fixed = "ALT", info = NA, geno = "AD"))
mt = readVcfAsVRanges(mtf, "hg19", param=ScanVcfParam(fixed = "ALT", info = NA, geno = "AD"))

message("Loaded ", basename(fbf))
message("Loaded ", basename(csf))

## Count overlap
hitscs = fb %in% cs

prop_cs = round(length(fb)/length(cs)*100,1)
prop_fb = round(sum(hitscs)/length(hitscs)*100,1)
message("Found ", length(cs), " variants in consensus callset")
message("Found ", length(fb), " freebayes calls (", prop_cs, "% of callset), amongst which ",
        sum(hitscs), " (", prop_fb, "%) matched calls from the consensus callset.")

hitsmt = fb %in% mt
prop_mt = round(sum(hitsmt)/length(mt)*100,1)
message("Found ", length(mt), " filtered Mutect2 calls, of which ", sum(hitsmt), " (",
        prop_mt, "%) exactly match calls from the consensus callset.")

## Subset freebayes output by only variants present in db
fb = fb[hitscs]

## Annotate M2-called variants
fb$called = hitsmt

## Create output dataframe
df = data.frame(aliquot_barcode = sample, chrom = seqnames(fb), start = start(fb), end = end(fb), alt = alt(fb),
                genotype = fb$GT, read_depth = totalDepth(fb), ref_count = refDepth(fb), alt_count = altDepth(fb),
                called = fb$called,
                stringsAsFactors = F)

## Clear some memory
rm(cs, fb, mt)

## Write to database
