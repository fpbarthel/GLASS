library(VariantAnnotation)
library(DBI)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg19)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")
time <- system.time(qres <- dbGetQuery(con, read_file(('sql/variant_status_leeds.sql'))))

ref_genome <- BSgenome.Hsapiens.UCSC.hg19
ref_organism <- GenomeInfoDb::organism(ref_genome)
ref_style <- seqlevelsStyle(ref_genome)

genome_name <- genome(ref_genome)[[1]]
seqlevelsStyle(ref_genome) = "NCBI"

vcf <- VCF(rowRanges = GRanges(seqnames = trimws(qres$chrom),
                        ranges = IRanges(start = as.integer(qres$start_pos),
                                         end = as.integer(qres$end_pos)),
                        seqinfo = seqinfo(ref_genome),
                        paramRangeID = rep(factor(NA),nrow(qres))),
    fixed = DataFrame(REF = DNAStringSet(qres$ref),
                      ALT = unname(split(DNAStringSet(qres$alt),1:length(qres$alt))),
                      QUAL = as.numeric(NA_integer_),
                      FILTER = 'PASS'),
    geno = SimpleList(GT = matrix(rep(".",nrow(qres), dim.names = c(1:nrow(qres), "TEST")))),
    colData = DataFrame(Samples = 1, row.names = c("TEST")))

vcfout <- split(vcf, sprintf("%s-%s", qres$tumor_pair_barcode, qres$variant_status))
#rm(vcf)

for (vcf_name in names(vcfout)) {
  vcf = vcfout[[vcf_name]]
  message("Writing ", vcf_name, " with ", nrow(vcf), " rows.")
  writeVcf(vcf, file = sprintf("results/mutect2/fractionated-vcf/%s.vcf", vcf_name), index = FALSE)
}