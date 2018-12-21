privateVsSharedMutationDbToGRanges <- function(con,
                                genome = "BSgenome.Hsapiens.UCSC.hg19",
                                verbose = FALSE)
{
  ref_genome <- base::get(genome)
  ref_organism <- GenomeInfoDb::organism(ref_genome)
  ref_style <- seqlevelsStyle(ref_genome)
  
  # Name the VCF's genome as the name of the genome build instead of
  # the BSgenome package name.
  genome_name <- genome(ref_genome)[[1]]
  seqlevelsStyle(ref_genome) = "NCBI"
  
  # Check the class of the reference genome
  if (!(class(ref_genome) == "BSgenome"))
    stop("Please provide the name of a BSgenome object.")
  
  q <- "WITH
        /*
        Define a list of 'selected tumor pairs': select a single primary/recurrent pair per subject/case after removing any aliquots from the blocklist and sorting samples by surgical interval,
        so that pairs with the highest surgical interval are prioritized, and pairs with WGS are prioritized over WXS when both are available
        */
        selected_tumor_pairs AS
        (
          SELECT
            tumor_pair_barcode,
            row_number() OVER (PARTITION BY case_barcode ORDER BY surgical_interval_mo DESC, portion_a ASC, portion_b ASC, substring(tumor_pair_barcode from 27 for 3) ASC) AS priority
          FROM analysis.tumor_pairs ps
          LEFT JOIN analysis.blocklist b1 ON b1.aliquot_barcode = ps.tumor_barcode_a
          LEFT JOIN analysis.blocklist b2 ON b2.aliquot_barcode = ps.tumor_barcode_b
          WHERE
            comparison_type = 'longitudinal' AND
            sample_type_b <> 'M1' AND
            --b1.cnv_exclusion = 'allow' AND b2.cnv_exclusion = 'allow' AND
            b1.coverage_exclusion = 'allow' AND b2.coverage_exclusion = 'allow' --AND
            --b1.clinical_exclusion = 'allow' AND b2.clinical_exclusion = 'allow'
        )
        /*
        Aggregate counts over tumor pairs.
        Performs an INNER JOIN with selected_tumor_pairs from above so that only those pairs are included in the analysis.
        Restrict to events with coverage >= 15 in both A and B
        */
        SELECT
          gtc.tumor_pair_barcode,
          gtc.chrom,
          lower(gtc.pos) AS start_pos,
          upper(gtc.pos) AS end_pos,
          snvs.ref,
          gtc.alt,
          (CASE WHEN mutect2_call_a AND mutect2_call_b THEN 'shared' WHEN mutect2_call_a AND NOT mutect2_call_b THEN 'primary' WHEN mutect2_call_b AND NOT mutect2_call_a THEN 'recurrent' END) AS status
        FROM analysis.master_genotype_comparison gtc
        INNER JOIN selected_tumor_pairs stp ON stp.tumor_pair_barcode = gtc.tumor_pair_barcode AND stp.priority = 1
        LEFT JOIN analysis.snvs snvs ON snvs.chrom = gtc.chrom AND snvs.pos = gtc.pos AND snvs.alt = gtc.alt
        WHERE (mutect2_call_a OR mutect2_call_b) AND read_depth_a >= 15 AND read_depth_b >= 15 AND snvs.variant_type = 'SNP'
        ORDER BY gtc.chrom,gtc.pos"
  
  if(verbose) message("Fetching VCF data from database...")
  time <- system.time(qres <- dbGetQuery(con, q))
  if(verbose) message("Fetch completed in ", round(time[3]), " seconds.")
  
  if(verbose) message("Building VCF object and extracting ranges.")
  vcf <- rowRanges(VCF(rowRanges = GRanges(seqnames = trimws(qres$chrom),
                                           ranges = IRanges(start = as.integer(qres$start_pos),
                                                            end = as.integer(qres$end_pos)-1),
                                           seqinfo = seqinfo(ref_genome),
                                           paramRangeID = rep(factor(NA),nrow(qres))),
                       fixed = DataFrame(REF = DNAStringSet(qres$ref),
                                         ALT = unname(split(DNAStringSet(qres$alt),1:length(qres$alt))),
                                         QUAL = as.numeric(NA_integer_),
                                         FILTER = 'PASS')))
  
  seqlevelsStyle(vcf) = "UCSC"
  
  groups <- c(extractSeqlevelsByGroup(species = ref_organism,
                                      style = ref_style,
                                      group = "auto"),
              extractSeqlevelsByGroup(species = ref_organism,
                                      style = ref_style,
                                      group = "sex"))
  
  groups <- intersect(groups, seqlevels(vcf))
  vcf <- keepSeqlevels(vcf, groups, pruning.mode = "tidy")
  
  if(verbose) message("Splitting VCF object by sample.")
  vcfout <- split(vcf, sprintf("%s-%s", qres$tumor_pair_barcode, qres$status))
  vcfout <- lapply(vcfout, function(gr) { names(gr) <- sprintf("%s:%s_%s/%s", gsub("chr","",as.character(seqnames(gr))), start(gr), gr$REF, unlist(gr$ALT)); return(gr) })
  vcfout <- GRangesList(vcfout)
  if(verbose) message("Done.")
  return(vcfout)
}

library(VariantAnnotation)
library(DBI)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg19)
library(MutationalPatterns)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

test1 <- read_vcfs_as_granges(c("results/mutect2/final/GLSS-SU-0005-TP-01-NB-01D-WXS.final.vcf",
                                "results/mutect2/final/GLSS-SU-0005-R2-01-NB-01D-WXS.final.vcf"),
                              c("GLSS-SU-0005-R1-01-NB-01D-WXS",
                                "GLSS-SU-0005-R2-01-NB-01D-WXS"),
                              genome = "BSgenome.Hsapiens.UCSC.hg19")

test2 <- privateVsSharedMutationDbToGRanges(con, genome = "BSgenome.Hsapiens.UCSC.hg19")

mat1 <- mut_matrix(test1, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")

tmp = test2[sapply(test2,length)>5]
mat2 <- mut_matrix(test2, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
