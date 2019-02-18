#Floris Table 1
#--------------------------------------------------

library(BSgenome.Hsapiens.UCSC.hg19)
library(MutationalPatterns)
library(ComplexHeatmap)
library(RColorBrewer)
library(DBI)
library(reshape)

rm(list=ls())
myinf1 <- "/projects/varnf/GLASS/analysis/signatures/GLASS_primary_recur_shared_mut_matrix.rda"

load(myinf1)		#Loads as mut_mat
mut_count <- apply(mut_mat,2,sum)

#Optimal contribution of known signatures (COSMIC)
#---------
#Download mutational signatures from the COSMIC website
sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/", 
		 "signatures_probabilities.txt", sep="")
cancer_signatures <- read.table(sp_url, sep ="\t", header=TRUE)
#Match the order of the mutational types to MutationalPatterns standard
new_order <- match(row.names(mut_mat), cancer_signatures[,"Somatic.Mutation.Type"])
#Reorder cancer signatures data.frame
cancer_signatures <- cancer_signatures[as.vector(new_order),]
#Add trinucleotide changes as row.names
row.names(cancer_signatures) <- cancer_signatures[,"Somatic.Mutation.Type"]
#Keep only 96 contributions of the signatures in matrix
cancer_signatures <- as.matrix(cancer_signatures[,4:33])

#Hierarchically cluster signatures
hclust_cosmic <- cluster_signatures(cancer_signatures, method = "average")
#Store signatures in new order
cosmic_order = colnames(cancer_signatures)[hclust_cosmic[["order"]]]

#Optimal contribution of COSMIC signatures to 96 mutational profiles
individual_fit <- fit_to_signatures(mut_mat, cancer_signatures)
individual_contribution <- individual_fit[["contribution"]]
write.table(individual_contribution, "/projects/varnf/GLASS/analysis/signatures/absolute_contributions_by_private_shared.txt",sep="\t",quote=F,row.names=T)
contribution_sums <- apply(individual_contribution,2,sum)
individual_contribution <- apply(individual_contribution,1,function(x)x/contribution_sums)

#Build Floris's first table
floris_table1 <- melt(individual_contribution)

mutation_status <- sapply(strsplit(as.character(floris_table1[,1]),"-"),function(x)x[9])
floris_table1 <- cbind(floris_table1,mutation_status)

floris_table1[,"mut_count"] <- mut_count[floris_table1[,1]]

floris_table1[,1] <- sapply(strsplit(as.character(floris_table1[,1]),"-"),function(x)paste(x[1:8],collapse="-"))
colnames(floris_table1) <- c("tumor_pair_barcode","signature","relative_contribution","mutation_status","mut_count")

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")
dbWriteTable(con, Id(schema="analysis",table="mutsig_private_vs_shared"), floris_table1, overwrite=TRUE)

write.table(floris_table1, "/projects/varnf/GLASS/analysis/signatures/floris_table1.txt",sep="\t",quote=F,row.names=F)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#Floris Table 2
#--------------------------------------------------

allMutationDbToGRanges <- function(con,
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
            tumor_barcode_a,
            tumor_barcode_b,
            row_number() OVER (PARTITION BY case_barcode ORDER BY surgical_interval_mo DESC, portion_a ASC, portion_b ASC, substring(tumor_pair_barcode from 27 for 3) ASC) AS priority
          FROM analysis.tumor_pairs ps
          LEFT JOIN analysis.blocklist b1 ON b1.aliquot_barcode = ps.tumor_barcode_a
          LEFT JOIN analysis.blocklist b2 ON b2.aliquot_barcode = ps.tumor_barcode_b
          WHERE
            comparison_type = 'longitudinal' AND
            sample_type_b <> 'M1' AND
            b1.coverage_exclusion = 'allow' AND b2.coverage_exclusion = 'allow'
            --b1.cnv_exclusion = 'allow' AND b2.cnv_exclusion = 'allow' AND
            --b1.clinical_exclusion = 'allow' AND b2.clinical_exclusion = 'allow'
        )
        /*
        Aggregate counts over tumor pairs.
        Performs an INNER JOIN with selected_tumor_pairs from above so that only those pairs are included in the analysis.
        Restrict to events with coverage >= 15 in both A and B
        */
        SELECT DISTINCT -- remove duplicate entries
          gt.aliquot_barcode,
          gt.chrom::character varying(2),
          lower(gt.pos) AS start_pos,
          upper(gt.pos) AS end_pos,
          snvs.ref,
          gt.alt
        FROM analysis.genotypes gt
        INNER JOIN selected_tumor_pairs stp ON stp.priority = 1 AND (stp.tumor_barcode_a = gt.aliquot_barcode OR stp.tumor_barcode_b = gt.aliquot_barcode)
        LEFT JOIN analysis.snvs ON snvs.chrom = gt.chrom AND snvs.pos = gt.pos AND snvs.alt = gt.alt
        WHERE mutect2_call AND snvs.variant_type = 'SNP'"
  
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
  vcfout <- split(vcf, qres$aliquot_barcode)
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

vcfs <- allMutationDbToGRanges(con, genome = "BSgenome.Hsapiens.UCSC.hg19", verbose = TRUE)
mut_mat <- mut_matrix(vcfs, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
mut_count <- apply(mut_mat,2,sum)

individual_fit <- fit_to_signatures(mut_mat, cancer_signatures)
individual_contribution <- individual_fit[["contribution"]]
write.table(individual_contribution, "/projects/varnf/GLASS/analysis/signatures/absolute_contributions_by_sample.txt",sep="\t",quote=F,row.names=T)
contribution_sums <- apply(individual_contribution,2,sum)
individual_contribution <- apply(individual_contribution,1,function(x)x/contribution_sums)

floris_table2 <- melt(individual_contribution)
colnames(floris_table2) <- c("tumor_pair_barcode","signature","relative_contribution")
floris_table2[,"mut_count"] <- mut_count[floris_table2[,"tumor_pair_barcode"]]

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")
dbWriteTable(con, Id(schema="analysis",table="mutsig"), floris_table2, overwrite=TRUE)

write.table(floris_table2, "/projects/varnf/GLASS/analysis/signatures/floris_table2.txt",sep="\t",quote=F,row.names=F)
