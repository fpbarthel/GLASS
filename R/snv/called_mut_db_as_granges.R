mutationDbToGRanges <- function(con,
                                aliquot_barcodes,
                                genome = "BSgenome.Hsapiens.UCSC.hg19",
                                verbose = FALSE,
                                analysis = c("analysis.called_genotypes", "analysis.called_genotypes_titan_cluster"))
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
  
  
  if (analysis[1] == "analysis.called_genotypes")
  {
    q <- "SELECT chrom, start, \"end\", ref, alt,
          row_number() OVER (PARTITION BY chrom, start, \"end\" ORDER BY alt_count DESC) numv
          FROM analysis.called_genotypes
          WHERE aliquot_barcode = ? AND variant_type = 'SNP';"
  }
  else if (analysis[1] == "analysis.called_genotypes_titan_cluster")
  {
    q <- "SELECT aliquot_barcode_cluster, chrom, start, \"end\", ref, alt,
          row_number() OVER (PARTITION BY chrom, start, \"end\") numv
          FROM analysis.called_genotypes_titan_cluster
          WHERE aliquot_barcode = ? AND variant_type = 'SNP';"
  }
  
  rs <- dbSendQuery(con, q)
  
  vcf_list <- lapply (seq_along(aliquot_barcodes), function (index)
  {
    aliquot_barcode <- aliquot_barcodes[index]
    
    if(verbose)
      message("Fetching ", aliquot_barcode)
    
    dbBind(rs, aliquot_barcode)
    
    qres <- dbFetch(rs)
    
    if(any(qres$numv > 1))
      qres <- subset(qres, qres$numv == 1)
    
    vcf <- rowRanges(VCF(rowRanges = GRanges(seqnames = trimws(qres$chrom),
                                             ranges = IRanges(start = as.integer(qres$start),
                                                              end = qres$end),
                                             seqinfo = seqinfo(ref_genome)),
                         fixed = DataFrame(REF = DNAStringSet(qres$ref),
                                           ALT = DNAStringSetList(sapply(qres$alt,DNAStringSet)))))
    
    #rownames(vcf) <- sprintf("%s:%s_%s/%s", qres$chrom, as.integer(qres$start), qres$ref, qres$alt)
    
    if (analysis[1] == "analysis.called_genotypes_titan_cluster")
    {
      if(any(is.na(qres$aliquot_barcode_cluster)))
      {
        if (verbose)
          message("Removing ", sum(is.na(qres$aliquot_barcode_cluster)), " variants without cluster")
        
        vcf <- subset(vcf, !is.na(qres$aliquot_barcode_cluster))
        qres <- subset(qres, !is.na(qres$aliquot_barcode_cluster))
        
        stopifnot(length(vcf) == nrow(qres))
      }
      
      if (verbose)
        message("Identified ", length(unique(qres$aliquot_barcode_cluster)), " clusters")
      
      vcfout <- split(vcf, qres$aliquot_barcode_cluster)
    }
    else
    {
      vcfout <- GRangesList(vcf)
      names(vcfout) = aliquot_barcode
    }

    return(vcfout)
  })
  
  dbClearResult(rs)
  
  vcf_list <- do.call(c, vcf_list)
  
  return(vcf_list)
}

library(VariantAnnotation)
library(DBI)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg19)

mycon <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

qres <- dbGetQuery(mycon, "SELECT aliquot_barcode FROM analysis.blocklist WHERE fingerprint_exclusion = 'allow' AND coverage_exclusion = 'allow';")
aliquot_barcodes = qres$aliquot_barcode[1:5]

test1 <- mutationDbToGRanges(mycon, aliquot_barcodes, verbose = TRUE, analysis = "analysis.called_genotypes")

test2 <- mutationDbToGRanges(mycon, aliquot_barcodes, verbose = TRUE, analysis = "analysis.called_genotypes_titan_cluster")
