library(DBI)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ICAMS)

# define functions for SplitOneMutectVCF and SplitListofMutectVCFs
SplitOneMutectVCF <- function(vcf.df) {
  # Mutect VCFs can represent multiple non-reference alleles at the
  # same site; the alleles are separated by commas in the ALT columm;
  # these are quite rare and often dubious, so we ignore them.
  multiple.alt <- grep(",", vcf.df$ALT, fixed = TRUE)
  
  multiple.alt.df <- vcf.df[multiple.alt, ]
  
  if(length(multiple.alt) > 0)
    df <- vcf.df[-multiple.alt, ]
  else
    df <- vcf.df
  
  rm(multiple.alt, vcf.df)
  
  SBS.df <- df[nchar(df$REF) == 1 & nchar(df$ALT) == 1, ]
  
  DBS.df <- df[nchar(df$REF) == 2 & nchar(df$ALT) == 2, ]
  
  other.df <- df[nchar(df$REF) > 2 & nchar(df$ALT) == nchar(df$REF), ]
  
  ID.df <- df[nchar(df$REF) != nchar(df$ALT), ]
  
  return(list(SBS = SBS.df, DBS = DBS.df, ID = ID.df,
              other=other.df, multiple.alt = multiple.alt.df))
  
}

#' Split each Mutect VCF into SBS, DBS, and ID VCFs (plus two
#' VCF-like data frame with left-over rows).
#'
#' @param list.of.vcfs List of VCFs as in-memory data.frames.
#'
#' @return A list with 3 in-memory VCFs and two left-over
#' VCF-like data frames with rows that were not incorporated
#' into the first 3 VCFs, as follows:
#'
#' \enumerate{
#'
#'  \item \code{SBS} VCF with only single base substitutions.
#'
#'  \item \code{DBS} VCF with only doublet base substitutions
#'   as called by Mutect.
#'
#'  \item \code{ID} VCF with only small insertions and deletions.
#'
#'  \item \code{other.subs} VCF like data.frame with
#'  rows for coordinate substitutions involving
#'  3 or more nucleotides, e.g. ACT > TGA or AACT > GGTA.
#'
#'  \item \code{multiple.alternative.alleles} VCF like data.frame with
#'  rows for variants with multiple alternative alleles, for example
#'  ACT mutated to both AGT and ACT at the same position.
#'
#' }
#'
#' @keywords internal
SplitListOfMutectVCFs <- function(list.of.vcfs) {
  v1 <- lapply(list.of.vcfs, SplitOneMutectVCF)
  SBS <- lapply(v1, function(x) x$SBS)
  DBS <- lapply(v1, function(x) x$DBS)
  ID  <- lapply(v1, function(x) x$ID)
  other.subs <- lapply(v1, function(x) x$other.df)
  multiple.alternative.alleles <-
    lapply(v1, function(x) x$multiple.alt)
  
  return(list(SBS = SBS, DBS = DBS, ID = ID,
              other.subs = other.subs,
              multiple.alternative.alleles
              = multiple.alternative.alleles
  ))
}

#load in data from database
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")
time <- system.time(qres <- dbGetQuery(con, read_file(('/Users/c-kocake/Box Sync/Projects/GLASS-Hypermutation/sql/passanno_join_passgeno.sql'))))

m2df <- data.frame(CHROM = ifelse(qres$chrom == 23, "X", as.character(qres$chrom)),
                   POS = qres$start_pos,
                   ID = ".",
                   REF = qres$ref,
                   ALT = qres$alt,
                   QUAL = ".",
                   FILTER = "PASS",
                   INFO = ".",
                   FORMAT = ".",
                   VAF = qres$af,
                   stringsAsFactors = FALSE)

m2dfl = split(m2df, qres$aliquot_barcode)

split.vcfs <- SplitListOfMutectVCFs(m2dfl)

ref.genome <- BSgenome.Hsapiens.1000genomes.hs37d5
trans.ranges <- NULL #trans.ranges.GRCh37
region <- "unknown"

catalogs <- c(VCFsToSBSCatalogs(split.vcfs$SBS, ref.genome, trans.ranges, region),
              VCFsToDBSCatalogs(split.vcfs$DBS, ref.genome, trans.ranges, region),
              list(catID = VCFsToIDCatalogs(split.vcfs$ID, ref.genome, region)))

PlotCatalogToPdf.SBS96Catalog(catalogs, file = "/Users/c-kocake/Box Sync/Projects/GLASS-Hypermutation/mutational_signatures_v3/test123.pdf", plot.SBS12 = NULL, cex = NULL,
                 grid = NULL, upper = NULL, xlabels = NULL)

# define catalogs
SBS <- as.data.frame(catalogs[2]) %>% rename_all(list(name = ~ (str_replace(., "catSBS1536.", ""))))
DBS <- as.data.frame(catalogs[3]) %>% rename_all(list(name = ~ (str_replace(., "catDBS78.", ""))))
ID <- as.data.frame(catalogs[5]) %>% rename_all(list(name = ~ (str_replace(., "catID.", ""))))

PAN <- rbind(SBS, DBS, ID)

save.image('/Users/c-kocake/Box Sync/Projects/GLASS-Hypermutation/mutational_signatures_v3/data/signatures_input')


