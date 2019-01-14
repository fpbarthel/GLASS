library(DBI)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg19)
library(VariantAnnotation)
library(MutationalPatterns)

setwd("/projects/verhaak-lab/GLASS-analysis")

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

genome = "BSgenome.Hsapiens.UCSC.hg19"

ref_genome <- base::get(genome)
ref_organism <- GenomeInfoDb::organism(ref_genome)
ref_style <- seqlevelsStyle(ref_genome)

genome_name <- genome(ref_genome)[[1]]
seqlevelsStyle(ref_genome) = "NCBI"

if (!(class(ref_genome) == "BSgenome"))
  stop("Please provide the name of a BSgenome object.")

message("Fetching VCF data from database...")
time <- system.time(qres <- dbGetQuery(con, read_file("sql/mutsig_input_by_fraction.sql")))

message("Fetch completed in ", round(time[3]), " seconds.")

message("Building VCF object and extracting ranges.")

vcf <- rowRanges(VCF(rowRanges = GRanges(seqnames = trimws(qres$chrom),
                                         ranges = IRanges(start = as.integer(qres$start_pos),
                                                          end = as.integer(qres$end_pos)),
                                         seqinfo = seqinfo(ref_genome),
                                         paramRangeID = rep(factor(NA),nrow(qres))),
                     fixed = DataFrame(REF = DNAStringSet(qres$ref),
                                       ALT = unname(split(DNAStringSet(qres$alt),1:length(qres$alt))),
                                       QUAL = as.numeric(NA_integer_),
                                       FILTER = 'PASS')))

vcf$fraction <- qres$fraction
vcf$tumor_pair_barcode <- qres$tumor_pair_barcode

seqlevelsStyle(vcf) = "UCSC"

groups <- c(extractSeqlevelsByGroup(species = ref_organism,
                                    style = ref_style,
                                    group = "auto"),
            extractSeqlevelsByGroup(species = ref_organism,
                                    style = ref_style,
                                    group = "sex"))

groups <- intersect(groups, seqlevels(vcf))
vcf <- keepSeqlevels(vcf, groups)

message("Splitting VCF object by fraction")
vcf_by_fraction <- split(vcf, sprintf("%s-%s", vcf$tumor_pair_barcode, vcf$fraction))
vcf_by_fraction <- lapply(vcf_by_fraction, function(gr) { names(gr) <- sprintf("%s:%s_%s/%s", gsub("chr","",as.character(seqnames(gr))), start(gr), gr$REF, unlist(gr$ALT)); return(gr) })
vcf_by_fraction <- GRangesList(vcf_by_fraction)

########## BY SAMPLE #########

message("Fetching VCF data from database...")
time <- system.time(qres <- dbGetQuery(con, read_file("sql/mutsig_input_by_sample.sql")))

message("Fetch completed in ", round(time[3]), " seconds.")

message("Building VCF object and extracting ranges.")

vcf <- rowRanges(VCF(rowRanges = GRanges(seqnames = trimws(qres$chrom),
                                         ranges = IRanges(start = as.integer(qres$start_pos),
                                                          end = as.integer(qres$end_pos)),
                                         seqinfo = seqinfo(ref_genome),
                                         paramRangeID = rep(factor(NA),nrow(qres))),
                     fixed = DataFrame(REF = DNAStringSet(qres$ref),
                                       ALT = unname(split(DNAStringSet(qres$alt),1:length(qres$alt))),
                                       QUAL = as.numeric(NA_integer_),
                                       FILTER = 'PASS')))

vcf$aliquot_barcode <- qres$aliquot_barcode

seqlevelsStyle(vcf) = "UCSC"

groups <- c(extractSeqlevelsByGroup(species = ref_organism,
                                    style = ref_style,
                                    group = "auto"),
            extractSeqlevelsByGroup(species = ref_organism,
                                    style = ref_style,
                                    group = "sex"))

groups <- intersect(groups, seqlevels(vcf))
vcf <- keepSeqlevels(vcf, groups)

message("Splitting VCF object by sample.")
vcf_by_sample <- split(vcf, vcf$aliquot_barcode)
vcf_by_sample <- lapply(vcf_by_sample, function(gr) { names(gr) <- sprintf("%s:%s_%s/%s", gsub("chr","",as.character(seqnames(gr))), start(gr), gr$REF, unlist(gr$ALT)); return(gr) })
vcf_by_sample <- GRangesList(vcf_by_sample)

message("Constructing mutation matrix - fraction")
mut_mat_fraction <- mut_matrix(vcf_by_fraction, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")

message("Constructing mutation matrix - sample")
mut_mat_sample <- mut_matrix(vcf_by_sample, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")

#Optimal contribution of known signatures (COSMIC)
#--------------------------------------------------
#Download mutational signatures from the COSMIC website
sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/", 
                "signatures_probabilities.txt", sep="")
cancer_signatures <- read.table(sp_url, sep ="\t", header=TRUE)
new_order <- match(row.names(mut_mat_sample), cancer_signatures[,"Somatic.Mutation.Type"])
cancer_signatures <- cancer_signatures[as.vector(new_order),]
row.names(cancer_signatures) <- cancer_signatures[,"Somatic.Mutation.Type"]
cancer_signatures <- as.matrix(cancer_signatures[,4:33])
hclust_cosmic <- cluster_signatures(cancer_signatures, method = "average")
cosmic_order = colnames(cancer_signatures)[hclust_cosmic[["order"]]]
individual_fit <- fit_to_signatures(mut_mat_sample, cancer_signatures)
absolute_contribution <- individual_fit[["contribution"]]									#Absolute contribution
contribution_sums <- apply(absolute_contribution,2,sum)
relative_contribution <- apply(absolute_contribution,1,function(x)x/contribution_sums)		#Relative contribution
absolute_contribution <- t(absolute_contribution)

abs_mut_df_sample <- absolute_contribution %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "aliquot_barcode") %>% 
  gather(key = "signature", value = "signature_score", -aliquot_barcode) %>%
  mutate(signature = as.numeric(gsub("Signature.","",signature)))

rel_mut_df_sample <- relative_contribution %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "aliquot_barcode") %>% 
  gather(key = "signature", value = "signature_score", -aliquot_barcode) %>%
  mutate(signature = as.numeric(gsub("Signature.","",signature)))

cancer_signatures <- read.table(sp_url, sep ="\t", header=TRUE)
new_order <- match(row.names(mut_mat_fraction), cancer_signatures[,"Somatic.Mutation.Type"])
cancer_signatures <- cancer_signatures[as.vector(new_order),]
row.names(cancer_signatures) <- cancer_signatures[,"Somatic.Mutation.Type"]
cancer_signatures <- as.matrix(cancer_signatures[,4:33])
hclust_cosmic <- cluster_signatures(cancer_signatures, method = "average")
cosmic_order = colnames(cancer_signatures)[hclust_cosmic[["order"]]]
individual_fit <- fit_to_signatures(mut_mat_fraction, cancer_signatures)
absolute_contribution <- individual_fit[["contribution"]]
contribution_sums <- apply(absolute_contribution,2,sum)
relative_contribution <- apply(absolute_contribution,1,function(x)x/contribution_sums)
absolute_contribution <- t(absolute_contribution)

abs_mut_df_fraction <- absolute_contribution %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "key") %>% 
  mutate(tumor_pair_barcode = substr(key,1,29),
         fraction = substr(key,31,31)) %>%
  dplyr::select(-key) %>%
  gather(key = "signature", value = "signature_score", -tumor_pair_barcode, -fraction) %>%
  mutate(signature = as.numeric(gsub("Signature.","",signature)))
  
rel_mut_df_fraction <- relative_contribution %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "key") %>% 
  mutate(tumor_pair_barcode = substr(key,1,29),
         fraction = substr(key,31,31)) %>%
  dplyr::select(-key) %>%
  gather(key = "signature", value = "signature_score", -tumor_pair_barcode, -fraction) %>%
  mutate(signature = as.numeric(gsub("Signature.","",signature)))

abs_mut_df_fraction$mut_count <- apply(mut_mat_fraction, 2, sum)
abs_mut_df_sample$mut_count <- apply(mut_mat_sample, 2, sum)
rel_mut_df_fraction$mut_count <- apply(mut_mat_fraction, 2, sum)
rel_mut_df_sample$mut_count <- apply(mut_mat_sample, 2, sum)

dbWriteTable(con, Id(schema="analysis",table="abs_mutsig_fraction"), abs_mut_df_fraction, append = FALSE)
dbWriteTable(con, Id(schema="analysis",table="abs_mutsig_sample"), abs_mut_df_sample, append = FALSE)
dbWriteTable(con, Id(schema="analysis",table="rel_mutsig_fraction"), rel_mut_df_fraction, append = FALSE)
dbWriteTable(con, Id(schema="analysis",table="rel_mutsig_sample"), rel_mut_df_sample, append = FALSE)