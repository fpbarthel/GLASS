library(DBI)
library(odbc)
library(MutationalPatterns)
library(BSgenome)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

#Select samples to include in analysis from the database

q = "SELECT pb.pair_barcode \
FROM analysis.pairs pb \
LEFT JOIN analysis.blocklist bl ON bl.aliquot_barcode = pb.tumor_barcode \
WHERE bl.fingerprint_exclusion = 'allow' AND bl.coverage_exclusion = 'allow'"

pair_inclusion <- dbGetQuery(con,q)[,"pair_barcode"]

#Select all available VCFs that should be included in the analysis
VCF_dir <- "/projects/verhaak-lab/GLASS-analysis/results/mutect2/final/"
VCF_files <- dir(VCF_dir)
VCF_files <- VCF_files[grep("-WXS.final.vcf$",VCF_files)]		#Subset on WXS here
VCF_samples <- gsub(".final.vcf","",VCF_files)
inclusion_index <- which(VCF_samples %in% pair_inclusion)
VCF_files <- VCF_files[inclusion_index]
VCF_samples <- VCF_samples[inclusion_index]
VCF_files <- paste(VCF_dir,VCF_files,sep="")

#Load reference genome using BSgenome
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"		#Ref genome for VCF: human_g1k_v37_decoy.fasta
library(ref_genome, character.only=TRUE)

#Run mutationalSignatures on one of the files
vcfs <- read_vcfs_as_granges(VCF_files, VCF_samples, ref_genome)