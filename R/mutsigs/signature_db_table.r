library(DBI)
library(odbc)
library(MutationalPatterns)
library(BSgenome)
library(ComplexHeatmap)
library(RColorBrewer)
library(reshape)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

#Data cleanup
#--------------------------------------------------
#Select samples to include in analysis from the database
q = "SELECT pb.pair_barcode \
FROM analysis.pairs pb \
LEFT JOIN analysis.blocklist bl ON bl.aliquot_barcode = pb.tumor_barcode \
WHERE bl.fingerprint_exclusion = 'allow' AND bl.coverage_exclusion = 'allow'"

pair_inclusion <- dbGetQuery(con,q)[,"pair_barcode"]

#Get annotations
q = "SELECT su.*, \
    CASE WHEN su.surgery_number > 1 \
    THEN \
        (SELECT treatment_tmz FROM clinical.surgeries WHERE case_barcode = su.case_barcode AND surgery_number = su.surgery_number - 1) \
    ELSE FALSE \
    END AS prior_tmz, \
    CASE WHEN su.surgery_number > 1 \
    THEN \
        (SELECT treatment_radiotherapy FROM clinical.surgeries WHERE case_barcode = su.case_barcode AND surgery_number = su.surgery_number - 1) \
    ELSE FALSE \
    END AS prior_radiation \
FROM clinical.surgeries su \
ORDER BY case_barcode, surgery_number"

tumor_info <- dbGetQuery(con,q)

#Select all available VCFs that should be included in the analysis
VCF_dir <- "/projects/verhaak-lab/GLASS-analysis/results/mutect2/final/"
VCF_files <- dir(VCF_dir)
VCF_files <- VCF_files[grep("-WXS.final.vcf$",VCF_files)]		#Subset on WXS here
VCF_samples <- gsub(".final.vcf","",VCF_files)
inclusion_index <- which(VCF_samples %in% pair_inclusion)
VCF_files <- VCF_files[inclusion_index]
VCF_samples <- VCF_samples[inclusion_index]
VCF_files <- paste(VCF_dir,VCF_files,sep="")

#Testing code with 10 samples
#subset_num <- 10
#VCF_files <- VCF_files[1:subset_num]
#VCF_samples <- VCF_samples[1:subset_num]


#Run mutationalSignatures on the files
#--------------------------------------------------
#Load reference genome using BSgenome
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"		#Ref genome for VCF: human_g1k_v37_decoy.fasta
library(ref_genome, character.only=TRUE)

#Read VCFs as granges
vcfs <- read_vcfs_as_granges(VCF_files, VCF_samples, ref_genome)

#Create a mutation matrix for each sample
#--------------------------------------------------
mut_mat <- mut_matrix(vcfs, ref_genome = ref_genome)

#Optimal contribution of known signatures (COSMIC)
#--------------------------------------------------
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

#individual
individual_fit <- fit_to_signatures(mut_mat, cancer_signatures)

#Make database compatible signature table for WXS:
db_table <- individual_fit[["contribution"]]
db_table <- apply(db_table,1,function(x)x/contribution_sums)
db_table <- melt(db_table)
db_table <- db_table[,c(2,1,3)]
colnames(db_table) <- c("signature","barcode","value")

write.table(db_table,"/projects/varnf/GLASS/analysis/signatures/mutSig_WXS_rel_contribution.txt",sep="\t",quote=F,row.names=F)