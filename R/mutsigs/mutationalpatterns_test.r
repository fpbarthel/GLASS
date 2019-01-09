library(DBI)
library(odbc)
library(MutationalPatterns)
library(BSgenome)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

#Data cleanup
#--------------------------------------------------
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

subset_num <- 10
VCF_files <- VCF_files[1:subset_num]
VCF_samples <- VCF_samples[1:subset_num]

primaries <- VCF_samples[grep("-TP-",VCF_samples)]
R1s <- VCF_samples[grep("-R1-",VCF_samples)]
R2s <- VCF_samples[grep("-R2-",VCF_samples)]
R3s <- VCF_samples[grep("-R3-",VCF_samples)]
R4s <- VCF_samples[grep("-R4-",VCF_samples)]
recurrences <- c(R1s,R2s,R3s,R4s)


#Run mutationalSignatures on ten of the files
#--------------------------------------------------
#Load reference genome using BSgenome
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"		#Ref genome for VCF: human_g1k_v37_decoy.fasta
library(ref_genome, character.only=TRUE)

#Read VCFs as granges
vcfs <- read_vcfs_as_granges(VCF_files, VCF_samples, ref_genome)

#Create a mutation matrix consisting of the average number of mutations in primary and recurrent samples
#--------------------------------------------------
mut_mat <- mut_matrix(vcfs, ref_genome = ref_genome)
primary_average <- apply(mut_mat[,primaries],1,mean)
recurrence_average <- apply(mut_mat[,recurrences],1,mean)
collapsed_mut_mat <- matrix(c(primary_average,recurrence_average),
					 nrow=nrow(mut_mat),ncol=2,
					 dimnames=list(rownames(mut_mat),c("Primary","Recurrent")))

#Primary vs recurrence 96 mutational profile
pdf("/projects/varnf/GLASS/Figures/signatures/test/mut96_profile_TPvR.pdf",width=7,height=3)
plot_96_profile(collapsed_mut_mat,condensed=TRUE)
dev.off()

#Optimal contribution of known signatures (COSMIC)
#--------------------------------------------------
#Download mutational signatures from the COSMIC website
sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/", 
		 "signatures_probabilities.txt", sep="")
cancer_signatures <- read.table(sp_url, sep ="\t", header=TRUE)
#Match the order of the mutational types to MutationalPatterns standard
new_order <- match(row.names(collapsed_mut_mat), cancer_signatures[,"Somatic.Mutation.Type"])
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

#Similarity between mutational profiles and Cosmic signatures- heatmap
individual_sample_sim_matrix <- cos_sim_matrix(mut_mat,cancer_signatures)
rownames(individual_sample_sim_matrix) <- sapply(strsplit(rownames(individual_sample_sim_matrix),"-"),function(x)paste(x[1:5],collapse="-"))
write.table(individual_sample_sim_matrix, 
"/projects/varnf/GLASS/analysis/signatures/test/individual_sample_similarity_matrix.txt",
sep="\t",quote=F)

pdf("/projects/varnf/GLASS/Figures/signatures/test/individual_signature_heatmap.pdf",width=7,height=5)
plot_cosine_heatmap(individual_sample_sim_matrix, col_order=cosmic_order, cluster_rows=TRUE)
dev.off()

primary_vs_recurrent_sim_matrix <- cos_sim_matrix(collapsed_mut_mat,cancer_signatures)
write.table(primary_vs_recurrent_sim_matrix, 
"/projects/varnf/GLASS/analysis/signatures/test/primary_vs_recurrent_similarity_matrix.txt",
sep="\t",quote=F)

pdf("/projects/varnf/GLASS/Figures/signatures/test/primary_vs_recurrent_signature_heatmap.pdf",width=7,height=5)
plot_cosine_heatmap(primary_vs_recurrent_sim_matrix, col_order=cosmic_order, cluster_rows=TRUE)
dev.off()

#Optimal contribution of COSMIC signatures to 96 mutational profiles

#individual
individual_fit <- fit_to_signatures(mut_mat, cancer_signatures)
colnames(individual_fit[["contribution"]]) <- sapply(strsplit(colnames(individual_fit[["contribution"]]),"-"),function(x)paste(x[1:5],collapse="-"))
colnames(individual_fit[["reconstructed"]]) <- sapply(strsplit(colnames(individual_fit[["reconstructed"]]),"-"),function(x)paste(x[1:5],collapse="-"))
#Select signatures with some contribution
select <- which(rowSums(individual_fit[["contribution"]])>10)

#Plot contribution barplot
pdf("/projects/varnf/GLASS/Figures/signatures/test/individual_relative_sig_contribution.pdf",width=7,height=5)
plot_contribution(individual_fit[["contribution"]][select,],
					cancer_signatures[,select],
					coord_flip=FALSE,
					mode="relative")
dev.off()
pdf("/projects/varnf/GLASS/Figures/signatures/test/individual_absolute_relative_sig_contribution.pdf",width=7,height=5)
plot_contribution(individual_fit[["contribution"]][select,],
					cancer_signatures[,select],
					coord_flip=FALSE,
					mode="absolute")
dev.off()
#Plot contribution heatmap
pdf("/projects/varnf/GLASS/Figures/signatures/test/individual_sig_contribution_heatmap.pdf",width=7,height=5)
plot_contribution_heatmap(individual_fit[["contribution"]],
						  cluster_samples=TRUE,
						  method="complete")
dev.off()

#grouped primary and recurrent
primary_vs_recurrent_fit <- fit_to_signatures(collapsed_mut_mat, cancer_signatures)
#Select signatures with some contribution
select <- which(rowSums(primary_vs_recurrent_fit[["contribution"]])>10)

#Plot contribution barplot
pdf("/projects/varnf/GLASS/Figures/signatures/test/primary_vs_recurrent_relative_sig_contribution.pdf",width=4,height=5)
plot_contribution(primary_vs_recurrent_fit[["contribution"]][select,],
					cancer_signatures[,select],
					coord_flip=FALSE,
					mode="relative")
dev.off()
pdf("/projects/varnf/GLASS/Figures/signatures/test/primary_vs_recurrent_absolute_relative_sig_contribution.pdf",width=4,height=5)
plot_contribution(primary_vs_recurrent_fit[["contribution"]][select,],
					cancer_signatures[,select],
					coord_flip=FALSE,
					mode="absolute")
dev.off()
#Plot contribution heatmap
pdf("/projects/varnf/GLASS/Figures/signatures/test/primary_vs_recurrent_sig_contribution_heatmap.pdf",width=7,height=5)
plot_contribution_heatmap(primary_vs_recurrent_fit[["contribution"]],
						  cluster_samples=TRUE,
						  method="complete")
dev.off()
