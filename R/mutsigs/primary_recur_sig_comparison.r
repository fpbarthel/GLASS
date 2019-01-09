library(DBI)
library(odbc)
library(MutationalPatterns)
library(BSgenome)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)

#Relative contribution differences between primary and recurrent samples
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

#Get primary-recurrent pair information
q = "SELECT * \
	FROM analysis.tumor_mut_comparison_anno tmc \
	LEFT JOIN analysis.blocklist bl ON  bl.aliquot_barcode = tmc.tumor_barcode_a \
	WHERE tmc.sample_type_a = 'TP'"
	
tp_rec_pairs <- dbGetQuery(con,q)

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
#subset_num <- 100
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
colnames(individual_fit[["contribution"]]) <- sapply(strsplit(colnames(individual_fit[["contribution"]]),"-"),function(x)paste(x[1:5],collapse="-"))
colnames(individual_fit[["reconstructed"]]) <- sapply(strsplit(colnames(individual_fit[["reconstructed"]]),"-"),function(x)paste(x[1:5],collapse="-"))
#Calculate relative contribution
contribution <- individual_fit[["contribution"]]
contribution_sums <- apply(contribution,2,sum)
rel_contribution <- apply(contribution,1,function(x)x/contribution_sums)

#Examine differences in relative signature contributions between primary and recurrent samples
tp_rec_pairs <- tp_rec_pairs[which(tp_rec_pairs[,"comparison_type"]=="longitudinal"),]
tp_rec_pairs <- tp_rec_pairs[grep("-WXS-",tp_rec_pairs[,"tumor_barcode_a"]),]
tp_rec_pairs <- tp_rec_pairs[which(tp_rec_pairs[,"sample_type_b"]=="R1"),]
shortened_barcode_a <- sapply(strsplit(tp_rec_pairs[,"tumor_barcode_a"],"-"),
					   function(x)paste(c(x[1:4],gsub("D","",x[5])),collapse="-"))
shortened_barcode_b <- sapply(strsplit(tp_rec_pairs[,"tumor_barcode_b"],"-"),
					   function(x)paste(c(x[1:4],gsub("D","",x[5])),collapse="-"))
tp_rec_pairs <- cbind(tp_rec_pairs,shortened_barcode_a,shortened_barcode_b)
primaries <- as.character(tp_rec_pairs[which(tp_rec_pairs[,"sample_type_a"] == "TP"), "shortened_barcode_a"])
recurrence <- as.character(tp_rec_pairs[which(tp_rec_pairs[,"sample_type_a"] == "TP"), "shortened_barcode_b"])

#mutational signatures pair matrix
tumor_barcode_a = tumor_barcode_b = rep("",nrow(tp_rec_pairs))
paired_contributions <- matrix(0,nrow=nrow(tp_rec_pairs),ncol=60)
colnames(paired_contributions) <- c(paste(colnames(rel_contribution),"_a",sep=""),paste(colnames(rel_contribution),"_b",sep=""))
for(i in 1:nrow(tp_rec_pairs))
{
	tumor_barcode_a[i] <- as.character(tp_rec_pairs[i,"shortened_barcode_a"])
	tumor_barcode_b[i] <- as.character(tp_rec_pairs[i,"shortened_barcode_b"])
	
	mutsigs_a <- if(tumor_barcode_a[i] %in% rownames(rel_contribution)){
		rel_contribution[tumor_barcode_a[i],]}else{
		rep(NA,30)}
	mutsigs_b <- if(tumor_barcode_b[i] %in% rownames(rel_contribution)){
		rel_contribution[tumor_barcode_b[i],]}else{
		rep(NA,30)}
	paired_contributions[i,] <- c(mutsigs_a,mutsigs_b)
}
rel_contribution_pairs <- data.frame(tumor_barcode_a,tumor_barcode_b,paired_contributions)

#add clinical info to the matrix (note: this plot uses recurrent tumor info)
sample_barcode <- sapply(strsplit(as.character(rel_contribution_pairs[,"tumor_barcode_b"]),"-"),function(x)paste(x[1:4],collapse="-"))
rel_contribution_pairs <- cbind(rel_contribution_pairs, sample_barcode)
rel_contribution_pairs <- merge(rel_contribution_pairs, tumor_info, by="sample_barcode")

#Remove everything without signatures
rel_contribution_pairs <- rel_contribution_pairs[-which(is.na(rel_contribution_pairs[,"Signature.1_a"])),]
rel_contribution_pairs <- rel_contribution_pairs[-which(is.na(rel_contribution_pairs[,"Signature.1_b"])),]

#Comparison 1: naive primary vs recurrent differences
whole_dataset_dif <- rel_contribution_pairs[,34:63] - rel_contribution_pairs[,4:33]
whole_dataset_means <- apply(whole_dataset_dif,2,mean)
whole_dataset_sds <- apply(whole_dataset_dif,2,sd)

signatures <- gsub("_b","",gsub("Signature.","",names(whole_dataset_means)))
plot_whole_dataset <- data.frame(whole_dataset_means,whole_dataset_sds,signatures)
plot_whole_dataset[,"signatures"] <- factor(plot_whole_dataset[,"signatures"],levels = signatures)

#Comparison 2: idh wt noncodels
idhwt_noncodels <- rel_contribution_pairs[which(rel_contribution_pairs[,"idh_codel_subtype"]=="IDHwt_noncodel"),]
idhwt_noncodels_dif <- idhwt_noncodels[,34:63] - idhwt_noncodels[,4:33]
idhwt_noncodels_means <- apply(idhwt_noncodels_dif,2,mean)
idhwt_noncodels_sds <- apply(idhwt_noncodels_dif,2,sd)

signatures <- gsub("_b","",gsub("Signature.","",names(idhwt_noncodels_means)))
plot_idhwt_noncodels <- data.frame(idhwt_noncodels_means,idhwt_noncodels_sds,signatures)
plot_idhwt_noncodels[,"signatures"] <- factor(plot_idhwt_noncodels[,"signatures"],levels = signatures)


#Comparison 3: idh mut noncodels
idhmut_noncodels <- rel_contribution_pairs[which(rel_contribution_pairs[,"idh_codel_subtype"]=="IDHmut_noncodel"),]
idhmut_noncodels_dif <- idhmut_noncodels[,34:63] - idhmut_noncodels[,4:33]
idhmut_noncodels_means <- apply(idhmut_noncodels_dif,2,mean)
idhmut_noncodels_sds <- apply(idhmut_noncodels_dif,2,sd)

signatures <- gsub("_b","",gsub("Signature.","",names(idhmut_noncodels_means)))
plot_idhmut_noncodels <- data.frame(idhmut_noncodels_means,idhmut_noncodels_sds,signatures)
plot_idhmut_noncodels[,"signatures"] <- factor(plot_idhmut_noncodels[,"signatures"],levels = signatures)


#Comparison 4: idh mut codels
idhmut_codels <- rel_contribution_pairs[which(rel_contribution_pairs[,"idh_codel_subtype"]=="IDHmut_codel"),]
idhmut_codels_dif <- idhmut_codels[,34:63] - idhmut_codels[,4:33]
idhmut_codels_means <- apply(idhmut_codels_dif,2,mean)
idhmut_codels_sds <- apply(idhmut_codels_dif,2,sd)

signatures <- gsub("_b","",gsub("Signature.","",names(idhmut_codels_means)))
plot_idhmut_codels <- data.frame(idhmut_codels_means,idhmut_codels_sds,signatures)
plot_idhmut_codels[,"signatures"] <- factor(plot_idhmut_codels[,"signatures"],levels = signatures)

#Faceted barplot
analysis_type <- c(rep("Whole dataset",nrow(plot_whole_dataset)),
				 rep("IDHwt noncodel",nrow(plot_idhwt_noncodels)),
				 rep("IDHmut noncodel",nrow(plot_idhmut_noncodels)),
				 rep("IDHmut codel",nrow(plot_idhmut_codels)))
colnames(plot_whole_dataset) = colnames(plot_idhwt_noncodels) = 
colnames(plot_idhmut_noncodels) = colnames(plot_idhmut_codels) = c("means","sds","Signatures")
plot_results <- rbind(plot_whole_dataset, plot_idhwt_noncodels, plot_idhmut_noncodels, plot_idhmut_codels)
plot_results <- cbind(plot_results,analysis_type)
plot_results[,"analysis_type"] <- factor(plot_results[,"analysis_type"],levels = c("Whole dataset","IDHwt noncodel","IDHmut noncodel","IDHmut codel"))

pdf("/projects/varnf/GLASS/Figures/signatures/test/pri_recur_sig_differences.pdf")
ggplot(data=plot_results, aes(x=Signatures,y=means,width=1)) +
geom_bar(stat="identity",position="identity",size=0.2)+
geom_errorbar(aes(ymin=means-sds, ymax=means+sds), width=.2,
                 position=position_dodge(.9)) +
facet_grid(analysis_type ~ .) +
ylab("Relative contribution") + guides(fill=FALSE) +
theme_bw()+
theme(axis.title.y = element_text(size = 12, 
	vjust = 1), axis.text.y = element_text(size = 10), 
	axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 10, 
	vjust = 0.4), strip.text.x = element_text(size = 12), 
    strip.text.y = element_text(size = 12), panel.grid.major.x = element_blank(), 
    panel.spacing.x = unit(0, "lines"))
dev.off()

#Signature 16 significance test in codels:
wilcox.test(idhmut_codels[,"Signature.16_a"],idhmut_codels[,"Signature.16_b"],paired=TRUE)		#0.08

#Cluster differences
whole_dataset_mat <- rel_contribution_pairs[,34:63] - rel_contribution_pairs[,4:33]
rownames(whole_dataset_mat) <- paste(rel_contribution_pairs[,"tumor_barcode_a"],rel_contribution_pairs[,"tumor_barcode_b"],sep="-")
colnames(whole_dataset_mat) <- gsub("_b","",colnames(whole_dataset_mat))
annotation_table <- rel_contribution_pairs[,c("grade","idh_codel_subtype","surgical_interval_mo","prior_tmz","prior_radiation")]

#Set colors for sidebars
myset = brewer.pal(9,"Set1")
sidebar_labels <- apply(annotation_table,2,function(x)unique(x[!is.na(x)]))
#sample_type <- c("black","white"); names(sample_type) = sidebar_labels[[1]]
grade = myset[1:3]; names(grade) = sidebar_labels[[1]]
idh_codel_subtype = myset[4:6]; names(idh_codel_subtype) = sidebar_labels[[2]]
prior_tmz = c("black","white"); names(prior_tmz) = sidebar_labels[[4]]
prior_radiation = c("black","white"); names(prior_radiation) = sidebar_labels[[5]]
annotation_colors = list(grade=grade, idh_codel_subtype=idh_codel_subtype,
					prior_tmz=prior_tmz, prior_radiation=prior_radiation)
ha = HeatmapAnnotation(df = annotation_table,which="row",
	 col=annotation_colors)
	 
pdf("/projects/varnf/GLASS/Figures/signatures/sig_differences_complex_heatmap.pdf",width=7,height=5)
Heatmap(whole_dataset_mat,
		col = rev(brewer.pal(11,"RdBu")),
		cluster_rows=TRUE,
		show_row_names=FALSE,
		) +
ha
dev.off()