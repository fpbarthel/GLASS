#Generate mutation matrix for all primary/recurred/shared mutations across all WXS/WGS samples

#Helper function to get VCFs from databaset (get most recent data)
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
        WHERE (mutect2_call_a OR mutect2_call_b) AND read_depth_a >= 15 AND read_depth_b >= 15 AND snvs.variant_type = 'SNP';"
  
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

vcfs <- privateVsSharedMutationDbToGRanges(con,  verbose = TRUE)

#Run mutationalSignatures on the files
#--------------------------------------------------
#Load reference genome using BSgenome
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"		#Ref genome for VCF: human_g1k_v37_decoy.fasta
library(ref_genome, character.only=TRUE)

#Create a mutation matrix 
#--------------------------------------------------
mut_mat <- mut_matrix(vcfs, ref_genome = ref_genome)
write.table(mut_mat, "/projects/varnf/GLASS/analysis/signatures/GLASS_primary_recur_shared_mut_matrix.txt",sep="\t",quote=F,row.names=F)
save(mut_mat, file= "/projects/varnf/GLASS/analysis/signatures/GLASS_primary_recur_shared_mut_matrix.rda")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#Analysis
library(BSgenome.Hsapiens.UCSC.hg19)
library(MutationalPatterns)
library(ComplexHeatmap)
library(RColorBrewer)
library(DBI)
library(reshape)

rm(list=ls())
myinf1 <- "/projects/varnf/GLASS/analysis/signatures/GLASS_primary_recur_shared_mut_matrix.rda"

load(myinf1)		#Loads as mut_mat

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

#shared
individual_fit <- fit_to_signatures(mut_mat, cancer_signatures)
individual_contribution <- individual_fit[["contribution"]]
contribution_sums <- apply(individual_contribution,2,sum)
individual_contribution <- apply(individual_contribution,1,function(x)x/contribution_sums)

shared_contribution <- individual_contribution[grep("-shared",rownames(individual_contribution)),]
primary_contribution <- individual_contribution[grep("-primary",rownames(individual_contribution)),]
recur_contribution <- individual_contribution[grep("-recurrent",rownames(individual_contribution)),]


#WXS heatmap
wxs_contribution <- individual_contribution[grep("-WXS-",rownames(individual_contribution)),]

#Heatmap of all samples
annotation <- data.frame(sapply(strsplit(rownames(wxs_contribution),"-"),function(x)x[9]))
colnames(annotation) <- "sample_type"

#Set colors for sidebars
sample_type_col <- rep("",nrow(annotation))
sample_type_col[which(annotation[,"sample_type"] == "primary")] <- "red"
sample_type_col[which(annotation[,"sample_type"] == "recurrent")] <- "blue"
sample_type_col[which(annotation[,"sample_type"] == "shared")] <- "green"
names(sample_type_col) <- annotation[,"sample_type"]
annotation_colors = list(sample_type=sample_type_col)
ha = HeatmapAnnotation(df = annotation,which="row",
	 col=annotation_colors)
	 
pdf("/projects/varnf/GLASS/Figures/signatures/primary_recur_shared_wxs_heatmap.pdf",width=7,height=5)
Heatmap(wxs_contribution,
		col = brewer.pal(9,"YlGnBu"),
		cluster_rows=TRUE,
		show_row_names=FALSE,
		) +
ha
dev.off()

#WGS heatmap
wgs_contribution <- individual_contribution[grep("-WGS-",rownames(individual_contribution)),]

#Heatmap of all samples
annotation <- data.frame(sapply(strsplit(rownames(wgs_contribution),"-"),function(x)x[9]))
colnames(annotation) <- "sample_type"

#Set colors for sidebars
sample_type_col <- rep("",nrow(annotation))
sample_type_col[which(annotation[,"sample_type"] == "primary")] <- "red"
sample_type_col[which(annotation[,"sample_type"] == "recurrent")] <- "blue"
sample_type_col[which(annotation[,"sample_type"] == "shared")] <- "green"
names(sample_type_col) <- annotation[,"sample_type"]
annotation_colors = list(sample_type=sample_type_col)
ha = HeatmapAnnotation(df = annotation,which="row",
	 col=annotation_colors)
	 
pdf("/projects/varnf/GLASS/Figures/signatures/primary_recur_shared_wgs_heatmap.pdf",width=7,height=5)
Heatmap(wgs_contribution,
		col = brewer.pal(9,"YlGnBu"),
		cluster_rows=TRUE,
		show_row_names=FALSE,
		) +
ha
dev.off()

#--------------------------------------------------

#Integrate clinical information
#Get annotations
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

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

#----------------------------

#names for each wxs sample
tumor_info <- tumor_info[-which(is.na(tumor_info[,"sample_barcode"])),]
rownames(tumor_info) <- tumor_info[,"sample_barcode"]
wxs_sample_names <- sapply(strsplit(rownames(wxs_contribution),"-"),function(x)paste(x[1:4],collapse="-"))

idh_codel_subtype <- tumor_info[wxs_sample_names, "idh_codel_subtype"]
names(idh_codel_subtype) <- wxs_sample_names
tmz_status <- tumor_info[wxs_sample_names, "prior_tmz"]
names(tmz_status) <- wxs_sample_names


wt_noncodel <- wxs_contribution[which(idh_codel_subtype=="IDHwt_noncodel"),]
wt_noncodel_p <- apply(wt_noncodel[grep("-primary",rownames(wt_noncodel)),],2,mean)
wt_noncodel_r <- apply(wt_noncodel[grep("-recurrent",rownames(wt_noncodel)),],2,mean)
wt_noncodel_s <- apply(wt_noncodel[grep("-shared",rownames(wt_noncodel)),],2,mean)

mut_noncodel <- wxs_contribution[which(idh_codel_subtype=="IDHmut_noncodel"),]
mut_noncodel_p <- apply(mut_noncodel[grep("-primary",rownames(mut_noncodel)),],2,mean)
mut_noncodel_r <- apply(mut_noncodel[grep("-recurrent",rownames(mut_noncodel)),],2,mean)
mut_noncodel_s <- apply(mut_noncodel[grep("-shared",rownames(mut_noncodel)),],2,mean)

mut_codel <- wxs_contribution[which(idh_codel_subtype=="IDHmut_codel"),]
mut_codel_p <- apply(mut_codel[grep("-primary",rownames(mut_codel)),],2,mean)
mut_codel_r <- apply(mut_codel[grep("-recurrent",rownames(mut_codel)),],2,mean)
mut_codel_s <- apply(mut_codel[grep("-shared",rownames(mut_codel)),],2,mean)

average_contributions <- rbind(wt_noncodel_p, wt_noncodel_r, wt_noncodel_s,
							   mut_noncodel_p, mut_noncodel_r, mut_noncodel_s,
							   mut_codel_p, mut_codel_r, mut_codel_s)

plot_contributions <- melt(average_contributions)
colnames(plot_contributions) <- c("ID","Signature","rel_contribution")
plot_contributions[,"tumor_status"] <- rep(c("primary","recurrent","shared"),nrow(plot_contributions)/3)
plot_contributions[,"subtype"] <- sapply(strsplit(as.character(plot_contributions[,1]),"_"),function(x)paste("IDH",x[1],x[2],sep="_"))


#Plot 1: Stacked barplots for all samples
#Plot contribution barplot
pdf("/projects/varnf/GLASS/Figures/signatures/primary_recur_shared_wxs_barplots.pdf",width=10,height=10)
ggplot(plot_contributions, aes(y=rel_contribution, x=subtype, fill=tumor_status)) +
geom_bar(stat="identity",position="dodge", colour="black") +
facet_wrap(Signature ~ .) +
theme_bw()+
theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), 
	panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),
	axis.title.y = element_text(size = 12, vjust = 1), axis.text.y = element_text(size = 10), 
	axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 10, 
	vjust = 0.4), strip.text.x = element_text(size = 12), 
    strip.text.y = element_text(size = 12), 
    panel.spacing.x = unit(0, "lines"))
dev.off()




#Plot 2: Stacked barplots for all samples that did not receive TMZ

