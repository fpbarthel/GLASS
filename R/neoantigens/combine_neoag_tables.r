library(DBI)
library(odbc)

rm(list=ls())
myDir1 <- "/projects/varnf/GLASS/GLASS/results/pvacseq/neoantigens/"
myintDir <- "/projects/varnf/GLASS/GLASS/results/pvacseq/tabular_results/labelled/"
myoutf <- "/projects/varnf/GLASS/GLASS/results/pvacseq/tabular_results/final/neoantigens_full.txt"

myrun <- dir(myDir1)
mysample <- sapply(strsplit(myrun,"_"),function(x)x[1])

myfile <- paste(myDir1,myrun,"/MHC_Class_I/",myrun,".final.tsv",sep="")

file.exist <- file.exists(myfile)
mysample <- mysample[which(file.exist)]
myfile <- myfile[which(file.exist)]

myintfile <- paste(myintDir,mysample,".labelled.txt",sep="")

for(i in 1:length(myfile))
{
	tumor_barcode <- sapply(strsplit(myfile[i],"/"),function(x)x[9])

	tmp.table <- read.delim(myfile[i],stringsAsFactor=FALSE,
	colClasses=c(Chromosome="character", Start="numeric",Stop="numeric",
				 Reference="character", Variant="character",
				 Transcript="character", Ensembl.Gene.ID="character",
				 Variant.Type="character",Mutation="character",
				 Protein.Position="character",Gene.Name="character",
				 HLA.Allele="character",Peptide.Length="numeric",
				 Sub.peptide.Position="numeric",Mutation.Position="numeric",
				 MT.Epitope.Seq="character",WT.Epitope.Seq="character",
				 Best.MT.Score.Method="character",Best.MT.Score="numeric",
				 Corresponding.WT.Score="numeric",Corresponding.Fold.Change="numeric",
				 Tumor.DNA.Depth="numeric",Tumor.DNA.VAF="numeric",
				 Tumor.RNA.Depth="numeric",Tumor.RNA.VAF="numeric",
				 Normal.Depth="numeric",Normal.VAF="numeric",
				 Gene.Expression="numeric",Transcript.Expression="numeric",
				 Median.MT.Score="numeric",Median.WT.Score="numeric",
				 Median.Fold.Change="numeric",NetMHCpan.WT.Score="numeric",
				 NetMHCpan.MT.Score="numeric"))
	tumor_barcode <- rep(tumor_barcode, nrow(tmp.table))
	
	tmp.table <- cbind(tmp.table,tumor_barcode)
	
	write.table(tmp.table, myintfile[i], sep="\t", quote=FALSE)
}



big_table <- do.call(rbind, lapply(myintfile, 
	function(x) read.delim(x, stringsAsFactors = FALSE,
	colClasses=c(Chromosome="character", Start="integer",Stop="integer",
	 Reference="character", Variant="character",
	 Transcript="character", Ensembl.Gene.ID="character",
	 Variant.Type="character",Mutation="character",
	 Protein.Position="character",Gene.Name="character",
	 HLA.Allele="character",Peptide.Length="integer",
	 Sub.peptide.Position="integer",Mutation.Position="integer",
	 MT.Epitope.Seq="character",WT.Epitope.Seq="character",
	 Best.MT.Score.Method="character",Best.MT.Score="numeric",
	 Corresponding.WT.Score="numeric",Corresponding.Fold.Change="numeric",
	 Tumor.DNA.Depth="numeric",Tumor.DNA.VAF="numeric",
	 Tumor.RNA.Depth="numeric",Tumor.RNA.VAF="numeric",
	 Normal.Depth="numeric",Normal.VAF="numeric",
	 Gene.Expression="numeric",Transcript.Expression="numeric",
	 Median.MT.Score="numeric",Median.WT.Score="numeric",
	 Median.Fold.Change="numeric",NetMHCpan.WT.Score="numeric",
	 NetMHCpan.MT.Score="numeric"))))

#Choose columns with information. Note: For this run Best and Median scores/fold-changes are the same because only NetMHCpan was used
big_table_edit <- big_table[,c("tumor_barcode","Chromosome","Start","Stop","Reference",
			"Variant","Transcript","Ensembl.Gene.ID","Variant.Type","Mutation","Protein.Position",
			"Gene.Name","HLA.Allele","Peptide.Length","Sub.peptide.Position","Mutation.Position",
			"MT.Epitope.Seq","WT.Epitope.Seq","NetMHCpan.MT.Score","NetMHCpan.WT.Score",
			"Corresponding.Fold.Change")]
colnames(big_table_edit) <- tolower(c("tumor_barcode","chrom","Start","Stop","ref",
			"alt","Transcript","Ensembl_Gene_ID","Variant_Type","Mutation","Protein_Position",
			"Gene_Name","HLA_Allele","Peptide_Length","Sub_peptide_Position","Mutation_Position",
			"MT_Epitope_Seq","WT_Epitope_Seq","NetMHCpan_MT_Score","NetMHCpan_WT_Score",
			"NetMHCpan_Fold_Change"))

#Change X chromosone to 23 for speed
big_table_edit[which(big_table_edit[,"chrom"]=='X'),"chrom"] = 23
big_table_edit[,"chrom"] = as.numeric(big_table_edit[,"chrom"])

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")
dbWriteTable(con, Id(schema="analysis",table="neoantigens_by_patient"), big_table_edit, overwrite=TRUE, row.names=FALSE)

write.table(big_table_edit, myoutf,sep="\t",quote=F,row.names=F)
