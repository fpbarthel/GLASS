library(VariantAnnotation)
library(stringr)

vcff = snakemake@input[["vcf"]]
tsvf = snakemake@output[["tsv"]]

vcf <- readVcf(vcff)

funcolumns <- unlist(strsplit(unlist(strsplit(info(header(vcf))['FUNCOTATION',3], '\\: '))[2],'\\|'))
funcotation <- as.data.frame(do.call('rbind', str_split(gsub("^\\[|\\]$","",as.character(info(vcf)[,'FUNCOTATION'])), "\\|")))
colnames(funcotation) <- funcolumns

df <- data.frame(chrom = gsub("^chr","",as.character(seqnames(vcf))),
                 pos = sprintf("[%s,%s]", start(vcf), end(vcf)),
                 ref = ref(vcf),
                 alt = unstrsplit(CharacterList(alt(vcf)), sep=","),
                 gene_symbol = funcotation$Gencode_19_hugoSymbol,
                 variant_classification = funcotation$Gencode_19_variantClassification,
                 secondary_variant_classification = funcotation$Gencode_19_secondaryVariantClassification,
                 variant_type = funcotation$Gencode_19_variantType,
                 genome_change = funcotation$Gencode_19_genomeChange,
                 transcript = funcotation$Gencode_19_annotationTranscript,
                 transcript_strand = funcotation$Gencode_19_transcriptStrand,
                 transcript_exon = funcotation$Gencode_19_transcriptExon,
                 transcript_position = funcotation$Gencode_19_transcriptPos,
                 cdna_change = funcotation$Gencode_19_cDnaChange,
                 cds_change = funcotation$Gencode_19_codonChange,
                 protein_change = funcotation$Gencode_19_proteinChange,
                 gc_content = funcotation$Gencode_19_gcContent, 
                 reference_context = funcotation$Gencode_19_referenceContext,
                 stringsAsFactors = FALSE)

write.table(df, file = tsvf, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
