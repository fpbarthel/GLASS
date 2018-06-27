#!/bin/bash
## Select variants from af-only-gnomad.raw.sites.b37.vcf.gz

module load bcftools

## Need to split multi-allelic sites across multiple lines, see
## https://gatkforums.broadinstitute.org/gatk/discussion/10975/use-select-variants-on-a-gnomad-vcf-for-mutect2-contamination-filtering
bcftools norm \
	-f "/projects/verhaak-lab/verhaak_ref/gatk-legacy-bundles/b37/human_g1k_v37_decoy.fasta" \
	-m- \
	-o "/projects/verhaak-lab/verhaak_ref/gatk-cnv/af-only-gnomad.raw.sites.b37.norm.vcf.gz" \
	-O z \
	--threads 6 \
	"/projects/verhaak-lab/verhaak_ref/gatk-legacy-bundles/Mutect2/af-only-gnomad.raw.sites.b37.vcf.gz"

## Had to remove contigs from VCF file because bcftools does not add length attribute to contigs
## and SelectVariants complains if they are missing
bcftools view \
	-h "/projects/verhaak-lab/verhaak_ref/gatk-cnv/af-only-gnomad.raw.sites.b37.norm.vcf.gz" | \
	sed '/^##contig/d' \
	> "/projects/verhaak-lab/verhaak_ref/gatk-cnv/newheader.txt"

bcftools reheader \
	-h "/projects/verhaak-lab/verhaak_ref/gatk-cnv/newheader.txt" \
	-o "/projects/verhaak-lab/verhaak_ref/gatk-cnv/af-only-gnomad.raw.sites.b37.norm.reheader.vcf.gz" \
	"/projects/verhaak-lab/verhaak_ref/gatk-cnv/af-only-gnomad.raw.sites.b37.norm.vcf.gz"

## Need to have an index because SelectVariants complains w/o index
gatk IndexFeatureFile \
	-F "/projects/verhaak-lab/verhaak_ref/gatk-cnv/af-only-gnomad.raw.sites.b37.norm.reheader.vcf.gz"

## Finally we can select variants
gatk SelectVariants \
	-V  "/projects/verhaak-lab/verhaak_ref/gatk-cnv/af-only-gnomad.raw.sites.b37.norm.reheader.vcf.gz" \
	-O "/projects/verhaak-lab/verhaak_ref/gatk-cnv/af-only-gnomad.raw.sites.b37.selected.vcf.gz" \
	-R "/projects/verhaak-lab/verhaak_ref/gatk-legacy-bundles/b37/human_g1k_v37_decoy.fasta" \
	--select "AF>0.05"