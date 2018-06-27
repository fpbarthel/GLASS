#!/bin/bash
## Select variants from af-only-gnomad.raw.sites.b37.vcf.gz

module load bcftools

bcftools norm \
	-f "/projects/verhaak-lab/verhaak_ref/gatk-legacy-bundles/b37/human_g1k_v37_decoy.fasta" \
	-m- \
	-o "/projects/verhaak-lab/verhaak_ref/gatk-cnv/af-only-gnomad.raw.sites.b37.norm.vcf.gz" \
	-O z \
	--threads 6 \
	"/projects/verhaak-lab/verhaak_ref/gatk-legacy-bundles/Mutect2/af-only-gnomad.raw.sites.b37.vcf.gz"

bcftools view \
	-h "/projects/verhaak-lab/verhaak_ref/gatk-cnv/af-only-gnomad.raw.sites.b37.norm.vcf.gz" | \
	sed '/^##contig/d' \
	> "/projects/verhaak-lab/verhaak_ref/gatk-cnv/newheader.txt"

bcftools reheader \
	-h "/projects/verhaak-lab/verhaak_ref/gatk-cnv/newheader.txt" \
	-o "/projects/verhaak-lab/verhaak_ref/gatk-cnv/af-only-gnomad.raw.sites.b37.norm.vcf.gz" \
	"/projects/verhaak-lab/verhaak_ref/gatk-cnv/af-only-gnomad.raw.sites.b37.norm.vcf.gz"

gatk IndexFeatureFile \
	-F "/projects/verhaak-lab/verhaak_ref/gatk-cnv/af-only-gnomad.raw.sites.b37.norm.vcf.gz"

gatk SelectVariants \
	-V  "/projects/verhaak-lab/verhaak_ref/gatk-cnv/af-only-gnomad.raw.sites.b37.norm.vcf.gz" \
	-L 1 \
	-O "/projects/verhaak-lab/verhaak_ref/gatk-cnv/af-only-gnomad.raw.sites.b37.selected.vcf.gz" \
	-R "/projects/verhaak-lab/verhaak_ref/gatk-legacy-bundles/b37/human_g1k_v37_decoy.fasta" \
	--select "AF>0.05"