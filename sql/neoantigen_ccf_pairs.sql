WITH pairs AS
(
	SELECT tp.tumor_pair_barcode,
		tp.case_barcode,
		tp.tumor_barcode_a,
		tp.tumor_barcode_b,
		tp.sample_type_a,
		tp.sample_type_b,
		tp.portion_a,
		tp.portion_b,
		tp.comparison_type,
		tp.surgical_interval_mo,
		n1.variant_id,
		n1.chrom,
		n1.pos,
		n1.alt,
		n1.transcript,
		n1.ensembl_gene_id,
		n1.pvacseq_variant_type,
		n1.mutation,
		n1.pvacseq_protein_position,
		n1.gene_name,
		n1.hla_allele,
		n1.peptide_length,
		n1.sub_peptide_position,
		n1.mutation_position,
		n1.mt_epitope_seq,
		n1.wt_epitope_seq,
		n1.netmhcpan_mt_score,
		n1.netmhcpan_wt_score,
		n1.netmhcpan_fold_change,
		n1.ad_ref AS ref_count_a,
		n2.ad_ref AS ref_count_b,
		n1.ad_alt AS alt_count_a,
		n2.ad_alt AS alt_count_b,
		n1.af AS af_a,
		n2.af AS af_b,
		n1.ssm2_pass_call AS mutect2_call_a,
		n2.ssm2_pass_call AS mutect2_call_b,
		n1.pyclone_ccf AS pyclone_ccf_a,
		n2.pyclone_ccf AS pyclone_ccf_b
	FROM analysis.neoantigens_by_aliquot n1
		JOIN analysis.neoantigens_by_aliquot n2 ON n1.case_barcode = n2.case_barcode AND n1.chrom = n2.chrom AND n1.pos = n2.pos AND n1.alt = n2.alt AND n1.hla_allele = n2.hla_allele AND n1.peptide_length = n2.peptide_length AND n1.sub_peptide_position = n2.sub_peptide_position
		JOIN analysis.tumor_pairs tp ON tp.tumor_barcode_a = n1.aliquot_barcode AND tp.tumor_barcode_b = n2.aliquot_barcode
)
SELECT tumor_pair_barcode,
(CASE WHEN (mutect2_call_a AND mutect2_call_b) THEN 'S' WHEN (mutect2_call_a AND NOT mutect2_call_b) THEN 'P' WHEN (mutect2_call_b AND NOT mutect2_call_a) THEN 'R' END) AS fraction,
chrom, pos, alt, transcript, ensembl_gene_id, pvacseq_variant_type, mutation, pvacseq_protein_position, gene_name, hla_allele, peptide_length, sub_peptide_position, mutation_position, mt_epitope_seq, wt_epitope_seq, netmhcpan_mt_score,netmhcpan_wt_score,netmhcpan_fold_change,
pyclone_ccf_a,pyclone_ccf_b
FROM pairs
WHERE ((mutect2_call_a  AND NOT mutect2_call_b AND (ref_count_a + alt_count_a) > 14) OR  (mutect2_call_b AND NOT mutect2_call_a AND (ref_count_b + alt_count_b) > 14) OR
(mutect2_call_a AND mutect2_call_b AND (ref_count_a + alt_count_a) > 14 AND (ref_count_b + alt_count_b) > 14))
AND (pyclone_ccf_a IS NOT NULL AND pyclone_ccf_b IS NOT NULL)
ORDER BY tumor_pair_barcode,fraction,chrom,pos

