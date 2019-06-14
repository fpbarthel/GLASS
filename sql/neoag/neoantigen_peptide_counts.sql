/*
Creates Supplementary Table 6 (?): List of all unique neoantigens in the GLASS cohort and the number of initial/recurrent tumors harboring each one
A separate R script is used to save the table for publication: R/neoag/analysis/SuppTable6_writetottext.r
*/

WITH neoag_by_ali AS
(
	SELECT aliquot_barcode, variant_id, gene_name, mutation, pvacseq_protein_position, peptide_length, sub_peptide_position, mt_epitope_seq 
	FROM analysis.neoantigens_by_aliquot neo
	WHERE ssm2_pass_call = TRUE
	GROUP BY aliquot_barcode, variant_id, gene_name, mutation, pvacseq_protein_position, peptide_length, sub_peptide_position, mt_epitope_seq
),
ini_counts AS
(
	SELECT neo.gene_name, neo.pvacseq_protein_position, neo.mutation, neo.mt_epitope_seq, COUNT(*) AS total
	FROM analysis.gold_set gs 
	LEFT JOIN neoag_by_ali neo ON neo.aliquot_barcode = gs.tumor_barcode_a
	GROUP BY neo.gene_name, neo.pvacseq_protein_position, neo.mutation, neo.mt_epitope_seq
	ORDER BY total DESC
),
rec_counts AS
(
	SELECT neo.gene_name, neo.pvacseq_protein_position, neo.mutation, neo.mt_epitope_seq, COUNT(*) AS total
	FROM analysis.gold_set gs 
	LEFT JOIN neoag_by_ali neo ON neo.aliquot_barcode = gs.tumor_barcode_b
	GROUP BY neo.gene_name, neo.pvacseq_protein_position, neo.mutation, neo.mt_epitope_seq
	ORDER BY total DESC
)
SELECT ini.gene_name, ini.pvacseq_protein_position, ini.mutation, ini.mt_epitope_seq, 
COALESCE(ini.total,0) AS initial_total, 
COALESCE(rec.total,0) AS recurrent_total, 
COALESCE(ini.total,0) + COALESCE(rec.total,0) AS total_tumors
FROM ini_counts ini
LEFT JOIN rec_counts rec ON rec.gene_name = ini.gene_name AND 
	rec.pvacseq_protein_position = ini.pvacseq_protein_position AND
	rec.mutation = ini.mutation AND
	rec.mt_epitope_seq = ini.mt_epitope_seq
ORDER BY total_tumors DESC
