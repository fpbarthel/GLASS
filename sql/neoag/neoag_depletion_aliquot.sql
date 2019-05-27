/*
-----------------------------------------------------------------------------------
Neoantigen depletion
PostgreSQL implementation of Rooney et all (Cell 2015) neoantigen depletion method
See https://www.ncbi.nlm.nih.gov/pubmed/25594174
Authors: Fred Varn, Floris barthel
-----------------------------------------------------------------------------------

## TERMS ##

Nbar: the expected number of non-silent mutations per silent mutation
Bbar: the expected number of high-affinity neo-peptide binders per non-silent mutation

** Here we compute Nbar and Bbar by summing all samples in the gold set **

Npred: the predicted number of non-silent mutations in sample s
Bpred: the predicted number of neo-peptide binders in sample s
Nobs: the observed number of non-silent mutations in sample s
Bobs: the observed number of neoepitope-generating SNVs in sample s

** Npred and Bpred are calculated using the observed silent mutations in sample s, Bbar and Nbar **

Rneo: the ratio between the observed and expected rate of neo-peptides

- Since the observed values and Nbar/Bbar were defined based on the same dataset, this ratio follows a normal distribution and is centered at one
- Lower values of this score can be interpreted as evidence of higher neoantigen depletion relative to other samples in the dataset

*/
WITH
selected_tumor_pairs AS
(
	SELECT * FROM analysis.gold_set
),
selected_aliquots AS
(
	SELECT tumor_barcode_a AS aliquot_barcode, case_barcode, 'P' AS sample_type FROM selected_tumor_pairs
	UNION
	SELECT tumor_barcode_b AS aliquot_barcode, case_barcode, 'R' AS sample_type  FROM selected_tumor_pairs
),
filtered_neoag AS
(
	SELECT aliquot_barcode,variant_id,COUNT(*)
	FROM analysis.neoantigens_by_aliquot nag
	WHERE ssm2_pass_call IS TRUE --AND peptide_length IN (9,10)
	GROUP BY 1,2
),
variant_context_counts AS
(	
	SELECT
		pg.aliquot_barcode,
		(CASE WHEN variant_classification IN ('MISSENSE','NONSENSE','START_CODON_SNP','NONSTOP') THEN 'non' WHEN variant_classification = 'SILENT' THEN 'syn' ELSE NULL END) AS variant_class,
		(CASE WHEN nag.count > 0 THEN 'imm' WHEN nag.count IS NULL OR nag.count <= 0 THEN 'non' ELSE NULL END) AS immune_fraction,
		trinucleotide_context, 
		pa.alt,
		COUNT(*) AS mut_n
	FROM variants.passgeno pg
	INNER JOIN selected_aliquots sa ON sa.aliquot_barcode = pg.aliquot_barcode
	INNER JOIN variants.passanno pa ON pa.variant_id = pg.variant_id
	LEFT JOIN filtered_neoag nag ON nag.aliquot_barcode = pg.aliquot_barcode AND nag.variant_id = pg.variant_id
	WHERE pg.ssm2_pass_call IS TRUE AND variant_type = 'SNP' AND pg.ad_alt + pg.ad_ref >= 15 AND variant_classification IN ('MISSENSE','NONSENSE','SILENT','START_CODON_SNP','NONSTOP')
	GROUP BY 1,2,3,4,5
),
mis AS
(
	SELECT trinucleotide_context, alt, SUM(mut_n) n
	FROM variant_context_counts
	WHERE variant_class = 'non'
	GROUP BY 1,2
),
syn AS
(
	SELECT trinucleotide_context, alt, SUM(mut_n) n
	FROM variant_context_counts
	WHERE variant_class = 'syn'
	GROUP BY 1,2
),
imm AS
(
	SELECT trinucleotide_context, alt, SUM(mut_n) n
	FROM variant_context_counts
	WHERE variant_class = 'non' AND immune_fraction = 'imm'
	GROUP BY 1,2
),
nbar AS
(
	SELECT mis.trinucleotide_context, mis.alt, mis.n / syn.n AS nbar
	FROM mis
	FULL JOIN syn ON syn.trinucleotide_context = mis.trinucleotide_context AND syn.alt = mis.alt
),
bbar AS
(
	SELECT imm.trinucleotide_context, imm.alt, imm.n / mis.n AS bbar
	FROM imm
	FULL JOIN mis ON mis.trinucleotide_context = imm.trinucleotide_context AND mis.alt = imm.alt
),
sil AS
(
	SELECT aliquot_barcode, trinucleotide_context, alt, SUM(mut_n) n
	FROM variant_context_counts
	WHERE variant_class = 'syn'
	GROUP BY 1,2,3
),
npred AS 
(
	SELECT aliquot_barcode, SUM(nbar * sil.n) AS npred
	FROM sil
	FULL JOIN nbar ON nbar.trinucleotide_context = sil.trinucleotide_context AND nbar.alt = sil.alt
	GROUP BY 1
),
bpred AS 
(
	SELECT aliquot_barcode, SUM(nbar * sil.n * bbar) AS bpred
	FROM sil
	FULL JOIN nbar ON nbar.trinucleotide_context = sil.trinucleotide_context AND nbar.alt = sil.alt
	FULL JOIN bbar ON bbar.trinucleotide_context = sil.trinucleotide_context AND bbar.alt = sil.alt
	GROUP BY 1
),
obs AS
(
	SELECT aliquot_barcode, SUM(CASE WHEN variant_class = 'non' AND immune_fraction = 'imm' THEN mut_n ELSE 0 END) AS bobs, SUM(CASE WHEN variant_class = 'non' THEN mut_n ELSE 0 END) AS nobs
	FROM variant_context_counts 
	GROUP BY 1
)
SELECT obs.aliquot_barcode, bobs, nobs, npred, bpred, (bobs/nobs) AS obs, (bpred/npred) AS exp, (bobs/nobs)/(bpred/npred) AS rneo
FROM obs
INNER JOIN bpred ON bpred.aliquot_barcode = obs.aliquot_barcode
INNER JOIN npred ON npred.aliquot_barcode = obs.aliquot_barcode
ORDER BY 8 DESC

-- END --