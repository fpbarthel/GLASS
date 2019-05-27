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
filtered_neoag AS
(
	SELECT aliquot_barcode,variant_id,COUNT(*)
	FROM analysis.neoantigens_by_aliquot nag
	--WHERE ssm2_pass_call IS TRUE --AND peptide_length IN (9,10)
	GROUP BY 1,2
),
variant_context_counts AS
(	
	SELECT
		pg.tumor_pair_barcode,
		(CASE
         WHEN pg.mutect2_call_a AND pg.mutect2_call_b THEN 'S'
         WHEN pg.mutect2_call_a AND NOT pg.mutect2_call_b THEN 'P'
         WHEN pg.mutect2_call_b AND NOT pg.mutect2_call_a THEN 'R'
         ELSE NULL::text
         END) AS fraction,
		(CASE WHEN pa.variant_classification IN ('MISSENSE','NONSENSE','START_CODON_SNP','NONSTOP') THEN 'non' WHEN pa.variant_classification = 'SILENT' THEN 'syn' ELSE NULL END) AS variant_class,
		(CASE WHEN nag.count > 0 THEN 'imm' WHEN nag.count IS NULL OR nag.count <= 0 THEN 'non' ELSE NULL END) AS immune_fraction,
		trinucleotide_context, 
		pa.alt,
		COUNT(*) AS mut_n
	FROM variants.pgeno pg
	INNER JOIN selected_tumor_pairs stp ON stp.tumor_pair_barcode = pg.tumor_pair_barcode
	INNER JOIN variants.passanno pa ON pa.variant_id = pg.variant_id
	LEFT JOIN filtered_neoag nag ON nag.aliquot_barcode = pg.tumor_barcode_a AND nag.variant_id = pg.variant_id
	WHERE 
		pg.variant_type = 'SNP' AND 
		(
			(pg.mutect2_call_a AND NOT pg.mutect2_call_b AND (pg.ref_count_a + pg.alt_count_a) > 14) OR 
			(pg.mutect2_call_b AND NOT pg.mutect2_call_a AND (pg.ref_count_b + pg.alt_count_b) > 14) OR 
			(pg.mutect2_call_a AND pg.mutect2_call_b AND (pg.ref_count_a + pg.alt_count_a) > 14 AND (pg.ref_count_b + pg.alt_count_b) > 14)
		) AND
		pa.variant_classification IN ('MISSENSE','NONSENSE','SILENT','START_CODON_SNP','NONSTOP')
	GROUP BY 1,2,3,4,5,6
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
	SELECT tumor_pair_barcode, fraction, trinucleotide_context, alt, SUM(mut_n) n
	FROM variant_context_counts
	WHERE variant_class = 'syn'
	GROUP BY 1,2,3,4
),
npred AS 
(
	SELECT tumor_pair_barcode, fraction, SUM(nbar * sil.n) AS npred
	FROM sil
	FULL JOIN nbar ON nbar.trinucleotide_context = sil.trinucleotide_context AND nbar.alt = sil.alt
	GROUP BY 1,2
),
bpred AS 
(
	SELECT tumor_pair_barcode, fraction, SUM(nbar * sil.n * bbar) AS bpred
	FROM sil
	FULL JOIN nbar ON nbar.trinucleotide_context = sil.trinucleotide_context AND nbar.alt = sil.alt
	FULL JOIN bbar ON bbar.trinucleotide_context = sil.trinucleotide_context AND bbar.alt = sil.alt
	GROUP BY 1,2
),
obs AS
(
	SELECT tumor_pair_barcode, fraction, SUM(CASE WHEN variant_class = 'non' AND immune_fraction = 'imm' THEN mut_n ELSE 0 END) AS bobs, SUM(CASE WHEN variant_class = 'non' THEN mut_n ELSE 0 END) AS nobs
	FROM variant_context_counts 
	GROUP BY 1,2
)
SELECT obs.tumor_pair_barcode, obs.fraction, bobs, nobs, npred, bpred, (bobs/NULLIF(nobs,0)) AS obs, (bpred/NULLIF(npred,0)) AS exp, (NULLIF(bobs,0)/NULLIF(nobs,0))/(NULLIF(bpred,0)/NULLIF(npred,0)) AS rneo
FROM obs
INNER JOIN bpred ON bpred.tumor_pair_barcode = obs.tumor_pair_barcode AND bpred.fraction = obs.fraction
INNER JOIN npred ON npred.tumor_pair_barcode = obs.tumor_pair_barcode AND npred.fraction = obs.fraction
ORDER BY 9 DESC

-- END --