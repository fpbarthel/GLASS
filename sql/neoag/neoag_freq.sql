/*
Create materialized view that gives the mutation count and the neoantigen-causing mutation count (analysis.neoag_freq)
- Uses the analysis.neoantigens_by_aliquot table for neoantigen counts and the variants.passgeno table for mutation counts
- Counts each mutation giving rise to a neoantigen in a patient only once (some mutations give rise to several putative neoantigens in a patient)
- Both neoantigen and mutation calls are filtered based on whether the ssm2_pass_call IS TRUE
- Filters mutations so that only those with HIGH and MODERATE variant classification impacts are included (trying to filter for coding mutations)
	- Currently this query compares neoantigens (by definition coding) as defined using VEP to coding mutations as defined using Funcotator
	- *** This needs to be modified once we have VEP calls for every mutation ***
- Coverage filters exist in case we want to also calculate coverage-adjusted mutation and neoantigen frequency (currently commented out)
- Calculates the proportion immunogenic (neoantigen count divided by mutation count) for Figure 4A
- Figure 4A is subset on the gold set before being made (no subsetting occurs here)
*/

WITH neoag_table AS 
(
	SELECT neoantigens_by_aliquot.aliquot_barcode,
		neoantigens_by_aliquot.variant_id,
		row_number() OVER (PARTITION BY neoantigens_by_aliquot.aliquot_barcode, neoantigens_by_aliquot.variant_id ORDER BY neoantigens_by_aliquot.netmhcpan_mt_score) AS neoag_number
	FROM analysis.neoantigens_by_aliquot
	WHERE neoantigens_by_aliquot.ssm2_pass_call IS TRUE
), 
neoag_count AS 
(
	SELECT neoag_table.aliquot_barcode,
		count(*) AS neoag_count
	FROM neoag_table
	WHERE neoag_table.neoag_number = 1
	GROUP BY neoag_table.aliquot_barcode
), 
mt_count AS 
(
	SELECT pg.aliquot_barcode,
		count(*) AS mt_count
	FROM variants.passgeno pg
	JOIN variants.passanno pa ON pa.variant_id = pg.variant_id
	JOIN variants.variant_classifications vc ON vc.variant_classification::text = pa.variant_classification::text
	WHERE pg.ssm2_pass_call IS TRUE AND (vc.variant_classification_impact::text = ANY (ARRAY['HIGH'::character varying::text, 'MODERATE'::character varying::text]))
	GROUP BY pg.aliquot_barcode
)
SELECT cov.aliquot_barcode,
	cov.cumulative_coverage,
	fmc.mt_count,
	fnc.neoag_count,
	--COALESCE(round(fmc.mt_count::numeric / cov.cumulative_coverage::numeric * '1000000'::numeric, 4), 0::numeric) AS coverage_adj_mut_freq,
	--COALESCE(round(fnc.neoag_count::numeric / cov.cumulative_coverage::numeric * '1000000'::numeric, 4), 0::numeric) AS coverage_adj_neoag_freq,
	fnc.neoag_count::numeric / fmc.mt_count::numeric AS prop_immunogenic
FROM analysis.coverage cov
LEFT JOIN analysis.blocklist bl ON bl.aliquot_barcode = cov.aliquot_barcode
LEFT JOIN neoag_count fnc ON fnc.aliquot_barcode = cov.aliquot_barcode
LEFT JOIN mt_count fmc ON fmc.aliquot_barcode = cov.aliquot_barcode
WHERE cov.coverage = 15 AND bl.fingerprint_exclusion = 'allow'::bpchar AND bl.coverage_exclusion = 'allow'::bpchar AND fmc.mt_count IS NOT NULL
ORDER BY fnc.aliquot_barcode