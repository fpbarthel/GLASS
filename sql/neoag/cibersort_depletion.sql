/*
Creates query that includes each aliquot's observed to expected neoantigen ratio (from analysis.neoantigen_depletion) to their CIBERSORT values from the Wang et al Cancer Cell paper (PMID: 28697342)
- Results from this query are directly used to make Extended Data Figure 12C
- Collapses the CIBERSORT values from 22 cells into 11 based on lineage (annotated below)
	- Collapsing was similar to PRECOG paper (PMID: 26193342)
*/
WITH deplete AS
(
	SELECT gs.* , nd1.rneo AS nd_a, nd2.rneo AS nd_b, clin.idh_codel_subtype AS subtype
	FROM analysis.gold_set gs
	LEFT JOIN analysis.neoantigen_depletion nd1 ON gs.tumor_barcode_a = nd1.aliquot_barcode
	LEFT JOIN analysis.neoantigen_depletion nd2 ON gs.tumor_barcode_b = nd2.aliquot_barcode
	LEFT JOIN clinical.subtypes clin ON clin.case_barcode = gs.case_barcode
	WHERE nd1.rneo IS NOT NULL AND nd2.rneo IS NOT NULL AND (nd1.nobs >= 3 AND nd2.nobs >= 3) 
	ORDER BY nd1.rneo
),
bcells AS
(
	SELECT sample_barcode, sum(value) AS proportion
	FROM analysis.cibersort cs
	WHERE cs.cell IN ('Bcellsnaive', 'Bcellsmemory')
	GROUP BY sample_barcode
),
plasmacells AS
(
	SELECT sample_barcode, sum(value) AS proportion
	FROM analysis.cibersort cs
	WHERE cs.cell = 'Plasmacells'
	GROUP BY sample_barcode
),
cd8tcells AS
(
	SELECT sample_barcode, sum(value) AS proportion
	FROM analysis.cibersort cs
	WHERE cs.cell = 'TcellsCD8'
	GROUP BY sample_barcode
),
cd4tcells AS
(
	SELECT sample_barcode, sum(value) AS proportion
	FROM analysis.cibersort cs
	WHERE cs.cell IN ('TcellsCD4naive','TcellsCD4memoryresting','TcellsCD4memoryactivated','Tcellsfollicularhelper','TcellsregulatoryTregs')
	GROUP BY sample_barcode
),
gdtcells AS
(
	SELECT sample_barcode, sum(value) AS proportion
	FROM analysis.cibersort cs
	WHERE cs.cell = 'Tcellsgammadelta'
	GROUP BY sample_barcode
),
nkcells AS
(
	SELECT sample_barcode, sum(value) AS proportion
	FROM analysis.cibersort cs
	WHERE cs.cell IN ('NKcellsresting','NKcellsactivated')
	GROUP BY sample_barcode
),
macrophages AS
(
	SELECT sample_barcode, sum(value) AS proportion
	FROM analysis.cibersort cs
	WHERE cs.cell IN ('Monocytes','MacrophagesM0','MacrophagesM1','MacrophagesM2')
	GROUP BY sample_barcode
),
dendriticcells AS
(
	SELECT sample_barcode, sum(value) AS proportion
	FROM analysis.cibersort cs
	WHERE cs.cell IN ('Dendriticcellsresting','Dendriticcellsactivated')
	GROUP BY sample_barcode
),
mastcells AS
(
	SELECT sample_barcode, sum(value) AS proportion
	FROM analysis.cibersort cs
	WHERE cs.cell IN ('Mastcellsresting','Mastcellsactivated')
	GROUP BY sample_barcode
),
eosinophils AS
(
	SELECT sample_barcode, sum(value) AS proportion
	FROM analysis.cibersort cs
	WHERE cs.cell = 'Eosinophils'
	GROUP BY sample_barcode
),
pmns AS
(
	SELECT sample_barcode, sum(value) AS proportion
	FROM analysis.cibersort cs
	WHERE cs.cell = 'Neutrophils'
	GROUP BY sample_barcode
)
SELECT nd.tumor_barcode_a AS aliquot_barcode, nd_a AS nd,
	ttn1.normal_contamination AS norm_contam, 
	CASE WHEN mf1.coverage_adj_mut_freq >= 10 THEN 1 WHEN mf1.coverage_adj_mut_freq < 10 THEN 0 END AS hm,
	'Initial' AS timepoint,
	bc1.proportion AS bcells, 
	pc1.proportion AS plasmacells, 
	t81.proportion AS cd8tcells, 
	t41.proportion AS cd4tcells, 
	gt1.proportion AS gdtcells, 
	nk1.proportion AS nkcells, 
	mac1.proportion AS macrophage, 
	dc1.proportion AS dendriticcells, 
	mc1.proportion AS mastcells, 
	eo1.proportion AS eosinophils, 
	pmn1.proportion AS pmn
FROM deplete nd
INNER JOIN biospecimen.aliquots al1 ON al1.aliquot_barcode = nd.tumor_barcode_a
LEFT JOIN analysis.mut_freq mf1 ON mf1.aliquot_barcode = nd.tumor_barcode_a
LEFT JOIN variants.titan_params ttn1 ON ttn1.aliquot_barcode = nd.tumor_barcode_a
INNER JOIN bcells bc1 ON bc1.sample_barcode = al1.sample_barcode
INNER JOIN plasmacells pc1 ON pc1.sample_barcode = al1.sample_barcode
INNER JOIN cd8tcells t81 ON t81.sample_barcode = al1.sample_barcode
INNER JOIN cd4tcells t41 ON t41.sample_barcode = al1.sample_barcode
INNER JOIN gdtcells gt1 ON gt1.sample_barcode = al1.sample_barcode
INNER JOIN nkcells nk1 ON nk1.sample_barcode = al1.sample_barcode
INNER JOIN macrophages mac1 ON mac1.sample_barcode = al1.sample_barcode
INNER JOIN dendriticcells dc1 ON dc1.sample_barcode = al1.sample_barcode
INNER JOIN mastcells mc1 ON mc1.sample_barcode = al1.sample_barcode
INNER JOIN eosinophils eo1 ON eo1.sample_barcode = al1.sample_barcode
INNER JOIN pmns pmn1 ON pmn1.sample_barcode = al1.sample_barcode

UNION

SELECT nd.tumor_barcode_b AS aliquot_barcode, nd_b AS nd,
	ttn2.normal_contamination AS norm_contam,
	CASE WHEN mf2.coverage_adj_mut_freq >= 10 THEN 1 WHEN mf2.coverage_adj_mut_freq < 10 THEN 0 END AS hm,
	'Recurrent' AS timepoint,
	bc2.proportion AS bcells, 
	pc2.proportion AS plasmacells, 
	t82.proportion AS cd8tcells, 
	t42.proportion AS cd4tcells, 
	gt2.proportion AS gdtcells, 
	nk2.proportion AS nkcells, 
	mac2.proportion AS macrophage, 
	dc2.proportion AS dendriticcells, 
	mc2.proportion AS mastcells, 
	eo2.proportion AS eosinophils, 
	pmn2.proportion AS pmn
FROM deplete nd
INNER JOIN biospecimen.aliquots al2 ON al2.aliquot_barcode = nd.tumor_barcode_b
LEFT JOIN analysis.mut_freq mf2 ON mf2.aliquot_barcode = nd.tumor_barcode_b
LEFT JOIN variants.titan_params ttn2 ON ttn2.aliquot_barcode = nd.tumor_barcode_b
INNER JOIN bcells bc2 ON bc2.sample_barcode = al2.sample_barcode
INNER JOIN plasmacells pc2 ON pc2.sample_barcode = al2.sample_barcode
INNER JOIN cd8tcells t82 ON t82.sample_barcode = al2.sample_barcode
INNER JOIN cd4tcells t42 ON t42.sample_barcode = al2.sample_barcode
INNER JOIN gdtcells gt2 ON gt2.sample_barcode = al2.sample_barcode
INNER JOIN nkcells nk2 ON nk2.sample_barcode = al2.sample_barcode
INNER JOIN macrophages mac2 ON mac2.sample_barcode = al2.sample_barcode
INNER JOIN dendriticcells dc2 ON dc2.sample_barcode = al2.sample_barcode
INNER JOIN mastcells mc2 ON mc2.sample_barcode = al2.sample_barcode
INNER JOIN eosinophils eo2 ON eo2.sample_barcode = al2.sample_barcode
INNER JOIN pmns pmn2 ON pmn2.sample_barcode = al2.sample_barcode
