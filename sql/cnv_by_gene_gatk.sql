/*
Perform gene level copy number calling using GATK CallSegments segmentation files
--
1. Define a list of genes on which to call copy number
2. Intersect gene start/stop with segments
3. Group results by gene and sample, and compute weighted average of copy number
4. Perform simple calling using -1 and 1 thresholds
*/
WITH
selected_genes AS
(
	SELECT * FROM ref.genes
),
gene_seg_intersect AS
(
    SELECT aliquot_barcode, gene_symbol, gs.chrom, upper(t0.pos * gs.pos) - lower(t0.pos * gs.pos) AS weight, log2_copy_ratio::decimal
    FROM analysis.gatk_seg gs
    INNER JOIN selected_genes t0 ON t0.chrom = gs.chrom AND t0.pos && gs.pos
),
gene_sample_call AS
(
    SELECT aliquot_barcode, gene_symbol, 
		sum(weight * log2_copy_ratio) / sum(weight) AS corrected_cn
    FROM gene_seg_intersect
    GROUP BY aliquot_barcode, gene_symbol
)
SELECT
	*,
	(CASE WHEN corrected_cn < -1 THEN -1 WHEN corrected_cn > 1 THEN 1 ELSE 0 END)::smallint cn_call
FROM gene_sample_call
--WHERE aliquot_barcode LIKE 'GLSS-SM-R109-%'
--AND gene_symbol IN ('PTEN','CDKN2A')
--ORDER BY corrected_cn ASC 