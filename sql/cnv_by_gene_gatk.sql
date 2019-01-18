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
	SELECT * FROM ref.genes WHERE gene_symbol IN ('CDKN2A','EGFR')
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
	gc.aliquot_barcode,
	gc.gene_symbol,
	gc.corrected_cn,
	(CASE
	 WHEN gc.corrected_cn < -1 THEN '-'
	 WHEN gc.corrected_cn > 1 THEN '+'
	 ELSE '0' END) hard_thresholded_call,
	(CASE
	 WHEN 2^gc.corrected_cn >= 0.9 AND 2^gc.corrected_cn <= 1.1 THEN '0'
	 WHEN (2^gc.corrected_cn - fwmean) < -2.0 * fwsd THEN '-'
	 WHEN (2^gc.corrected_cn - fwmean) > 2.0 * fwsd THEN '+'
	 ELSE '0'
	 END) zscore_call_2sd,
	titan_call,
	(2^gc.corrected_cn - fwmean) / fwsd AS zscore,
	(CASE
	 WHEN 2^gc.corrected_cn >= 0.9 AND 2^gc.corrected_cn <= 1.1 THEN 0
	 WHEN 2^gc.corrected_cn < loss_wmean - 2 * loss_wsd THEN -2
	 WHEN (2^gc.corrected_cn - fwmean) < -2.0 * fwsd THEN -1
	 WHEN 2^gc.corrected_cn > gain_wmean + 2 * gain_wsd THEN 2
	 WHEN (2^gc.corrected_cn - fwmean) > 2.0 * fwsd THEN 1
	 ELSE 0
	 END) hlvl_call
FROM gene_sample_call gc
LEFT JOIN analysis.filtered_seg_wmean_wsd rs ON rs.aliquot_barcode = gc.aliquot_barcode
LEFT JOIN analysis.pairs tp ON tp.tumor_barcode = gc.aliquot_barcode
LEFT JOIN analysis.cnv_by_gene cng ON cng.pair_barcode = tp.pair_barcode AND cng.gene_symbol = gc.gene_symbol
LEFT JOIN analysis.taylor_aneuploidy_v2 ta ON ta.aliquot_barcode = gc.aliquot_barcode
WHERE fwsd > 0