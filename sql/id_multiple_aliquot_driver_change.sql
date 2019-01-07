/*
Using driver stable/change annotation, identify patients with:
- Many available samples
- Driver change
*/
WITH
t1 AS
(
	SELECT case_barcode,aliquot_analysis_type,COUNT(*) AS num_aliquots
	FROM biospecimen.aliquots al
	LEFT JOIN biospecimen.samples sa ON al.sample_barcode = sa.sample_barcode
	GROUP BY case_barcode,aliquot_analysis_type
	ORDER BY 3 DESC
)
SELECT ds.case_barcode,driver_count,driver_status,target,num_aliquots
FROM analysis.driver_status ds
LEFT JOIN t1 ON ds.case_barcode = t1.case_barcode
ORDER BY 2 DESC