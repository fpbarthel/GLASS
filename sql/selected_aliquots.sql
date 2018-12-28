/*
Select aliquots
- Stringent blocklist filtering (diamond set)
	* Fingerprinting
	* Coverage
	* CNV
- Create sample short names using surgery number and portion
- Drop cases with less than 2 aliquots
*/
WITH
selected_aliquots AS
(
	SELECT case_barcode, aliquot_analysis_type, al.aliquot_barcode, ROUND(purity::decimal,2) AS purity, case_barcode || '-' || aliquot_analysis_type AS short_name, COUNT(*) OVER (PARTITION BY case_barcode, aliquot_analysis_type) AS num_samples
	FROM biospecimen.aliquots al
	LEFT JOIN analysis.blocklist bl ON bl.aliquot_barcode = al.aliquot_barcode
	LEFT JOIN clinical.surgeries su ON al.sample_barcode = su.sample_barcode
	LEFT JOIN analysis.pairs pa ON al.aliquot_barcode = pa.tumor_barcode
	LEFT JOIN analysis.titan_params tp ON tp.pair_barcode = pa.pair_barcode
	WHERE 
		bl.fingerprint_exclusion = 'allow' AND
		bl.coverage_exclusion = 'allow' AND
		bl.cnv_exclusion = 'allow'
	ORDER BY 1, su.surgery_number, al.aliquot_portion
)
SELECT * FROM selected_aliquots WHERE num_samples > 1