SELECT
	cp.tumor_pair_barcode,
	cp.case_barcode,
	(CASE WHEN recurrence_location = 'Distal' THEN 1 WHEN recurrence_location = 'Local' THEN 0 ELSE NULL END) location_distal,
	(CASE WHEN grade_change = 'Grade up' THEN 1 WHEN grade_change IN ('Grade up', 'Grade stable') THEN 0 ELSE NULL END) grade_change,
	(CASE WHEN received_alk = '1' THEN 1 WHEN received_alk = '0' THEN 0 ELSE NULL END) received_alk,
	(CASE WHEN received_rt = '1' THEN 1 WHEN received_rt = '0' THEN 0 ELSE NULL END) received_rt,
	(CASE WHEN hypermutator_status = '1' THEN 1 WHEN hypermutator_status = '0' THEN 0 ELSE NULL END) is_hypermutator,
	idh_codel_subtype
FROM analysis.tumor_clinical_comparison cp
INNER JOIN analysis.gold_set ss ON ss.tumor_pair_barcode = cp.tumor_pair_barcode
LEFT JOIN clinical.subtypes st ON st.case_barcode = cp.case_barcode