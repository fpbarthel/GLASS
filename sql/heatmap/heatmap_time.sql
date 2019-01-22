WITH t1 AS
(
	SELECT
		case_barcode,
		surgery_number::varchar(255),
		surgical_interval_mo AS time_mo,
		idh_codel_subtype
	FROM clinical.surgeries
	WHERE idh_codel_subtype IS NOT NULL

	UNION

	SELECT
		cc.case_barcode,
		case_vital_status::varchar(255) AS surgery_number,
		case_overall_survival_mo AS time_mo,
		idh_codel_subtype
	FROM clinical.cases cc
	LEFT JOIN clinical.subtypes cs ON cc.case_barcode = cs.case_barcode
)
SELECT t1.*, case_source
FROM t1
INNER JOIN analysis.silver_set ss ON t1.case_barcode = ss.case_barcode
LEFT JOIN clinical.cases cc ON t1.case_barcode = cc.case_barcode