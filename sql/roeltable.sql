 SELECT crosstab.case_source_description,
    crosstab.case_project,
    crosstab.aliquot_analysis_type,
    crosstab."Primary",
    crosstab."1st Recurrence",
    crosstab."2nd Recurrence",
    crosstab."3rd Recurrence",
    crosstab."4th Recurrence"
   FROM crosstab('
	SELECT case_source_description, case_project, aliquot_analysis_type, sample_type, COUNT( DISTINCT ca.case_barcode )
	FROM biospecimen.aliquots AS al
	INNER JOIN biospecimen.samples AS sa ON sa.sample_barcode = al.sample_barcode
	INNER JOIN clinical.cases AS ca ON ca.case_barcode = sa.case_barcode
	INNER JOIN clinical.case_sources AS cs ON ca.case_source = cs.case_source
	WHERE sa.sample_type IN (''TP'', ''R1'', ''R2'', ''R3'', ''R4'')
	GROUP BY case_source_description, case_project, aliquot_analysis_type, sample_type
	ORDER BY 2,3,1
	'::text, '
	SELECT sample_type FROM biospecimen.sample_types WHERE sample_type IN (''TP'', ''R1'', ''R2'', ''R3'', ''R4'')
	'::text) crosstab(case_source_description character varying, case_project character(4), aliquot_analysis_type character(3), "Primary" integer, "1st Recurrence" integer, "2nd Recurrence" integer, "3rd Recurrence" integer, "4th Recurrence" integer);