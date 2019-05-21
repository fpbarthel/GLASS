/*
GLASS database version 2 clinical_tumor_pairs.
*/

WITH roman_to_int (grade, grade_int)
AS
(
     VALUES ('I', 1),
            ('II', 2),
            ('III', 3),
            ('IV', 4)
)
SELECT
    tumor_pair_barcode,
    tp.case_barcode,
    tumor_barcode_a,
    tumor_barcode_b,
    (CASE
     WHEN s1.surgery_location = s2.surgery_location AND (s1.surgery_laterality = s2.surgery_laterality OR (s1.surgery_laterality IS NULL AND s2.surgery_laterality IS NULL)) THEN 'Local'
     WHEN s1.surgery_location <> s2.surgery_location OR s1.surgery_laterality <> s2.surgery_laterality THEN 'Distal'
     END) AS recurrence_location,
    (SELECT bool_or(treatment_tmz) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) AS received_tmz,
    (SELECT sum(treatment_tmz_cycles) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number)::integer AS received_tmz_sum_cycles,
    (SELECT bool_or(treatment_radiotherapy) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) AS received_rt,
    (SELECT sum(treatment_radiation_dose_gy) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) AS received_rt_sum_gy,
    (SELECT bool_or(treatment_alkylating_agent) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) AS received_alk,
    (SELECT bool_or(treatment_alkylating_agent or treatment_radiotherapy) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) AS received_treatment,
    (CASE
     WHEN r2.grade_int > r1.grade_int THEN 'Grade up'
     WHEN r2.grade_int = r1.grade_int THEN 'Grade stable'
     WHEN r2.grade_int < r1.grade_int THEN 'Grade down'
     END) AS grade_change,
    s2.surgical_interval_mo - s1.surgical_interval_mo AS surgical_interval,
    (CASE
     WHEN mf2.coverage_adj_mut_freq > mf1.coverage_adj_mut_freq AND mf2.coverage_adj_mut_freq > 10 THEN true 
     WHEN mf2.coverage_adj_mut_freq < 10 THEN false
     ELSE NULL
     END) hypermutator_status,
    (CASE
     WHEN mf2.coverage_adj_mut_freq > mf1.coverage_adj_mut_freq AND mf2.coverage_adj_mut_freq > 10 AND (SELECT bool_or(treatment_alkylating_agent) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) IS TRUE THEN true 
     WHEN mf2.coverage_adj_mut_freq < 10 AND (SELECT bool_or(treatment_alkylating_agent) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) IS TRUE THEN false
     ELSE NULL
     END) alk_assoc_hypermutator_status
FROM analysis.tumor_pairs tp
LEFT JOIN biospecimen.aliquots a1 ON a1.aliquot_barcode = tp.tumor_barcode_a
LEFT JOIN biospecimen.aliquots a2 ON a2.aliquot_barcode = tp.tumor_barcode_b
LEFT JOIN clinical.surgeries s1 ON s1.sample_barcode = a1.sample_barcode
LEFT JOIN clinical.surgeries s2 ON s2.sample_barcode = a2.sample_barcode
LEFT JOIN roman_to_int r1 ON r1.grade = s1.grade
LEFT JOIN roman_to_int r2 ON r2.grade = s2.grade
LEFT JOIN analysis.mut_freq mf1 ON mf1.aliquot_barcode = tp.tumor_barcode_a
LEFT JOIN analysis.mut_freq mf2 ON mf2.aliquot_barcode = tp.tumor_barcode_b

-- END --