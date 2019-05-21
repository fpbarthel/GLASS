SELECT gs.case_barcode, case_age_diagnosis_years AS age, surgical_interval, idh_codel_subtype, hypermutator_status::integer, fraction, signature, mut_n, abs_score, rel_score, RANK() OVER (PARTITION BY gs.case_barcode,fraction ORDER BY rel_score DESC) AS rnk, COUNT(*) OVER (PARTITION BY gs.case_barcode, signature) AS all_fractions_counts
FROM analysis.mut_sig_fraction_limited ms
INNER JOIN analysis.gold_set gs ON gs.tumor_pair_barcode = ms.tumor_pair_barcode
INNER JOIN analysis.tumor_clinical_comparison tcc ON tcc.tumor_pair_barcode = ms.tumor_pair_barcode
INNER JOIN clinical.subtypes st ON st.case_barcode = gs.case_barcode
INNER JOIN clinical.cases ca ON ca.case_barcode = gs.case_barcode