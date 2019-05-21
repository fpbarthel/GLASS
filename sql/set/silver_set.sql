WITH
selected_tumor_pairs AS
(
    SELECT
        ps.tumor_pair_barcode,
        ps.case_barcode,
        ps.tumor_barcode_a,
        ps.tumor_barcode_b,
        row_number() OVER (PARTITION BY ps.case_barcode ORDER BY ps.surgical_interval_mo DESC, ps.portion_a, ps.portion_b, ("substring"(ps.tumor_pair_barcode, 27, 3))) AS priority
    FROM analysis.tumor_pairs ps
    LEFT JOIN analysis.blocklist b1 ON b1.aliquot_barcode = ps.tumor_barcode_a
    LEFT JOIN analysis.blocklist b2 ON b2.aliquot_barcode = ps.tumor_barcode_b
    WHERE ps.comparison_type = 'longitudinal' AND ps.sample_type_b <> 'M1' AND b1.coverage_exclusion = 'allow' AND b2.coverage_exclusion = 'allow'
)
SELECT
    tumor_pair_barcode,
    case_barcode,
    tumor_barcode_a,
    tumor_barcode_b
FROM selected_tumor_pairs
WHERE selected_tumor_pairs.priority = 1