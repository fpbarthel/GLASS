/*
This is the initial definition but was deprecated out of confusion
This definition could include cases (patients) that are also in the silver set, but using a different combionation of primary and recurrence
The new definition (not commented, below) instead takes the subset of the silver set
Note that for GISTIC we instead took the deprecated gold set to define a set of unique primaries and unique recurrences
===
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
    WHERE ps.comparison_type = 'longitudinal' AND ps.sample_type_b <> 'M1' AND b1.coverage_exclusion = 'allow' AND b2.coverage_exclusion = 'allow' AND b1.cnv_exclusion != 'block' AND b2.cnv_exclusion != 'block'
)
SELECT
    tumor_pair_barcode,
    case_barcode,
    tumor_barcode_a,
    tumor_barcode_b
FROM selected_tumor_pairs
WHERE selected_tumor_pairs.priority = 1*/

SELECT *
FROM analysis.silver_set ss
INNER JOIN analysis.blocklist bl1 ON bl1.aliquot_barcode = ss.tumor_barcode_a
INNER JOIN analysis.blocklist bl2 ON bl2.aliquot_barcode = ss.tumor_barcode_b
WHERE bl1.cnv_exclusion != 'block' AND bl2.cnv_exclusion != 'block'