WITH
selected_tumor_pairs AS
(
    SELECT tumor_pair_barcode, tumor_barcode_a, tumor_barcode_b, ss.case_barcode, idh_codel_subtype
    FROM analysis.gold_set ss
    INNER JOIN clinical.subtypes st ON st.case_barcode = ss.case_barcode
),
selected_aliquots AS
(
    SELECT tumor_barcode_a AS aliquot_barcode, idh_codel_subtype, case_barcode, 'P' AS sample_type FROM selected_tumor_pairs
    UNION
    SELECT tumor_barcode_b AS aliquot_barcode, idh_codel_subtype, case_barcode, 'R' AS sample_type FROM selected_tumor_pairs
),
pyclone_clusters AS
(
    SELECT sa.case_barcode, sa.idh_codel_subtype, pc.cluster_id, COUNT(*) AS num_samples, min(size) as cluster_size, MIN(mean) as min_ccf, MAX(mean) AS max_ccf, sum(mean)/COUNT(mean) AS mean_ccf
    FROM variants.pyclone_cluster pc
    RIGHT JOIN selected_aliquots sa ON sa.aliquot_barcode = pc.aliquot_barcode
    GROUP BY 1,2,3
    HAVING MIN(mean) > 0.1 OR MAX(mean) > 0.1 OR bool_and(mean IS NULL)
    ORDER BY 5 DESC
)
--SELECT DISTINCT case_barcode FROM pyclone_clusters
SELECT * FROM pyclone_clusters ORDER BY 1,2