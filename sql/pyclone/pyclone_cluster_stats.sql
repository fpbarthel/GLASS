WITH
selected_tumor_pairs AS
(
    SELECT tumor_pair_barcode, tumor_barcode_a, tumor_barcode_b, ss.case_barcode, idh_codel_subtype
    FROM analysis.gold_set ss
    INNER JOIN clinical.subtypes st ON st.case_barcode = ss.case_barcode
),
selected_aliquots AS
(
    SELECT tumor_barcode_a AS aliquot_barcode, idh_codel_subtype AS subtype, case_barcode, 'P' AS sample_type FROM selected_tumor_pairs
    UNION
    SELECT tumor_barcode_b AS aliquot_barcode, idh_codel_subtype AS subtype, case_barcode, 'R' AS sample_type FROM selected_tumor_pairs
),
driver_counts AS
(
    SELECT case_barcode, cluster_id, COUNT(*) AS num_drivers
    FROM variants.passanno pa
    INNER JOIN ref.driver_genes dg ON dg.gene_symbol = pa.gene_symbol
    INNER JOIN variants.pyclone_loci pl ON pl.variant_id = pa.variant_id
    INNER JOIN biospecimen.aliquots al ON al.aliquot_barcode = pl.aliquot_barcode
    INNER JOIN biospecimen.samples sa ON sa.sample_barcode = al.sample_barcode
    WHERE has_mut IS TRUE AND variant_allele_frequency > 0
    GROUP BY 1,2    
),
pyclone_clusters AS
(
    SELECT sa.case_barcode, pc.cluster_id, COUNT(*) AS num_samples, min(size) as cluster_size, MIN(mean) as min_ccf, MAX(mean) AS max_ccf, sum(mean)/COUNT(mean) AS mean_ccf
    FROM variants.pyclone_cluster pc
    --INNER JOIN biospecimen.aliquots al ON al.aliquot_barcode = pc.aliquot_barcode
    INNER JOIN selected_aliquots sa ON sa.aliquot_barcode = pc.aliquot_barcode
    --INNER JOIN biospecimen.samples sa ON sa.sample_barcode = al.sample_barcode
    GROUP BY 1,2
    ORDER BY 5 DESC
    HAVING MIN(mean) > 0.1 OR MAX(mean) > 0.1
),
t2 AS
(
    SELECT t1.case_barcode,
        COUNT(t1.cluster_id) AS num_clust,
        COUNT(CASE WHEN min_ccf > 0.25 THEN 1 END) AS num_clonal,
        bool_or(min_ccf > 0.25 AND num_drivers > 0)::integer AS any_clonal_drivers,
        bool_or(min_ccf <= 0.25 AND num_drivers > 0)::integer AS any_subclonal_drivers
        --COUNT(CASE WHEN min_ccf > 0.25 THEN num_drivers END) AS num_drivers_clonal,
        --COUNT(CASE WHEN min_ccf <= 0.25 THEN num_drivers END) AS num_drivers_subclonal
    FROM pyclone_clusters t1
    LEFT JOIN driver_counts t0 ON t0.case_barcode = t1.case_barcode AND t0.cluster_id = t1.cluster_id
    GROUP BY 1
    ORDER BY 3 DESC
)
SELECT * FROM t2