WITH
selected_tumor_pairs AS
(
    SELECT tumor_pair_barcode, tumor_barcode_a, tumor_barcode_b, ss.case_barcode, idh_codel_subtype
    FROM analysis.gold_set ss
    INNER JOIN clinical.subtypes st ON st.case_barcode = ss.case_barcode
),
pyclone_clusters AS
(
    SELECT stp.case_barcode, stp.idh_codel_subtype, pc1.cluster_id, pc1.size AS size, pc1.mean AS ccf_a, pc2.mean AS ccf_b,
    	(RANK() OVER (PARTITION BY stp.case_barcode ORDER BY pc1.mean DESC))::integer AS rank_a,
    	(RANK() OVER (PARTITION BY stp.case_barcode ORDER BY pc2.mean DESC))::integer AS rank_b
    FROM selected_tumor_pairs stp
    INNER JOIN variants.pyclone_cluster pc1 ON pc1.aliquot_barcode = stp.tumor_barcode_a
    INNER JOIN variants.pyclone_cluster pc2 ON pc2.aliquot_barcode = stp.tumor_barcode_b AND pc2.cluster_id = pc1.cluster_id
    WHERE pc1.size > 1 AND (pc1.mean > 0.1 OR pc2.mean > 0.1)
)
SELECT * FROM pyclone_clusters