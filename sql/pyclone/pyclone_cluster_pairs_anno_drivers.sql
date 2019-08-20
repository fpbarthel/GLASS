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
),
selected_genes AS
(
    SELECT DISTINCT sn.gene_symbol, variant_id, chrom, pos, alt, sn.variant_classification, variant_classification_priority, protein_change
    FROM variants.passanno sn
    INNER JOIN ref.driver_genes ds ON ds.gene_symbol = sn.gene_symbol
    LEFT JOIN variants.variant_classifications vc ON sn.variant_classification = vc.variant_classification
    WHERE
        has_mut IS TRUE AND
        ((sn.gene_symbol NOT IN ('TERT','IDH1') AND variant_classification_priority IS NOT NULL) OR
        (sn.gene_symbol = 'TERT' AND sn.variant_classification = '5''Flank' AND lower(sn.pos) IN (1295228,1295250)) OR
        (sn.gene_symbol = 'IDH1' AND sn.protein_change IN ('p.R132C','p.R132G','p.R132H','p.R132S')))
),
selected_genes_geno AS
(
    SELECT DISTINCT case_barcode, cluster_id, string_agg(DISTINCT gene_symbol, ', ') AS drivers
    FROM selected_genes sg
    INNER JOIN variants.passgeno pg ON pg.variant_id = sg.variant_id
    INNER JOIN variants.pyclone_loci pl ON pl.variant_id = sg.variant_id AND pl.aliquot_barcode = pg.aliquot_barcode
    WHERE ssm2_pass_call IS TRUE
    GROUP BY 1,2
)
--SELECT * FROM selected_genes_geno ORDER BY 3 DESC
SELECT pc.case_barcode, idh_codel_subtype, pc.cluster_id, size, ccf_a, ccf_b, rank_a, rank_b, drivers
FROM pyclone_clusters pc
LEFT JOIN selected_genes_geno sgg ON sgg.case_barcode = pc.case_barcode AND sgg.cluster_id = pc.cluster_id