WITH
paired_classes AS
(
    SELECT
        tp.tumor_pair_barcode,
        tumor_barcode_a,
        tumor_barcode_b,
        class
    FROM analysis.tumor_pairs tp, classes
),
variant_contexts AS
(
    SELECT DISTINCT
        sp.ref_context AS trinucleotide_context,
        sp.alt
    FROM ref.signature_proba sp
    
    /*UNION
    
    SELECT DISTINCT
        reverse_complement(sp.ref_context) AS trinucleotide_context,
        sp.alt
    FROM ref.signature_proba sp*/
),
variant_contexts_pairs AS
(
    SELECT
        tumor_pair_barcode,
        tumor_barcode_a,
        tumor_barcode_b,
        class,
        trinucleotide_context,
        alt
    FROM paired_classes, variant_contexts
),
variant_contexts_counts AS
(
    SELECT
        pg.tumor_pair_barcode,
        pg.tumor_barcode_a,
        pg.tumor_barcode_b,
        (CASE
         WHEN vc.variant_classification_impact = 'LOW' THEN 'syn'
         WHEN vc.variant_classification_impact = 'MODIFIER' THEN 'nc'
         ELSE 'non'
         END) AS class,
        pa.trinucleotide_context,
        pa.alt,
        count(*) AS mut_n
    FROM variants.pgeno pg
    INNER JOIN variants.passanno pa ON pa.variant_id = pg.variant_id
    INNER JOIN variants.variant_classifications vc ON vc.variant_classification = pg.variant_classification
    WHERE
        pg.comparison_type = 'longitudinal' AND 
        pg.variant_type = 'SNP' AND 
        (
            (pg.mutect2_call_a AND NOT pg.mutect2_call_b AND (pg.ref_count_a + pg.alt_count_a) > 14) OR
            (pg.mutect2_call_b AND NOT pg.mutect2_call_a AND (pg.ref_count_b + pg.alt_count_b) > 14) OR 
            (pg.mutect2_call_a AND pg.mutect2_call_b AND (pg.ref_count_a + pg.alt_count_a) > 14 AND (pg.ref_count_b + pg.alt_count_b) > 14)
        )
    GROUP BY 1,2,3,4,5,6
),
variant_contexts_counts_pairs AS
(
    SELECT
        vcp.tumor_pair_barcode,
        vcp.class,
        vcp.trinucleotide_context,
        vcp.alt,
        COALESCE(vcc1.mut_n,0) + COALESCE(vcc2.mut_n, 0) AS mut_n,
        sum(COALESCE(vcc1.mut_n,0) + COALESCE(vcc2.mut_n, 0)) OVER (PARTITION BY vcp.tumor_pair_barcode, vcp.class) AS mut_n_total
    FROM variant_contexts_pairs vcp
    LEFT JOIN variant_contexts_counts vcc1 ON vcc1.tumor_barcode_a = vcp.tumor_barcode_a AND vcc1.tumor_barcode_b = vcp.tumor_barcode_b AND vcc1.class = vcp.class AND vcc1.trinucleotide_context = vcp.trinucleotide_context AND vcc1.alt = vcp.alt
    LEFT JOIN variant_contexts_counts vcc2 ON vcc2.tumor_barcode_a = vcp.tumor_barcode_a AND vcc2.tumor_barcode_b = vcp.tumor_barcode_b AND vcc2.class = vcp.class AND vcc2.trinucleotide_context = reverse_complement(vcp.trinucleotide_context) AND vcc2.alt = reverse_complement(vcp.alt)
),
ref_context_array AS
(
    SELECT array_agg(t.a ORDER BY t.signature) AS ref_context_arr
    FROM (SELECT sp.signature, array_agg(sp.proba ORDER BY sp.ref_context, sp.alt) AS a FROM ref.signature_proba sp GROUP BY sp.signature) t
),
context_reconstruction AS
(
    SELECT
        tumor_pair_barcode,
        class,
        ref_context_arr,
        sum(mut_n) AS mut_n,
        array_agg(mut_n ORDER BY trinucleotide_context, alt) AS array_agg,
        lsqnonneg(ref_context_arr, array_agg(mut_n ORDER BY trinucleotide_context, alt)::double precision[]) AS mut_sigs
    FROM variant_contexts_counts_pairs, ref_context_array
    WHERE mut_n_total > 1
    GROUP BY 1,2,3
)
SELECT
    tumor_pair_barcode,
    class,
    generate_series(1, 30) AS signature,
    mut_n,
    unnest(mut_sigs) AS abs_score,
    unnest(mut_sigs) / (( SELECT sum(s.s) AS sum FROM unnest(mut_sigs) s(s))) AS rel_score
FROM context_reconstruction
ORDER BY 1,2,3

--SELECT * FROM variant_contexts_counts_pairs --DISTINCT trinucleotide_context, alt FROM variant_contexts_pairs

-- END --