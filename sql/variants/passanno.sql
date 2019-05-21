/*
    - - - - - - - - - -
    variants.passanno
    - - - - - - - - - -
    Limited variant annotations (only variants that PASS filters)
    Define a list of variants for which we want to preserve annotations
    In this case meaning all PASS variants and IDH/TERT variants
    We have to specifically retain IDH/TERT due to GATK 4.1.0.0 bug with force-calling
*/
WITH t1 AS (
    SELECT DISTINCT info.variant_id
    FROM variants.info
    WHERE info.filter = 'PASS' OR 
        (info.chrom = 2 AND lower(info.pos) IN (209113112, 209113113)) OR 
        (info.chrom = 5 AND lower(info.pos) IN (1295169, 1295228, 1295242, 1295250)) OR
        (info.chrom = 15 AND lower(info.pos) IN (90631837, 90631838, 90631839))
)
SELECT
    anno.variant_id,
    chrom,
    pos,
    ref,
    alt,
    gene_symbol,
    variant_classification,
    secondary_variant_classification,
    variant_type,
    genome_change,
    transcript,
    transcript_strand,
    transcript_exon,
    transcript_position,
    cdna_change,
    cds_change,
    protein_change,
    gc_content,
    reference_context,
    "substring"(reference_context::text, 10, 3) AS trinucleotide_context
FROM variants.anno
INNER JOIN t1 ON t1.variant_id = anno.variant_id