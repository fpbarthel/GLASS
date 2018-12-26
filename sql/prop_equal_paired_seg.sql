/*
---
Table: `titan_seg_paired_comparison`
---
For each `tumor_pair_barcode` with paired segmentation data in the `titan_seg_paired_delta` table:
- Quantify the sum of all segment sizes (this should roughly add up to the total genome size)
- Quantify the sum of all segments that did not change in call between (a) and (b) as `delta_eq_size`
- Quantify the proportion of segments that did not change between (a) and (b)
	* a proportion of 1 means that the copy number profile of (a) and (b) are equal
	* a proportion of 0 means that there is no overlap in copy number profile between (a) and (b)

In addition: for each `tumor_pair_barcode` in the `tumor_pairs` table:
- Compare number of segments between (a) and (b)
- Compare proportion of genome that is heterozyous between (a) and (b)
	* a difference of -1 means that (a) is largely heterozygous, whereas (b) is entirely non-heterozygous
	* a difference of 0 indicates (a) and (b) have equal sized heterozygous region
	* a difference of +1 indicates that (b) is heterozygous, whereas (a) is not

Now updated to subset using `selected_tumor_pairs` to filter based on a platinum set including optimal CNV samples only
*/
WITH
selected_tumor_pairs AS
(
	SELECT
		tumor_pair_barcode,
		row_number() OVER (PARTITION BY case_barcode ORDER BY surgical_interval_mo DESC, portion_a ASC, portion_b ASC, substring(tumor_pair_barcode from 27 for 3) ASC) AS priority
	FROM analysis.tumor_pairs ps
	LEFT JOIN analysis.blocklist b1 ON b1.aliquot_barcode = ps.tumor_barcode_a
	LEFT JOIN analysis.blocklist b2 ON b2.aliquot_barcode = ps.tumor_barcode_b
	WHERE
		comparison_type = 'longitudinal' AND
		sample_type_b <> 'M1' AND 													-- exclude metastatic samples here because this is outside the scope of our study
		b1.cnv_exclusion = 'allow' AND b2.cnv_exclusion = 'allow' AND
		b1.coverage_exclusion = 'allow' AND b2.coverage_exclusion = 'allow'
),
paired_seg_quant AS
(
	SELECT
		tumor_pair_barcode,
		tumor_barcode_a,
		tumor_barcode_b,
		sum(upper(pos_intersect) - lower(pos_intersect)) AS delta_seg_size,
		sum(CASE WHEN delta_minor = 0 AND delta_major = 0 THEN upper(pos_intersect) - lower(pos_intersect) ELSE 0 END) AS delta_eq_size
	FROM analysis.titan_seg_paired_delta
	GROUP BY tumor_pair_barcode, tumor_barcode_a, tumor_barcode_b
)
SELECT
	pq.tumor_pair_barcode,
	tumor_barcode_a,
	tumor_barcode_b,
	(EXISTS ( SELECT stp.tumor_pair_barcode, stp.priority FROM selected_tumor_pairs stp WHERE stp.tumor_pair_barcode = pq.tumor_pair_barcode AND stp.priority = 1)) AS diamond_set,
	s1.num_seg AS num_seg_a,
	s2.num_seg AS num_seg_b,
	s1.prop_het AS prop_het_a,
	s2.prop_het AS prop_het_b,
	s2.num_seg - s1.num_seg AS delta_num_seg,
	s2.prop_het - s1.prop_het AS delta_prop_het,
	delta_eq_size,
	delta_seg_size,
	round(delta_eq_size::decimal/delta_seg_size,4) AS prop_delta_eq
FROM paired_seg_quant pq
-- We no longer want to INNER JOIN here because we are using EXISTS in the SELECT clause to determine whether said row is a part of the diamond set
-- INNER JOIN selected_tumor_pairs stp ON stp.tumor_pair_barcode = pq.tumor_pair_barcode AND stp.priority = 1
LEFT JOIN analysis.pairs p1 ON p1.tumor_barcode = pq.tumor_barcode_a
LEFT JOIN analysis.pairs p2 ON p2.tumor_barcode = pq.tumor_barcode_b
LEFT JOIN analysis.titan_seg_prop_het s1 ON s1.pair_barcode = p1.pair_barcode
LEFT JOIN analysis.titan_seg_prop_het s2 ON s2.pair_barcode = p2.pair_barcode
ORDER BY 9 -- or: 12