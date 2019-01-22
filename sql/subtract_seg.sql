/*
For each `tumor_pair_barcode`, consisting of `tumor_barcode_a` (a) and `tumor_barcode_b` (b) in the `tumor_pairs` table:
- Intersect all segments from (a) and (b) --> `pos_intersect`
- Determine change in `minor_cn` between (a) and (b) as (b-a)
- Determine change in `major_cn` between (a) and (b) as (b-a)
	* For the two items above, a value of 0 indicates no change, positive numbers indicate (b > a) and negative numbers (a > b)
- Determine the difference in log2(tumor/normal) copy number between (a) and (b) as (b-a)
*/
SELECT
	pa.tumor_pair_barcode, 
	pa.case_barcode, 
	pa.tumor_barcode_a, 
	pa.tumor_barcode_b, 
	s1.chrom,
	s1.pos * s2.pos AS pos, 
	s2.minor_cn - s1.minor_cn AS delta_minor, 
	s2.major_cn - s1.major_cn AS delta_major,
	t1.ploidy AS ploidy_a,
	t2.ploidy AS ploidy_b,
	t1.purity AS purity_a,
	t2.purity AS purity_b,
	s2.logr_copy_number - s1.logr_copy_number AS delta_cn,
	(s2.logr_copy_number/t2.ploidy) - (s1.logr_copy_number/t1.ploidy) AS delta_cn_ploidy_adj
FROM analysis.tumor_pairs pa
LEFT JOIN analysis.pairs p1 ON p1.tumor_barcode = pa.tumor_barcode_a
LEFT JOIN analysis.pairs p2 ON p2.tumor_barcode = pa.tumor_barcode_b
LEFT JOIN analysis.titan_seg s1 ON s1.pair_barcode = p1.pair_barcode
LEFT JOIN analysis.titan_seg s2 ON s2.pair_barcode = p2.pair_barcode AND s1.chrom = s2.chrom AND s1.pos && s2.pos
LEFT JOIN analysis.titan_params t1 ON t1.pair_barcode = p1.pair_barcode
LEFT JOIN analysis.titan_params t2 ON t2.pair_barcode = p2.pair_barcode