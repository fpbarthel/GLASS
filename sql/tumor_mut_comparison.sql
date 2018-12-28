SELECT
	tp.tumor_pair_barcode,
	tp.case_barcode,
	tp.tumor_barcode_a,
	tp.tumor_barcode_b,
	tp.sample_type_a,
	tp.sample_type_b,
	tp.portion_a,
	tp.portion_b,
	tp.comparison_type,
	tp.surgical_interval_mo,
	
	( 	SELECT count(*) AS count
		FROM analysis.called_genotypes gt
		WHERE gt.aliquot_barcode = tp.tumor_barcode_a AND gt.read_depth > 14) AS count_a,
	
	( 	SELECT count(*) AS count
		FROM analysis.called_genotypes gt
		WHERE gt.aliquot_barcode = tp.tumor_barcode_b AND gt.read_depth > 14) AS count_b,
		
	(	SELECT count(*) AS count
		FROM (	SELECT
					gt.chrom,
					gt.start,
					gt."end",
					gt.alt
				FROM analysis.called_genotypes gt
				WHERE gt.aliquot_barcode = tp.tumor_barcode_a AND gt.read_depth > 14
				UNION
				SELECT
					gt.chrom,
					gt.start,
					gt."end",
					gt.alt
				FROM analysis.called_genotypes gt
				WHERE gt.aliquot_barcode = tp.tumor_barcode_b AND gt.read_depth > 14) t) AS union_ab,
		
	(	SELECT count(*) AS count
		FROM (	SELECT
					gt.chrom,
					gt.start,
					gt."end",
					gt.alt
				FROM analysis.called_genotypes gt
				WHERE gt.aliquot_barcode = tp.tumor_barcode_a AND gt.read_depth > 14
				INTERSECT
				SELECT
					gt.chrom,
					gt.start,
					gt."end",
					gt.alt
				FROM analysis.called_genotypes gt
				WHERE gt.aliquot_barcode = tp.tumor_barcode_b AND gt.read_depth > 14) t) AS intersection_ab,
		
	(	SELECT count(*) AS count
		FROM (	SELECT
					gt.chrom,
					gt.start,
					gt."end",
					gt.alt
				FROM analysis.called_genotypes gt
				WHERE gt.aliquot_barcode = tp.tumor_barcode_a AND gt.read_depth > 14
				EXCEPT
				SELECT
					gt.chrom,
					gt.start,
					gt."end",
					gt.alt
				FROM analysis.called_genotypes gt
				WHERE gt.aliquot_barcode = tp.tumor_barcode_b AND gt.read_depth > 14) t) AS setdiff_a,
	
	(	SELECT count(*) AS count
		FROM ( SELECT
					gt.chrom,
					gt.start,
					gt."end",
					gt.alt
				FROM analysis.called_genotypes gt
				WHERE gt.aliquot_barcode = tp.tumor_barcode_b AND gt.read_depth > 14
				EXCEPT
				SELECT
					gt.chrom,
					gt.start,
					gt."end",
					gt.alt
				FROM analysis.called_genotypes gt
				WHERE gt.aliquot_barcode = tp.tumor_barcode_a AND gt.read_depth > 14) t) AS setdiff_b
	 
FROM analysis.tumor_pairs tp
LEFT JOIN analysis.blocklist b1 ON b1.aliquot_barcode = tp.tumor_barcode_a
LEFT JOIN analysis.blocklist b2 ON b2.aliquot_barcode = tp.tumor_barcode_b
WHERE b1.coverage_exclusion = 'allow'::bpchar AND b2.coverage_exclusion = 'allow'::bpchar;
