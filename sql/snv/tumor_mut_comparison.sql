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
		FROM variants.passgeno gt
		WHERE gt.aliquot_barcode = tp.tumor_barcode_a AND gt.ad_ref + gt.ad_alt > 14 AND ssm2_pass_call IS TRUE) AS count_a,
	
	( 	SELECT count(*) AS count
		FROM variants.passgeno gt
		WHERE gt.aliquot_barcode = tp.tumor_barcode_b AND gt.ad_ref + gt.ad_alt > 14 AND ssm2_pass_call IS TRUE) AS count_b,
		
	(	SELECT count(*) AS count
		FROM (	SELECT variant_id
				FROM variants.passgeno gt
				WHERE gt.aliquot_barcode = tp.tumor_barcode_a AND gt.ad_ref + gt.ad_alt > 14 AND ssm2_pass_call IS TRUE
				UNION
				SELECT variant_id
				FROM variants.passgeno gt
				WHERE gt.aliquot_barcode = tp.tumor_barcode_b AND gt.ad_ref + gt.ad_alt > 14 AND ssm2_pass_call IS TRUE) t) AS union_ab,
		
	(	SELECT count(*) AS count
		FROM (	SELECT variant_id
				FROM variants.passgeno gt
				WHERE gt.aliquot_barcode = tp.tumor_barcode_a AND gt.ad_ref + gt.ad_alt > 14 AND ssm2_pass_call IS TRUE
				INTERSECT
				SELECT variant_id
				FROM variants.passgeno gt
				WHERE gt.aliquot_barcode = tp.tumor_barcode_b AND gt.ad_ref + gt.ad_alt > 14 AND ssm2_pass_call IS TRUE) t) AS intersection_ab,
		
	(	SELECT count(*) AS count
		FROM (	SELECT variant_id
				FROM variants.passgeno gt
				WHERE gt.aliquot_barcode = tp.tumor_barcode_a AND gt.ad_ref + gt.ad_alt > 14 AND ssm2_pass_call IS TRUE
				EXCEPT
				SELECT variant_id
				FROM variants.passgeno gt
				WHERE gt.aliquot_barcode = tp.tumor_barcode_b AND gt.ad_ref + gt.ad_alt > 14 AND ssm2_pass_call IS TRUE) t) AS setdiff_a,
	
	(	SELECT count(*) AS count
		FROM ( SELECT variant_id
				FROM variants.passgeno gt
				WHERE gt.aliquot_barcode = tp.tumor_barcode_b AND gt.ad_ref + gt.ad_alt > 14 AND ssm2_pass_call IS TRUE
				EXCEPT
				SELECT variant_id
				FROM variants.passgeno gt
				WHERE gt.aliquot_barcode = tp.tumor_barcode_a AND gt.ad_ref + gt.ad_alt > 14 AND ssm2_pass_call IS TRUE) t) AS setdiff_b
	 
FROM analysis.tumor_pairs tp