/*
Determine arm level driver changes between selected primary and recurrences
*/
WITH
selected_tumor_pairs AS
(
	SELECT
		tumor_pair_barcode,
		ps.case_barcode,
		idh_codel_subtype,
		tumor_barcode_a,
		tumor_barcode_b
	FROM analysis.tumor_pairs ps
	LEFT JOIN analysis.blocklist b1 ON b1.aliquot_barcode = ps.tumor_barcode_a
	LEFT JOIN analysis.blocklist b2 ON b2.aliquot_barcode = ps.tumor_barcode_b
	LEFT JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
	WHERE
		sample_type_b <> 'M1' AND
		b1.coverage_exclusion = 'allow' AND b2.coverage_exclusion = 'allow' 
),
selected_arms AS
(
	SELECT
		arm,
		chrom,
		pos
	FROM ref.chr_arms ca
	WHERE --ca.chrom NOT IN ('X','Y') AND ca.arm <> '21p'
		ca.chrom IN ('7','10','19','20') OR ca.arm IN ('1p','19q','22q','13q') -- define drivers for arm level
),
selected_arms_pairs AS
(
	SELECT
		tumor_pair_barcode,
		tumor_barcode_a,
		tumor_barcode_b,
		case_barcode,
		idh_codel_subtype,
		arm,
		chrom,
		pos
	FROM selected_tumor_pairs, selected_arms
),
cnv_by_pair_arm AS
(
	SELECT
		sgs.tumor_pair_barcode,
		sgs.case_barcode,
		sgs.idh_codel_subtype,
		sgs.tumor_barcode_a,
		sgs.tumor_barcode_b,
		sgs.chrom,
		sgs.arm,
		first_value(c1.arm_call) OVER w AS arm_a,
		first_value(c2.arm_call) OVER w AS arm_b,
		first_value(CASE
		 WHEN c1.arm_call = c2.arm_call THEN 'S'
		 WHEN c1.arm_call <> 0 AND c2.arm_call = 0 THEN 'P'
		 WHEN c1.arm_call = 0 AND c2.arm_call <> 0 THEN 'R'
		 WHEN c1.arm_call <> 0 AND c2.arm_call <> 0 THEN 'PR'
		 ELSE NULL
		 END) OVER w AS cnv_change,
		first_value(CASE
		 WHEN (c1.arm_call = c2.arm_call OR c1.arm_call = 0 OR c2.arm_call = 0) AND (c1.arm_call = -1 OR c2.arm_call = -1) THEN 'del'
		 WHEN (c1.arm_call = c2.arm_call OR c1.arm_call = 0 OR c2.arm_call = 0) AND (c1.arm_call = 1 OR c2.arm_call = 1) THEN 'amp'
		 WHEN c1.arm_call <> 0 AND c2.arm_call <> 0 THEN 'del/amp'
		 ELSE NULL
		 END) OVER w AS effect,
		row_number() OVER w AS rn,
		COUNT(*) OVER w AS numr
	FROM selected_arms_pairs sgs
	LEFT JOIN analysis.cnv_by_arm_gatk c1 ON c1.aliquot_barcode = sgs.tumor_barcode_a AND c1.arm = sgs.arm
	LEFT JOIN analysis.cnv_by_arm_gatk c2 ON c2.aliquot_barcode = sgs.tumor_barcode_b AND c2.arm = sgs.arm
	--LEFT JOIN seg_stats_optimized ss1 ON ss1.aliquot_barcode = sgs.tumor_barcode_a
	--LEFT JOIN seg_stats_optimized ss2 ON ss2.aliquot_barcode = sgs.tumor_barcode_b
	WHERE
		(abs(c1.arm_call) = 1 OR abs(c2.arm_call) = 1) AND
		c1.arm_call IS NOT NULL AND
		c2.arm_call IS NOT NULL
	WINDOW w AS (PARTITION BY tumor_pair_barcode, sgs.chrom, c1.arm_call, c2.arm_call)
),
cnv_by_pair AS
(
	SELECT
		cpa.case_barcode,
		cpa.tumor_pair_barcode,
		cpa.tumor_barcode_a,
		cpa.tumor_barcode_b,
		cpa.idh_codel_subtype,
		bool_and(cnv_change = 'S') AS shared,
		bool_or(cnv_change = 'P') AS private_a,
		bool_or(cnv_change = 'R') AS private_b,
		bool_or(cnv_change = 'PR') AS private_ab,
		COUNT(DISTINCT cpa.arm) AS driver_count,
		COUNT(CASE WHEN cnv_change = 'S' THEN cpa.arm END) AS driver_count_shared,
		COUNT(CASE WHEN cnv_change = 'P' THEN cpa.arm END) AS driver_count_private_a,
		COUNT(CASE WHEN cnv_change = 'R' THEN cpa.arm END) AS driver_count_private_b,
		COUNT(CASE WHEN cnv_change = 'PR' THEN cpa.arm END) AS driver_count_private_ab,
		trim(BOTH ', ' FROM string_agg(CASE WHEN cnv_change = 'S' THEN (CASE WHEN numr > 1 THEN cpa.chrom ELSE cpa.arm END) || ' ' || effect || ', ' ELSE '' END, '')) AS driver_shared,
		trim(BOTH ', ' FROM string_agg(CASE WHEN cnv_change = 'R' OR cnv_change = 'PR' THEN '+' || (CASE WHEN numr > 1 THEN cpa.chrom ELSE cpa.arm END) || ' ' || effect || ', ' ELSE '' END, '')) AS target_a,
		trim(BOTH ', ' FROM string_agg(CASE WHEN cnv_change = 'P' OR cnv_change = 'PR' THEN '-' || (CASE WHEN numr > 1 THEN cpa.chrom ELSE cpa.arm END) || ' ' || effect || ', ' ELSE '' END, '')) AS target_b
	FROM cnv_by_pair_arm cpa
	WHERE rn = 1
	GROUP BY 1,2,3,4,5
),
cnv_driver_status AS
(
	SELECT
		case_barcode,
		tumor_pair_barcode,
		tumor_barcode_a,
		tumor_barcode_b,
		driver_count,
		driver_count_shared,
		driver_count_private_a,
		driver_count_private_b,
		driver_shared,
		(CASE WHEN shared THEN 'Driver stable'
		 WHEN NOT shared AND private_a AND NOT private_b THEN 'Driver loss'
		 WHEN NOT shared AND private_b AND NOT private_a THEN 'Driver gain'
		 WHEN NOT shared AND private_b AND private_b THEN 'Driver switch' END) driver_status,
		--TRIM(BOTH ', ' FROM target_a || ', ' || target_b) AS target
		target_a,
		target_b
	FROM cnv_by_pair
),
cnv_driver_stability AS
(
	SELECT
		stp.tumor_pair_barcode,
		stp.case_barcode,
		stp.tumor_barcode_a,
		stp.tumor_barcode_b,
		idh_codel_subtype,
		(CASE WHEN driver_count > 0 THEN driver_count ELSE 0 END) AS arm_driver_count,
		driver_count_shared AS arm_driver_count_shared,
		driver_count_private_a AS arm_driver_count_private_a,
		driver_count_private_b AS arm_driver_count_private_b,
		(CASE WHEN driver_shared <> '' THEN driver_shared ELSE NULL END) arm_driver_shared,
		(CASE WHEN driver_status IS NOT NULL THEN driver_status ELSE 'Driver null' END) arm_driver_status,
		(CASE
		 WHEN driver_status IN ('Driver switch','Driver loss','Driver gain') THEN 'Driver unstable'
		 WHEN driver_status IN ('Driver stable') OR driver_status IS NULL THEN 'Driver stable' END) arm_driver_stability,
		(CASE WHEN target_a <> '' THEN target_a ELSE NULL END) arm_driver_change_a,
		(CASE WHEN target_b <> '' THEN target_b ELSE NULL END) arm_driver_change_b
	FROM cnv_driver_status ds
	RIGHT JOIN selected_tumor_pairs stp ON stp.tumor_pair_barcode = ds.tumor_pair_barcode
)
SELECT * FROM cnv_driver_stability ORDER BY 1