WITH
selected_tumor_pairs AS
(
	SELECT ss.tumor_pair_barcode, ss.tumor_barcode_a, ss.tumor_barcode_b, ss.case_barcode, idh_codel_subtype, (CASE WHEN gs.tumor_pair_barcode IS NULL THEN 'Silver set' ELSE 'Gold set' END) AS gold_set
	FROM analysis.gold_set ss
	LEFT JOIN analysis.gold_set gs ON gs.tumor_pair_barcode = ss.tumor_pair_barcode
	INNER JOIN clinical.subtypes st ON st.case_barcode = ss.case_barcode
),
selected_arms AS
(
	SELECT chrom,arm,direction FROM ref.arm_drivers_subtype
),
cnv_by_pair_arm AS
(
	SELECT
		stp.tumor_pair_barcode,
		stp.case_barcode,
		stp.idh_codel_subtype,
		stp.tumor_barcode_a,
		stp.tumor_barcode_b,
		sa.chrom,
		sa.arm,
		c1.arm_call AS arm_a,
		c2.arm_call AS arm_b,
		(CASE
		 WHEN sa.direction = -1 AND (c1.arm_call = -1 OR c2.arm_call = -1) THEN 'del'
		 WHEN sa.direction = 1 AND (c1.arm_call = 1 OR c2.arm_call = 1) THEN 'amp'
		 WHEN sa.direction = -1 AND (c1.arm_call = 1 OR c2.arm_call = 1) THEN 'neut'
		 WHEN sa.direction = 1 AND (c1.arm_call = -1 OR c2.arm_call = -1) THEN 'neut'
		 WHEN (c1.arm_call = 0 OR c2.arm_call = 0) THEN 'neut'
		 ELSE NULL
		 END) cnv_state,
		(CASE
		 WHEN direction = -1 AND (c1.arm_call = -1 OR c2.arm_call = -1) AND c1.arm_call < c2.arm_call THEN 'P'
		 WHEN direction = -1 AND (c1.arm_call = -1 OR c2.arm_call = -1) AND (c1.arm_call = c2.arm_call OR (c1.arm_call IS NULL OR c2.arm_call IS NULL)) THEN 'S'
		 WHEN direction = -1 AND (c1.arm_call = -1 OR c2.arm_call = -1) AND c1.arm_call > c2.arm_call THEN 'R'
		 WHEN direction = 1 AND (c1.arm_call = 1 OR c2.arm_call = 1) AND c1.arm_call > c2.arm_call THEN 'P'
		 WHEN direction = 1 AND (c1.arm_call = 1 OR c2.arm_call = 1) AND (c1.arm_call = c2.arm_call OR (c1.arm_call IS NULL OR c2.arm_call IS NULL)) THEN 'S'
		 WHEN direction = 1 AND (c1.arm_call = 1 OR c2.arm_call = 1) AND c1.arm_call < c2.arm_call THEN 'R'
		 WHEN (c1.arm_call = 0 OR c2.arm_call = 0) AND (c1.arm_call = c2.arm_call OR (c1.arm_call IS NULL OR c2.arm_call IS NULL)) THEN 'S'
		 WHEN direction = -1 AND (c1.arm_call = 1 OR c2.arm_call = 1) AND (c1.arm_call = c2.arm_call OR (c1.arm_call IS NULL OR c2.arm_call IS NULL)) THEN 'S'
		 WHEN direction = 1 AND (c1.arm_call = -1 OR c2.arm_call = -1) AND (c1.arm_call = c2.arm_call OR (c1.arm_call IS NULL OR c2.arm_call IS NULL)) THEN 'S'
		 ELSE NULL
		 END) cnv_change
	FROM selected_tumor_pairs stp
	CROSS JOIN selected_arms sa
	LEFT JOIN analysis.gatk_cnv_by_arm c1 ON c1.aliquot_barcode = stp.tumor_barcode_a AND c1.arm = sa.arm
	LEFT JOIN analysis.gatk_cnv_by_arm c2 ON c2.aliquot_barcode = stp.tumor_barcode_b AND c2.arm = sa.arm
)
SELECT * FROM cnv_by_pair_arm ORDER BY 1