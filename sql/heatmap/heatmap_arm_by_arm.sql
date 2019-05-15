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
		 WHEN direction = -1 AND c1.arm_call = -1 THEN true
		 WHEN direction = 1 AND c1.arm_call = 1 THEN true
		 ELSE false
		 END) selected_call_a,
		(CASE
		 WHEN direction = -1 AND c2.arm_call = -1 THEN true
		 WHEN direction = 1 AND c2.arm_call = 1 THEN true
		 ELSE false
		 END) selected_call_b,
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
),
arm_by_arm AS
(
	SELECT
		arm,
		SUM(CASE WHEN selected_call_a THEN 1 ELSE 0 END)::integer		AS count_a,											-- count of mutations in A
		SUM(CASE WHEN selected_call_b THEN 1 ELSE 0 END)::integer		AS count_b,											-- count of mutations in B
		SUM(CASE WHEN selected_call_a AND selected_call_b THEN 1 ELSE 0 END)::integer 		AS shared, 						-- if both A and B are true then a variant is shared
		SUM(CASE WHEN selected_call_a AND NOT selected_call_b THEN 1 ELSE 0 END)::integer 	AS private_a,					-- if A is true and B is not then it is unique to A
		SUM(CASE WHEN selected_call_b AND NOT selected_call_a THEN 1 ELSE 0 END)::integer 	AS private_b,					-- if B is true and A is not then it is unique to B,
		SUM(CASE WHEN selected_call_b OR selected_call_a THEN 1 ELSE 0 END)::integer 		AS total,						-- if a is true OR B is true, sums to total variants between A and B
		ROUND(SUM(CASE WHEN selected_call_a AND selected_call_b THEN 1 ELSE 0 END) / 
			  SUM(CASE WHEN selected_call_b OR selected_call_a THEN 1 ELSE 0 END)::decimal,2) AS prop_shared,		-- divides `shared` by total to get the proportion shared
		ROUND(SUM(CASE WHEN selected_call_a AND NOT selected_call_b THEN 1 ELSE 0 END) /
			  SUM(CASE WHEN selected_call_b OR selected_call_a THEN 1 ELSE 0 END)::decimal,2) AS prop_private_a,	-- divides `private_a` by total to get the proportion private to A
		ROUND(SUM(CASE WHEN selected_call_b AND NOT selected_call_a THEN 1 ELSE 0 END) /
			  SUM(CASE WHEN selected_call_b OR selected_call_a THEN 1 ELSE 0 END)::decimal,2) AS prop_private_b		-- divides `private_b` by total to get the proportion private to B
	FROM cnv_by_pair_arm
	GROUP BY 1	-- aggregate across genes
	ORDER BY 7 DESC -- order by column 7 (= total)
),
arm_by_arm_subtype AS
(
	SELECT
		idh_codel_subtype,
		arm,	
		SUM(CASE WHEN selected_call_a THEN 1 ELSE 0 END)::integer		AS count_a,											-- count of mutations in A
		SUM(CASE WHEN selected_call_b THEN 1 ELSE 0 END)::integer		AS count_b,											-- count of mutations in B
		SUM(CASE WHEN selected_call_a AND selected_call_b THEN 1 ELSE 0 END)::integer 		AS shared, 						-- if both A and B are true then a variant is shared
		SUM(CASE WHEN selected_call_a AND NOT selected_call_b THEN 1 ELSE 0 END)::integer 	AS private_a,					-- if A is true and B is not then it is unique to A
		SUM(CASE WHEN selected_call_b AND NOT selected_call_a THEN 1 ELSE 0 END)::integer 	AS private_b,					-- if B is true and A is not then it is unique to B,
		SUM(CASE WHEN selected_call_b OR selected_call_a THEN 1 ELSE 0 END)::integer 		AS total,						-- if a is true OR B is true, sums to total variants between A and B
		ROUND(SUM(CASE WHEN selected_call_a AND selected_call_b THEN 1 ELSE 0 END) / 
			  NULLIF(SUM(CASE WHEN selected_call_b OR selected_call_a THEN 1 ELSE 0 END),0)::decimal,2) AS prop_shared,		-- divides `shared` by total to get the proportion shared
		ROUND(SUM(CASE WHEN selected_call_a AND NOT selected_call_b THEN 1 ELSE 0 END) /
			  NULLIF(SUM(CASE WHEN selected_call_b OR selected_call_a THEN 1 ELSE 0 END),0)::decimal,2) AS prop_private_a,	-- divides `private_a` by total to get the proportion private to A
		ROUND(SUM(CASE WHEN selected_call_b AND NOT selected_call_a THEN 1 ELSE 0 END) /
			  NULLIF(SUM(CASE WHEN selected_call_b OR selected_call_a THEN 1 ELSE 0 END),0)::decimal,2) AS prop_private_b		-- divides `private_b` by total to get the proportion private to B
	FROM cnv_by_pair_arm vg
	GROUP BY 1,2	-- aggregate across genes
	ORDER BY 7 DESC -- order by column 7 (= total)
)
(SELECT 'all' AS idh_codel_subtype, vg.*
FROM arm_by_arm vg
		 
UNION
		 
SELECT vs.*
FROM arm_by_arm_subtype vs)
ORDER BY 1 ASC, 8 DESC