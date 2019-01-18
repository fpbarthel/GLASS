/*
Taylor method to determine aneuploidy
- GATK4 CallSegments is used to assign + (amplification), - (deletion) or 0 (neutral) call to each segment
- Allosomal segments and segments from chr 21p are removed to match the published method
- Adjacent segments on the same chromosome arm with identical calls are merged/reduced into a single segment with a single call
- This results in one segment per chromosome arm at the minimum, 
	eg. if the input contains one segment for the entire chromosome it will be split into two segments for each arm, 
	but if the input contains >=3 segments one of the arms will be split into multiple segments
- For each **reduced/merged** segment, we compute the proportion of the chromosome arm it spans
- Segments greater than 80% of the arm length result in a call of either -1 (loss), 0 (neutral) or +1 (amp) to the entire arm
- For each aliquot count the number of arms with a +1 or -1 event as aneuploidy
- The resulting aneuploidy_score ranges from 0-39
*/
WITH
arrange_segments_by_arm AS
(
    SELECT
		aliquot_barcode,
		gs.chrom,
		ca.arm,
		ca.pos * gs.pos AS pos,
		2^log2_copy_ratio::decimal AS cr,
		gs.cnv_call AS cnv,
		(upper(ca.pos * gs.pos)::decimal - lower(ca.pos * gs.pos) -1) as seg_size,
		(sum((upper(ca.pos * gs.pos)::decimal - lower(ca.pos * gs.pos) -1)) OVER w) as arm_size,
		row_number() OVER w2 - row_number() OVER w3 AS grp
    FROM analysis.gatk_seg_call gs
	INNER JOIN ref.chr_arms ca ON ca.chrom = gs.chrom AND ca.pos && gs.pos
	WHERE gs.chrom NOT IN ('X','Y') AND ca.arm <> '21p' AND aliquot_barcode LIKE 'GLSS-CU-%'
	WINDOW w AS (PARTITION BY gs.aliquot_barcode, ca.arm), w2 AS (PARTITION BY gs.aliquot_barcode, ca.arm ORDER BY ca.pos * gs.pos), w3 AS (PARTITION BY gs.aliquot_barcode, ca.arm, gs.cnv_call ORDER BY ca.pos * gs.pos)
),
join_adjacent_segments AS
(
	SELECT
		aliquot_barcode,
		chrom,
		arm,
		grp,
		int4range(min(lower(pos)),max(upper(pos))) AS pos,
		cnv,
		COUNT(*) AS num_seg,
		sum(seg_size * cr) / sum(seg_size) AS wcr,
		(CASE WHEN COUNT(*) > 1 THEN sqrt( (sum(seg_size * cr^2) - (sum(seg_size * cr)^2)/sum(seg_size)) / (sum(seg_size)-1) ) ELSE 0 END) AS wsd,
		sum(seg_size) AS seg_size,
		min(arm_size) AS arm_size
	FROM arrange_segments_by_arm t1
	WHERE cr > 0 AND seg_size > 0
	GROUP BY 1,2,3,4,6
	ORDER BY 1,2,3,5
),
sort_by_longest_segment AS
(
	SELECT
		aliquot_barcode,
		chrom,
		arm,
		grp,
		cnv,
		row_number() OVER w AS partion_id,
		seg_size/arm_size AS prop_arm,
		num_seg,
		wcr,
		wsd
	FROM join_adjacent_segments t2
	WINDOW w AS (PARTITION BY aliquot_barcode, chrom, arm ORDER BY seg_size/arm_size DESC)
),
call_arm AS
(
	SELECT
		aliquot_barcode,
		chrom,
		arm,
		(CASE
		 WHEN prop_arm > 0.80 AND cnv = 1 THEN 1
		 WHEN prop_arm > 0.80 AND cnv = -1 THEN -1
		 WHEN prop_arm > 0.80 AND cnv = 0 THEN 0
		 ELSE NULL END) AS arm_call,
		num_seg AS arm_num_seg,
		wcr AS arm_cr_wmean,
		wsd AS arm_cr_wsd
	FROM sort_by_longest_segment t3
	WHERE partion_id = 1
),
call_aliquot AS
(
	SELECT
		aliquot_barcode,
		sum(CASE WHEN arm_call <> 0 THEN 1 ELSE 0 END) OVER w AS aneuploidy_score,
	arm_cr_wmean,
		ROUND(first_value(arm_cr_wmean) OVER w2,6) AS loss_wmean,
		ROUND(first_value(arm_cr_wsd) OVER w2,6) AS loss_wsd,
		ROUND(first_value(arm_cr_wmean) OVER w3,6) AS gain_wmean,
		ROUND(first_value(arm_cr_wsd) OVER w3,6) AS gain_wsd,
		row_number() OVER w AS rn
	--FROM call_arm t4
	FROM analysis.cnv_by_arm_gatk_v2
	WINDOW
		w  AS (PARTITION BY aliquot_barcode),
		w2 AS (PARTITION BY aliquot_barcode ORDER BY arm_cr_wmean ASC),
		w3 AS (PARTITION BY aliquot_barcode ORDER BY arm_cr_wmean DESC)
)
SELECT * FROM call_aliquot WHERE rn = 1--LIMIT 1000
--SELECT t5.aliquot_barcode, t5.aneuploidy_score AS newscore, ta.aneuploidy_score AS oldscore FROM t5
--LEFT JOIN analysis.taylor_aneuploidy ta ON ta.aliquot_barcode = t5.aliquot_barcode ORDER BY 2