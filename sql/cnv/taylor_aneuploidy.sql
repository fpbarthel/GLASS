/*
---
Modified Taylor method to determine aneuploidy
Postgres implementation of the method described in Taylor et al, Cancer Cell 2018
URL: https://www.cell.com/cancer-cell/fulltext/S1535-6108(18)30111-9
Additionally determines weighted mean and sd of continuous segments so these values can be used as a cutoff for gene-level copy number calling
---
1. Postgres port of ReCapSeg (modelled after GATK4 CallCopyRatioSegments) is used to assign a +1 (amplification), -1 (deletion) or 0 (neutral) call to each segment
2. Allosomal segments and segments from chr 21p are removed to match the published method
3. Adjacent segments on the same chromosome arm with identical calls are merged/reduced into a single segment with a single call
   This results in one segment per chromosome arm at the minimum, 
	 eg. if the input contains one segment for the entire chromosome it will be split into two segments for each arm, 
	 but if the input contains >=3 segments one of the arms will be split into multiple segments
4. For each **reduced/merged** segment, we compute:
	- The proportion of the chromosome arm it spans
	- The weighted mean (wmean) and standard deviation (wsd) non-log2 copy ratio
	- The wmean and wsd are calculated again after temporarily removing any segments not within 2WSD 
5. Segments greater than 80% of the arm length result in a call of either -1 (loss), 0 (neutral) or +1 (amp) to the entire arm
	- Additionally each arm is assigned a weighted mean/sd based on its longest segment
6. For each aliquot count the number of arms with a +1 or -1 event as aneuploidy
	- The resulting aneuploidy_score ranges from 0-39
	- For each sample, note the weighted mean and sd of the most deeply deleted and most deeply amplified chromosome arm 
	- This metric can be used to threshold deep amplifications and deep deletions
---
This code is the source for multiple tables:
o The `call_arm` table described here is the source for `cnv_by_arm_gatk`
o The `call_aliquot` table described here summarizes the arm level table and is saved as `taylor_aneuploidy`
---
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
	WHERE gs.chrom NOT IN ('X','Y') AND ca.arm <> '21p'
	WINDOW w AS (PARTITION BY gs.aliquot_barcode, ca.arm), w2 AS (PARTITION BY gs.aliquot_barcode, ca.arm ORDER BY ca.pos * gs.pos), w3 AS (PARTITION BY gs.aliquot_barcode, ca.arm, gs.cnv_call ORDER BY ca.pos * gs.pos)
),
unfiltered_join_adjacent_segments AS
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
filtered_join_adjacent_segments AS
(
	SELECT
		t1.aliquot_barcode,
		t1.chrom,
		t1.arm,
		t1.grp,
		int4range(min(lower(t1.pos)),max(upper(t1.pos))) AS pos,
		t1.cnv,
		COUNT(*) AS fnum_seg,
		sum(t1.seg_size * cr) / sum(t1.seg_size) AS fwcr,
		(CASE WHEN COUNT(*) > 1 THEN sqrt( (sum(t1.seg_size * cr^2) - (sum(t1.seg_size * cr)^2)/sum(t1.seg_size)) / (sum(t1.seg_size)-1) ) ELSE 0 END) AS fwsd,
		sum(t1.seg_size) AS seg_size,
		min(t1.arm_size) AS arm_size
	FROM arrange_segments_by_arm t1
	INNER JOIN unfiltered_join_adjacent_segments us ON us.aliquot_barcode = t1.aliquot_barcode AND us.chrom = t1.chrom AND us.arm = t1.arm AND us.grp = t1.grp AND us.cnv = t1.cnv
	WHERE
		cr > 0 AND
		t1.seg_size > 0 AND
		(cr - wcr) > -2.0 * wsd AND
		(cr - wcr) < 2.0 * wsd
	GROUP BY 1,2,3,4,6
	ORDER BY 1,2,3,5
),
sort_by_longest_segment AS
(
	SELECT
		t2.aliquot_barcode,
		t2.chrom,
		t2.arm,
		t2.grp,
		t2.cnv,
		row_number() OVER w AS partion_id,
		t2.seg_size/t2.arm_size AS prop_arm,
		num_seg,
		wcr,
		wsd,
		fss.fnum_seg,
		fss.fwcr,
		fss.fwsd
	FROM unfiltered_join_adjacent_segments t2
	LEFT JOIN filtered_join_adjacent_segments fss ON fss.aliquot_barcode = t2.aliquot_barcode AND fss.chrom = t2.chrom AND fss.arm = t2.arm AND fss.grp = t2.grp AND fss.cnv = t2.cnv
	WINDOW w AS (PARTITION BY t2.aliquot_barcode, t2.chrom, t2.arm ORDER BY t2.seg_size/t2.arm_size DESC)
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
		fnum_seg AS arm_num_seg,
		fwcr AS arm_cr_wmean,
		fwsd AS arm_cr_wsd
	FROM sort_by_longest_segment t3
	WHERE partion_id = 1
),
call_aliquot AS
(
	SELECT
		aliquot_barcode,
		sum(CASE WHEN arm_call <> 0 THEN 1 ELSE 0 END) OVER w AS aneuploidy_score,
		sum(CASE WHEN arm_call = 1 THEN 1 ELSE 0 END) OVER w AS aneuploidy_amp_score,
		sum(CASE WHEN arm_call = -1 THEN 1 ELSE 0 END) OVER w AS aneuploidy_del_score,
		first_value(arm_num_seg) OVER w2 AS max_loss_arm_n,
		first_value(arm_cr_wmean) OVER w2 AS max_loss_arm_wmean,
		first_value(arm_cr_wsd) OVER w2 AS max_loss_arm_wsd,
		first_value(arm_num_seg) OVER w3 AS max_gain_arm_n,
		first_value(arm_cr_wmean) OVER w3 AS max_gain_arm_wmean,
		first_value(arm_cr_wsd) OVER w3 AS max_gain_arm_wsd,
		row_number() OVER w AS rn
	FROM call_arm
	WINDOW
		w  AS (PARTITION BY aliquot_barcode),
		w2 AS (PARTITION BY aliquot_barcode ORDER BY arm_cr_wmean ASC NULLS LAST),
		w3 AS (PARTITION BY aliquot_barcode ORDER BY arm_cr_wmean DESC NULLS LAST)
)
SELECT aliquot_barcode,aneuploidy_score,aneuploidy_amp_score,aneuploidy_del_score,max_loss_arm_n,max_loss_arm_wmean,max_loss_arm_wsd,max_gain_arm_n,max_gain_arm_wmean,max_gain_arm_wsd FROM call_aliquot WHERE rn = 1
--SELECT * FROM sort_by_longest_segment