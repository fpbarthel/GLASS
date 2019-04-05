/*
Use GATK-based arm level calls (Taylor method) to call the glioma event combined amplification of chr 7 and chr 10 deletion
*/
WITH
t1 AS
(
	SELECT
		aliquot_barcode,
		chrom,
		arm,
		arm_call,
		(CASE
		 WHEN chrom = 7 THEN arm_call = 1
		 WHEN chrom = 10 THEN arm_call = -1
		 END) AS bool_call
	FROM analysis.gatk_cnv_by_arm
	WHERE chrom IN (7,10)
),
t2 AS
(
	SELECT
		aliquot_barcode,
		chrom,
		COUNT(CASE WHEN bool_call IS TRUE THEN 1 END) AS count_true, -- number of chromosome arms with event
		COUNT(CASE WHEN bool_call IS FALSE THEN 1 END) AS count_false, -- number of chromosome arms lacking event
		COUNT(CASE WHEN bool_call IS NULL THEN 1 END) AS count_null -- number of chromosome arms unknown
	FROM t1
	GROUP BY 1, 2
),
t3 AS
(
	SELECT
		aliquot_barcode,
		(CASE WHEN bool_or(count_null = 2) THEN NULL ELSE bool_and(count_true > 0) AND bool_and(count_false = 0) END) AS c710
	FROM t2
	GROUP BY 1
)
SELECT * FROM t3