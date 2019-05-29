WITH
snv_pairs AS
(
	SELECT t1.gene_symbol AS gene_symbol_a, t2.gene_symbol AS gene_symbol_b, t1.idh_codel_subtype AS idh_codel_subtype 
	FROM ref.snv_drivers_subtype t1
	INNER JOIN ref.snv_drivers_subtype t2 ON t1.idh_codel_subtype = t2.idh_codel_subtype
	WHERE t1.gene_symbol < t2.gene_symbol
	ORDER BY 3,1,2
),
cnv_pairs AS
(
	SELECT t1.gene_symbol AS gene_symbol_a, t1.direction AS direction_a, t2.gene_symbol AS gene_symbol_b, t2.direction AS direction_b, t1.idh_codel_subtype AS idh_codel_subtype
	FROM ref.cnv_drivers_subtype t1
	INNER JOIN ref.cnv_drivers_subtype t2 ON t1.idh_codel_subtype = t2.idh_codel_subtype
	WHERE t1.gene_symbol < t2.gene_symbol
	ORDER BY 5,1,2
),
arm_pairs AS
(
	SELECT t1.arm AS arm_a, t1.direction AS direction_a, t2.arm AS arm_b, t2.direction AS direction_b, t1.idh_codel_subtype AS idh_codel_subtype
	FROM ref.arm_drivers_subtype t1
	INNER JOIN ref.arm_drivers_subtype t2 ON t1.idh_codel_subtype = t2.idh_codel_subtype
	WHERE t1.arm < t2.arm
	ORDER BY 5,1,2
),
pairs AS
(
	SELECT gene_symbol_a  || ' mut' AS evnt_a, gene_symbol_b  || ' mut' AS evnt_b, idh_codel_subtype FROM snv_pairs
	UNION
	SELECT gene_symbol_a || (CASE direction_a WHEN -2 THEN ' del' WHEN 2 THEN ' amp' ELSE NULL END) AS evnt_a, gene_symbol_b || (CASE direction_b WHEN -2 THEN ' del' WHEN 2 THEN ' amp' ELSE NULL END) AS evnt_b, idh_codel_subtype FROM cnv_pairs
	UNION
	SELECT arm_a || (CASE direction_a WHEN -1 THEN ' del' WHEN 1 THEN ' amp' ELSE NULL END) AS evnt_a, arm_b || (CASE direction_b WHEN -1 THEN ' del' WHEN 1 THEN ' amp' ELSE NULL END) AS evnt_b, idh_codel_subtype FROM arm_pairs
)
SELECT * FROM pairs
ORDER BY 3,1,2