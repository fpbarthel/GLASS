/*
Build a squared gene x subject matrix for heatmap plotting
1. Get selected tumor pairs as described elsewhere (`selected_tumor_pairs`)
	- reduce this list to a list of aliquots, combining (a) and (b) into a single column (`selected_aliquots`). This list is only used for coverage
2. Select a list of genes and specific variants
	- Take genes deemed significant using the dNdS-CV method
	- Manually adding a genes that are not significant by dNdS but are known glioma genes (handpicked from Anzhela's list)
	- Filter the list of variants by subsetting hotspot variants only for those genes where they are known
3. Take the cartesian product between the genes and aliquots table, generating a table with (genes x aliquots) numbers of rows (`selected_genes_samples`)
	- Get coverage statistics from all these loci across all aliquots and compute median for each gene and aliquot (`gene_sample_coverage`)
	- Use the tumor pairs list from step 1 to align tumor pairs to genes and coverage statistics for sample (a) and (b) from each pair (`gene_pair_coverage`)
4. Take variants from each tumor_barcode, subsetting selected tumor pairs from step 1 and by variants from step 2 (`variants_by_case_and_gene`)
	- For each pair, there exists an optimal first tumor (a) and subsequent tumor sample (b)
	- Aggregate over variants by case_barcode and gene_symbol, selecting only the top variant for each subject/gene combination
	- Restrict to events with coverage >= 5 in both (a) and (b)
	- Variants for each tumor pair/gene combination are ordered according to variant_classification_priority (see new table analysis.variant_classifications) and whether or not the mutation was called in a/b and finally based on read_depth
	- Manually calling all variants based on a 0.05 VAF threshold, since we already applied stringent filtering prior to this step (eg. only variants with a mutect call in either (a) or (b) are included in the data at this point)
5. Right join the limited list of variants from step 4 with the extensive coverage statistics from step 3, generating a squared table with all genes and all subjects (`squared_variants`)
6. Perform a subject level call of whether the variant was called in sample A only (initial only), B only (recurrence only), or both (shared), or in neither mark as wild-type.
	- If there is not enough coverage in either A or B, mark the gene/subject combination as NULL to indicate insufficient coverage.
*/
WITH
selected_tumor_pairs AS
(
	SELECT
		tumor_pair_barcode,
		case_barcode,
		tumor_barcode_a,
		tumor_barcode_b,
		row_number() OVER (PARTITION BY case_barcode ORDER BY surgical_interval_mo DESC, portion_a ASC, portion_b ASC, substring(tumor_pair_barcode from 27 for 3) ASC) AS priority
	FROM analysis.tumor_pairs ps
	LEFT JOIN analysis.blocklist b1 ON b1.aliquot_barcode = ps.tumor_barcode_a
	LEFT JOIN analysis.blocklist b2 ON b2.aliquot_barcode = ps.tumor_barcode_b
	WHERE
		comparison_type = 'longitudinal' AND
		sample_type_b <> 'M1' AND 													-- exclude metastatic samples here because this is outside the scope of our study
		b1.coverage_exclusion = 'allow' AND b2.coverage_exclusion = 'allow' 
),
selected_aliquots AS
(
	SELECT tumor_barcode_a AS aliquot_barcode, case_barcode FROM selected_tumor_pairs WHERE priority = 1
	UNION
	SELECT tumor_barcode_b AS aliquot_barcode, case_barcode FROM selected_tumor_pairs WHERE priority = 1
),
selected_genes AS
(
	SELECT sn.gene_symbol, chrom, pos, alt, sn.variant_classification, variant_classification_priority, hgvs_p
	FROM analysis.snvs sn
	LEFT JOIN analysis.variant_classifications vc ON sn.variant_classification = vc.variant_classification
	WHERE
		variant_classification_priority IS NOT NULL AND 
		(sn.gene_symbol IN ('ATM','ATR','RPA1','RPA2','RPA3','RPA4','BRCA1','BRCA2','RAD51','RFC1','RFC2','RFC3','RFC4','RFC5','XRCC1','PCNA','PARP1','ERCC1','MSH3','PMS2','MLH1','MSH6','MSH2','MLH3','EXO1','NBN','RAD50','CHEK2','FANCI','FANCD2','FANCA','FANCC','FANCE','FANCL','FANCG','FANCM','ERCC4','ERCC2','ERCC5','PARP2','APEX1','FEN1','XPC','ERCC6','GTF2H2','ERCC3','XPA','RAD23B','PALB2','RAD51C','XRCC6','XRCC5','PRKDC','XRCC4','Lig4','FANCB','FANCF','FAAP24','CHEK1','BRIP1','SLX4','FAN1','MUS81','EME1','POLE','POLD1','MRE11A','RAD51D','RAD52','RAD51B','PMS1','RAD23A','LIG3','MGMT','OGG1','UNG','SMUG1','MBD4','TDG','MUTYH','NTHL1','MPG','NEIL1','NEIL2','NEIL3','APEX2','PNKP','APLF','PARP3','ALKBH2','ALKBH3','MSH4','MSH5','PMS2P3','CETN2','DDB1','DDB2','GTF2H1','GTF2H3','GTF2H4','GTF2H5','CDK7','CCNH','MNAT1','LIG1','ERCC8','UVSSA','XAB2','MMS19','DMC1','XRCC2','XRCC3','RAD54L','RAD54B','SHFM1','RBBP8','SLX1A','SLX1B','GEN1','FAAP20','DCLRE1C','NHEJ1','PAXIP1','BLM','MLL3','CRIP1','CDK12','BAP1','BARD1','WRN','BUB1','CENPE','ZW10','TTK','KNTC1','AURKB','POLB','POLH','POLQ','TDP1','TDP2','NUDT1','DUT','RRM2B','POLG','REV3L','MAD2L2','REV1','POLI','POLK','POLL','POLM','POLN','TREX1','TREX2','APTX','SPO11','ENDOV','UBE2A','UBE2B','RAD18','SHPRH','HLTF','RNF168','SPRTN','RNF8','RNF4','UBE2V2','UBE2N','H2AFX','CHAF1A','SETMAR','RECQL4','MPLKIP','DCLRE1A','DCLRE1B','PRPF19','RECQL','RECQL5','HELQ','RDM1','NABP2','ATRIP','MDC1','RAD1','RAD9A','HUS1','RAD17','TP53','TP53BP1','TOPBP1','CLK2','PER1'))
),
selected_genes_samples AS
(
	SELECT aliquot_barcode, case_barcode, gene_symbol, chrom, lower(pos) AS start_pos, upper(pos)-1 AS end_pos, alt
	FROM selected_aliquots, selected_genes
),
variants_by_case_and_gene AS
(
	SELECT
		gtc.gene_symbol,
		gtc.case_barcode,
		gtc.tumor_pair_barcode,
		gtc.chrom,
		gtc.pos,
		gtc.alt,
		gtc.variant_classification,
		sg.hgvs_p,
		mutect2_call_a AS selected_call_a, --(alt_count_a::decimal / (alt_count_a+ref_count_a) > 0.05) AS selected_call_a,
		mutect2_call_b AS selected_call_b, --(alt_count_b::decimal / (alt_count_b+ref_count_b) > 0.05) AS selected_call_b,
		row_number() OVER (PARTITION BY gtc.gene_symbol, gtc.case_barcode ORDER BY mutect2_call_a::integer + mutect2_call_b::integer = 2, variant_classification_priority, mutect2_call_a::integer + mutect2_call_b::integer DESC, read_depth_a + read_depth_b DESC) AS priority
	FROM analysis.master_genotype_comparison gtc
	INNER JOIN selected_tumor_pairs stp ON stp.tumor_pair_barcode = gtc.tumor_pair_barcode AND stp.priority = 1
	INNER JOIN selected_genes sg ON sg.chrom = gtc.chrom AND sg.pos = gtc.pos AND sg.alt = gtc.alt
	WHERE
		(mutect2_call_a OR mutect2_call_b) AND read_depth_a >= 5 AND read_depth_b >= 5
),
squared_variants AS
(
	SELECT
		vcg.gene_symbol,
		vcg.case_barcode,
		vcg.chrom,
		vcg.pos,
		vcg.alt,
		vcg.variant_classification,
		vcg.hgvs_p,
		vcg.selected_call_a,
		vcg.selected_call_b
	FROM (SELECT * FROM variants_by_case_and_gene WHERE priority = 1) vcg
)
SELECT
	gene_symbol,
	var.case_barcode,
	variant_classification,
	hgvs_p,
	selected_call_a,
	selected_call_b,
	(CASE
	 WHEN selected_call_a AND selected_call_b 				THEN 'Shared'
	 WHEN selected_call_a AND NOT selected_call_b 			THEN 'Shed'
	 WHEN selected_call_b AND NOT selected_call_a 			THEN 'Acquired'
	 ELSE NULL END) AS variant_call
FROM squared_variants var