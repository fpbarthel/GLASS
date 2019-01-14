WITH
/*
Define a list of 'selected tumor pairs': select a single primary/recurrent pair per subject/case after removing any aliquots from the blocklist and sorting samples by surgical interval,
so that pairs with the highest surgical interval are prioritized, and pairs with WGS are prioritized over WXS when both are available
*/
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
		b1.cnv_exclusion <> 'block' AND b2.cnv_exclusion <> 'block' 
),
selected_genes AS
(
	SELECT gene_symbol::varchar(255)
	FROM (VALUES ('ATM'),('ATR'),('RPA1'),('RPA2'),('RPA3'),('RPA4'),('BRCA1'),('BRCA2'),('RAD51'),('RFC1'),('RFC2'),('RFC3'),('RFC4'),('RFC5'),('XRCC1'),('PCNA'),('PARP1'),('ERCC1'),('MSH3'),('PMS2'),('MLH1'),('MSH6'),('MSH2'),('MLH3'),('EXO1'),('NBN'),('RAD50'),('CHEK2'),('FANCI'),('FANCD2'),('FANCA'),('FANCC'),('FANCE'),('FANCL'),('FANCG'),('FANCM'),('ERCC4'),('ERCC2'),('ERCC5'),('PARP2'),('APEX1'),('FEN1'),('XPC'),('ERCC6'),('GTF2H2'),('ERCC3'),('XPA'),('RAD23B'),('PALB2'),('RAD51C'),('XRCC6'),('XRCC5'),('PRKDC'),('XRCC4'),('Lig4'),('FANCB'),('FANCF'),('FAAP24'),('CHEK1'),('BRIP1'),('SLX4'),('FAN1'),('MUS81'),('EME1'),('POLE'),('POLD1'),('MRE11A'),('RAD51D'),('RAD52'),('RAD51B'),('PMS1'),('RAD23A'),('LIG3'),('MGMT'),('OGG1'),('UNG'),('SMUG1'),('MBD4'),('TDG'),('MUTYH'),('NTHL1'),('MPG'),('NEIL1'),('NEIL2'),('NEIL3'),('APEX2'),('PNKP'),('APLF'),('PARP3'),('ALKBH2'),('ALKBH3'),('MSH4'),('MSH5'),('PMS2P3'),('CETN2'),('DDB1'),('DDB2'),('GTF2H1'),('GTF2H3'),('GTF2H4'),('GTF2H5'),('CDK7'),('CCNH'),('MNAT1'),('LIG1'),('ERCC8'),('UVSSA'),('XAB2'),('MMS19'),('DMC1'),('XRCC2'),('XRCC3'),('RAD54L'),('RAD54B'),('SHFM1'),('RBBP8'),('SLX1A'),('SLX1B'),('GEN1'),('FAAP20'),('DCLRE1C'),('NHEJ1'),('PAXIP1'),('BLM'),('MLL3'),('CRIP1'),('CDK12'),('BAP1'),('BARD1'),('WRN'),('BUB1'),('CENPE'),('ZW10'),('TTK'),('KNTC1'),('AURKB'),('POLB'),('POLH'),('POLQ'),('TDP1'),('TDP2'),('NUDT1'),('DUT'),('RRM2B'),('POLG'),('REV3L'),('MAD2L2'),('REV1'),('POLI'),('POLK'),('POLL'),('POLM'),('POLN'),('TREX1'),('TREX2'),('APTX'),('SPO11'),('ENDOV'),('UBE2A'),('UBE2B'),('RAD18'),('SHPRH'),('HLTF'),('RNF168'),('SPRTN'),('RNF8'),('RNF4'),('UBE2V2'),('UBE2N'),('H2AFX'),('CHAF1A'),('SETMAR'),('RECQL4'),('MPLKIP'),('DCLRE1A'),('DCLRE1B'),('PRPF19'),('RECQL'),('RECQL5'),('HELQ'),('RDM1'),('NABP2'),('ATRIP'),('MDC1'),('RAD1'),('RAD9A'),('HUS1'),('RAD17'),('TP53'),('TP53BP1'),('TOPBP1'),('CLK2'),('PER1')) AS v (gene_symbol)
),
selected_genes_samples AS
(
	SELECT * FROM selected_tumor_pairs, selected_genes WHERE priority = 1
),
cnv_by_pair_gene AS
(
	SELECT
		sgs.case_barcode,
		sgs.gene_symbol,
		c1.cn_call AS cn_a,
		c2.cn_call AS cn_b
	FROM selected_genes_samples sgs
	LEFT JOIN analysis.cnv_by_gene_gatk c1 ON c1.aliquot_barcode = sgs.tumor_barcode_a AND c1.gene_symbol = sgs.gene_symbol
	LEFT JOIN analysis.cnv_by_gene_gatk c2 ON c2.aliquot_barcode = sgs.tumor_barcode_b AND c2.gene_symbol = sgs.gene_symbol
)
SELECT case_barcode, gene_symbol, cnv_retention, cnv_class
FROM cnv_by_pair_gene cn
LEFT JOIN analysis.gatk_seg_retention_states gs ON gs.cn_a = cn.cn_a AND gs.cn_b = cn.cn_b