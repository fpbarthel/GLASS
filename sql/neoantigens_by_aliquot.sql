WITH var_neo AS
(
	SELECT
		tumor_barcode::character(29) AS case_barcode,
		snv.variant_id,
		neo.chrom::integer,
		snv.pos,
		snv.alt::character varying(510),
		neo.transcript::character varying(255),
		ensembl_gene_id::character varying(255),
		neo.variant_type::character varying(255) AS pvacseq_variant_type,
		mutation::character varying(255),
		neo.protein_position::character varying(255) AS pvacseq_protein_position,
		gene_name::character varying(255),
		hla_allele::character varying(255),
		peptide_length::integer,
		sub_peptide_position::integer,
		mutation_position::integer,
		mt_epitope_seq::character varying(255),
		wt_epitope_seq::character varying(255),
		netmhcpan_mt_score::double precision,
		netmhcpan_wt_score::double precision,
		netmhcpan_fold_change::double precision
	FROM analysis.neoantigens_by_patient neo
	LEFT JOIN variants.anno snv ON snv.chrom = neo.chrom AND snv.pos = int4range(neo.stop,neo.stop+1) AND snv.alt = neo.alt AND snv.ref = neo.ref
	WHERE snv.variant_type IS NOT NULL AND snv.variant_type = 'SNP'

	UNION

	SELECT
		tumor_barcode::character(29),
		snv.variant_id,
		neo.chrom::integer,
		snv.pos,
		snv.alt::character varying(510),
		neo.transcript::character varying(255),
		ensembl_gene_id::character varying(255),
		neo.variant_type::character varying(255) AS pvacseq_variant_type,
		mutation::character varying(255),
		neo.protein_position::character varying(255) AS pvacseq_protein_position,
		gene_name::character varying(255),
		hla_allele::character varying(255),
		peptide_length::integer,
		sub_peptide_position::integer,
		mutation_position::integer,
		mt_epitope_seq::character varying(255),
		wt_epitope_seq::character varying(255),
		netmhcpan_mt_score::double precision,
		netmhcpan_wt_score::double precision,
		netmhcpan_fold_change::double precision
	FROM analysis.neoantigens_by_patient neo
	LEFT JOIN variants.anno snv ON snv.chrom = neo.chrom AND snv.pos = int4range(neo.start,neo.start+(length(neo.ref)-1),'[]') AND snv.alt = neo.alt AND snv.ref = neo.ref
	WHERE snv.variant_type IS NOT NULL AND snv.variant_type IN ('INS','DEL','DNP','TNP')
),
expanded AS
(
	SELECT al.aliquot_barcode,vn.*
	FROM var_neo vn
	LEFT JOIN biospecimen.aliquots al ON substr(al.sample_barcode,1,12) = vn.case_barcode
	ORDER BY case_barcode,vn.chrom,vn.pos,hla_allele,peptide_length,sub_peptide_position,aliquot_barcode
),
aliquot_var AS
(
	SELECT ex.*,
	pg.ad_ref,pg.ad_alt,pg.af,pg.ssm2_pass_call,pg.ssm2_saaf_none,pg.major_cn,pg.minor_cn,pg.titan_cluster,pg.titan_ccf,pg.pyclone_cluster,pg.pyclone_ccf,pg.pyclone_vaf
	FROM expanded ex
	INNER JOIN variants.passgeno pg ON pg.aliquot_barcode = ex.aliquot_barcode AND pg.variant_id = ex.variant_id
	ORDER BY ex.case_barcode,ex.chrom,ex.pos,hla_allele,peptide_length,sub_peptide_position,aliquot_barcode
)
SELECT *
INTO analysis.neoantigens_by_aliquot
FROM aliquot_var

