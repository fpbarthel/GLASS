## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## pVACseq pipeline for calling neoantigens
## Author: Frederick Varn
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

LENGTHS = [8,9,10,11]

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Decompress zipped VCF files
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule gunzip:
    input:
        ancient("results/mutect2/m2filter/{case_barcode}.filtered.vcf.gz")
    output:
        temp("results/pvacseq/gunzip/{case_barcode}.filtered.vcf")
    log:
        "logs/gunzip/{case_barcode}.log"
    message:
        "Decompressing zipped VCF files \n"
        "Sample: {wildcards.case_barcode}"
    shell:
        "(gunzip -c {input} > {output}) 2>{log}"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Select Mutect2 calls that have passed filters
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule vcf_pass:
    input:
        "results/pvacseq/gunzip/{case_barcode}.filtered.vcf"
    output:
        temp("results/pvacseq/passed/{case_barcode}.filtered.passed.vcf")
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/passed/{case_barcode}.log"
    message:
        "Filtering VCF \n"
        "Sample: {wildcards.case_barcode}"
    shell:
    	"(gatk SelectVariants \
    	--exclude-filtered TRUE \
    	-V {input} \
    	-O {output}) 2>{log}"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run VEP with Downstream and Wildtype plugins
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule vep_plugins:
    input:
        "results/pvacseq/passed/{case_barcode}.filtered.passed.vcf"
    output:
        "results/pvacseq/vep/{case_barcode}.filtered.passed.vep.vcf"
    params:
        mem = 4
    conda:
        "../envs/vep.yaml"
    log:
        "logs/vep/{case_barcode}.log"
    message:
        "Running VEP with Downstream and Wildtype plugins \n"
        "Sample: {wildcards.case_barcode}"
    shell:
        "(vep \
		-i {input} \
		-o {output} \
		--cache \
		--offline \
		--assembly GRCh37 \
		--dir_cache {config[vepplugindata]} \
		--format vcf \
		--vcf \
		--symbol \
		--pick \
		--fork {params.mem} \
		--plugin Downstream \
		--plugin Wildtype \
		--terms SO \
		--dir_plugins {config[vepplugins]}) 2>{log}"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Make each VCF contain only one sample (need to do this for pVACseq compatibility)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule vcf_extract:
    input:
        "results/pvacseq/vep/{case_barcode}.filtered.passed.vep.vcf"
    output:
        "results/pvacseq/vep/{case_barcode}.filtered.passed.vep.onesamp.vcf"
    params:
    	sample = lambda wildcards: manifest.getTumorByCase(wildcards.case_barcode)[0]
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/vcf_extract/{case_barcode}.log"
    message:
        "Unifying VCF file \n"
        "Sample: {wildcards.case_barcode}"
    shell:
        "(bcftools view --samples {params.sample} {input} > {output}) 2>{log}"
		
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Call neoantigens with pVACseq
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule pvacseq:
    input:
        "results/pvacseq/vep/{case_barcode}.filtered.passed.vep.onesamp.vcf"
    output:
        "results/pvacseq/neoag_frag/{case_barcode}_{epitope_lengths}/MHC_Class_I/{case_barcode}_{epitope_lengths}.final.tsv"
    params:
    	normal				= lambda wildcards: manifest.getNormalByCase(wildcards.case_barcode)[0],
    	sample_name 		= "{case_barcode}",
        algorithm 			= "NetMHCpan",
    	output_dir 			= "results/pvacseq/neoag_frag/{case_barcode}_{epitope_lengths}/",
        min_fold_change 	= 1,
        top_score_metric 	= "lowest",
        downstream_length 	= "full"
    conda:
        "../envs/pvacseq.yaml"
    log:
        "logs/pvacseq/{case_barcode}.{epitope_lengths}.log"
    message:
        "Calling neoantigens with pVACseq \n"
        "Sample: {wildcards.case_barcode} \n"
        "Epitope length: {wildcards.epitope_lengths}"
    shell:
    	"""
    	module load python/2.7.13
    	allele=`cut -f2- results/optitype/HLA_calls/{params.normal}/{params.normal}_result.tsv | rev | cut -f3- | rev | sed -n '2p' | awk '{{for(i=1;i<=NF;i++)sub("^", "HLA-", $i)}}; 1' | sed 'y/ /,/'`
		(pvacseq run \
		{input} \
		{params.sample_name}_{wildcards.epitope_lengths} \
		$allele \
		{params.algorithm} \
		{params.output_dir} \
		-e {wildcards.epitope_lengths} \
		--iedb-install-directory {config[iedb]} \
		-c {params.min_fold_change} \
		--top-score-metric={params.top_score_metric} \
		-d {params.downstream_length}) 2>{log}
		"""
		
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Merge pVACseq outputs into one file
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule mergepvac:
    input:
        expand("results/pvacseq/neoag_frag/{{case_barcode}}_{epitope_lengths}/MHC_Class_I/{{case_barcode}}_{epitope_lengths}.final.tsv", epitope_lengths = LENGTHS)
    output:
        protected("results/pvacseq/neoantigens/{case_barcode}/MHC_Class_I/{case_barcode}.final.tsv")
    log:
        "logs/merge/{case_barcode}.log"
    message:
        "Merging pVACseq output \n"
        "Sample: {wildcards.case_barcode}"
    shell:
    	"""
    	set +o pipefail; cat {input} | head -1 > {output}
    	tail -n+2 -q {input} >> {output} 2>{log}
    	"""
