## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## pVACseq pipeline for calling neoantigens
## Author: Frederick Varn
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

LENGTHS = [8,9,10,11]

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Select Mutect2 calls that have passed filters:
## Filter 1: PASSED calls
## Filter 2: Force-called hotspot mutations (IDH, TERT)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule vcf_pass:
    input:
        "results/mutect2/dropgt/{case_barcode}.filtered.normalized.sorted.vcf.gz"
    output:
        temp("results/pvacseq/passed/{case_barcode}.filtered.normalized.sorted.passed.vcf")
    params:
    	mem = CLUSTER_META["vcf_pass"]["mem"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/select/{case_barcode}.pass.log"
    message:
        "Filtering VCF for PASSed calls \n"
        "Sample: {wildcards.case_barcode}"
    shell:
    	"(gatk --java-options -Xmx{params.mem}g SelectVariants \
    	--exclude-filtered TRUE \
    	-V {input} \
    	-O {output}) 2>{log}"

rule vcf_hotspot:
    input:
        "results/mutect2/dropgt/{case_barcode}.filtered.normalized.sorted.vcf.gz"
    output:
        temp("results/pvacseq/passed/{case_barcode}.filtered.normalized.sorted.hotspot.vcf")
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/select/{case_barcode}.hotspot.log"
    message:
        "Filtering VCF for force-called hotspot mutations \n"
        "Sample: {wildcards.case_barcode}"
    shell:
    	"(gatk SelectVariants \
    	-L 2:209113112 \
    	-L 2:209113113 \
    	-L 5:1295169 \
    	-L 5:1295228 \
    	-L 5:1295242 \
    	-L 5:1295250 \
    	-L 15:90631837 \
    	-L 15:90631838 \
    	-L 15:90631839 \
    	-V {input} \
    	-O {output}) 2>{log}"

rule vcf_merge:
    input:
        I1 = "results/pvacseq/passed/{case_barcode}.filtered.normalized.sorted.passed.vcf",
        I2 = "results/pvacseq/passed/{case_barcode}.filtered.normalized.sorted.hotspot.vcf"
    output:
        "results/pvacseq/passed/{case_barcode}.filtered.normalized.sorted.select.vcf"
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/select/{case_barcode}.merge.log"
    message:
        "Merging filtered VCF files \n"
        "Sample: {wildcards.case_barcode}"
    shell: 
    	"(gatk MergeVcfs \
    	-I {input.I1} \
    	-I {input.I2} \
    	-O {output}) 2>{log}"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run VEP with Downstream and Wildtype plugins
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule vep_plugins:
    input:
        "results/pvacseq/passed/{case_barcode}.filtered.normalized.sorted.select.vcf"
    output:
        "results/pvacseq/vep/{case_barcode}.filtered.normalized.sorted.select.vep.vcf"
    params:
        mem = 8
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
        "results/pvacseq/vep/{case_barcode}.filtered.normalized.sorted.select.vep.vcf"
    output:
        "results/pvacseq/vep/{case_barcode}.filtered.normalized.sorted.select.vep.onesamp.vcf"
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
        vcf = "results/pvacseq/vep/{case_barcode}.filtered.normalized.sorted.select.vep.onesamp.vcf",
        hla = lambda wildcards: expand("results/optitype/HLA_calls/{normals}/{normals}_result.tsv", normals = manifest.getNormalByCase(wildcards.case_barcode))
    output:
        "results/pvacseq/neoag_frag/{case_barcode}_{epitope_lengths}/MHC_Class_I/{case_barcode}_{epitope_lengths}.final.tsv"
    params:
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
    	allele=`tail -n+2 -q {input.hla} | cat | rev | sort -k2r | cut -f3- | rev | cut -f2- | head -1 | awk '{{for(i=1;i<=NF;i++)sub("^", "HLA-", $i)}}; 1' | sed 'y/ /,/'`
		(pvacseq run \
		{input.vcf} \
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
