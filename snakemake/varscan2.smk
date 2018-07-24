## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Snakefile for Varscan 2 pipeline
## Authors: Floris Barthel
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run VarScan2 to call SNVs on a tumor/normal pair
## See: http://dkoboldt.github.io/varscan/somatic-calling.html
## Using M2 interval_list (1-based) but converted to bed file (0-based)
## 07/23/2018:
## Added ancient() flag to input because BAM input files are frequently copied and
## timestamps change
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

#VARSCAN_SET = [ "snp.Somatic.hc", "snp.Somatic.lc", "snp.Germline", "snp.LOH" ]

rule varscan:
    input:
        tumor = lambda wildcards: ancient("results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["tumor_aliquot_id"])),
        normal = lambda wildcards: ancient("results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["normal_aliquot_id"])),
        intervalbed = lambda wildcards: "{dir}/{interval}/scattered.bed".format(dir=config["wgs_scatterdir"], interval=wildcards.interval)
    output:
        temp("results/varscan2/vs2-scatter/{pair_id}.{interval}.snp.vcf"),
        temp("results/varscan2/vs2-scatter/{pair_id}.{interval}.indel.vcf")
    params:
        mem = CLUSTER_META["varscan"]["mem"],
        outputprefix = "results/varscan2/vs2-scatter/{pair_id}.{interval}"
    threads:
        CLUSTER_META["varscan"]["ppn"]
    log:
        "logs/varscan/{pair_id}.{interval}.log"
    benchmark:
        "benchmarks/varscan/{pair_id}.{interval}.txt"
    message:
        "Calling SNVs (VarScan2)\n"
        "Pair: {wildcards.pair_id}\n"
        "Interval: {wildcards.interval}"
    shell:
        "TUMOR_MPILEUP=$(printf 'samtools mpileup -q 1 -f {config[reference_fasta]} -l {input.intervalbed} {input.tumor}');"
        "NORMAL_MPILEUP=$(printf 'samtools mpileup -q 1 -f {config[reference_fasta]} -l {input.intervalbed} {input.normal}');"
        "java -Xmx{params.mem}g -Djava.io.tmpdir={config[tempdir]} \
            -jar jar/VarScan.v2.4.3.jar somatic \
            <($NORMAL_MPILEUP) \
            <($TUMOR_MPILEUP) \
            {params.outputprefix} \
            --min-coverage 8 \
            --min-coverage-normal 6 \
            --min-coverage-tumor 8 \
            --min-var-freq 0.10 \
            --min-freq-for-hom 0.75 \
            --tumor-purity 1.0 \
            --strand-filter 1 \
            --somatic-p-value 0.05 \
            --output-vcf 1 \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Fix Varscan header
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Adds contig lines to VS2 header so VCF files can be merged
## 07/23 added awk one-liner to fix IUPAC codes in input VCF modified from:
## http://seqanswers.com/forums/showthread.php?t=39054
## Added sed one-liner to remove empty lines, taken from:
## https://serverfault.com/questions/252921/how-to-remove-empty-blank-lines-in-a-file-in-unix
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule fixvs2header:
    input:
        snp = "results/varscan2/vs2-scatter/{pair_id}.{interval}.snp.vcf",
        indel = "results/varscan2/vs2-scatter/{pair_id}.{interval}.indel.vcf"
    output:
        snp1 = temp("results/varscan2/vs2-fixheader/{pair_id}.{interval}.snp.fixedIUPAC.vcf"),
        indel1 = temp("results/varscan2/vs2-fixheader/{pair_id}.{interval}.indel.fixedIUPAC.vcf"),
        snp2 = temp("results/varscan2/vs2-fixheader/{pair_id}.{interval}.snp.fixed.vcf"),
        indel2 = temp("results/varscan2/vs2-fixheader/{pair_id}.{interval}.indel.fixed.vcf")
    params:
        mem = CLUSTER_META["fixvs2header"]["mem"]
    threads:
        CLUSTER_META["fixvs2header"]["ppn"]
    log:
        "logs/fixvs2header/{pair_id}.{interval}.log"
    benchmark:
        "benchmarks/fixvs2header/{pair_id}.{interval}.txt"
    message:
        "Fixing VCF header\n"
        "Pair: {wildcards.pair_id}\n"
        "Interval: {wildcards.interval}"
    shell:
    	"awk '{{gsub(/\\y[W|K|Y|R|S|M]\\y/,\"N\",$4); OFS = \"\\t\"; print}}' {input.snp} | sed '/^$/d' > {output.snp1}; "
    	"gatk --java-options -Xmx{params.mem}g UpdateVCFSequenceDictionary \
            -V {output.snp1} \
            --source-dictionary {config[reference_dict]} \
            --replace true \
            -O {output.snp2} \
            > {log} 2>&1; "
        "awk '{{gsub(/\\y[W|K|Y|R|S|M]\\y/,\"N\",$4); OFS = \"\\t\"; print}}' {input.indel} | sed '/^$/d' > {output.indel1}; "
        "gatk --java-options -Xmx{params.mem}g UpdateVCFSequenceDictionary \
            -V {output.indel1} \
            --source-dictionary {config[reference_dict]} \
            --replace true \
            -O {output.indel2} \
            > {log} 2>&1"         

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Merge Varscan
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Copied and edited from M2-merge SNV rule
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule mergevarscan:
    input:
        snp = lambda wildcards: expand("results/varscan2/vs2-fixheader/{pair_id}.{interval}.snp.fixed.vcf", pair_id=wildcards.pair_id, interval=WGS_SCATTERLIST),
        indel = lambda wildcards: expand("results/varscan2/vs2-fixheader/{pair_id}.{interval}.indel.fixed.vcf", pair_id=wildcards.pair_id, interval=WGS_SCATTERLIST)
    output:
        snp = protected("results/varscan2/vcf/{pair_id}.snp.vcf.gz"),
        indel = protected("results/varscan2/vcf/{pair_id}.indel.vcf.gz")
    params:
        mem = CLUSTER_META["mergevarscan"]["mem"]
    threads:
        CLUSTER_META["mergevarscan"]["ppn"]
    log:
        "logs/mergevarscan/{pair_id}.log"
    benchmark:
        "benchmarks/mergevarscan/{pair_id}.txt"
    message:
        "Merging VCF files (VS2)\n"
        "Pair: {wildcards.pair_id}"
    run:
        input_snps = " ".join(["-I " + s for s in input['snp']])
        input_indels = " ".join(["-I " + s for s in input['indel']])
        shell("gatk --java-options -Xmx{params.mem}g MergeVcfs \
            {input_snps} \
            -O {output.snp} \
            --CREATE_INDEX false \
            > {log} 2>&1")
        shell("gatk --java-options -Xmx{params.mem}g MergeVcfs \
            {input_indels} \
            -O {output.indel} \
            --CREATE_INDEX false \
            > {log} 2>&1")

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Process Somatic (Varscan)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Isolate Germline/LOH/Somatic calls from output and identifies high confidence calls
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule processsomatic:
    input:
        snp = "results/varscan2/vcf/{pair_id}.snp.vcf.gz",
        indel = "results/varscan2/vcf/{pair_id}.indel.vcf.gz"
    output:
    	tmpsnp = temp("results/varscan2/vs2-processed/{pair_id}/{pair_id}.snp.vcf"),
    	tmpindel = temp("results/varscan2/vs2-processed/{pair_id}/{pair_id}.indel.vcf"),
        snpgermlinehc = temp("results/varscan2/vs2-processed/{pair_id}/{pair_id}.snp.Germline.hc.vcf"),
        snpgermline = temp("results/varscan2/vs2-processed/{pair_id}/{pair_id}.snp.Germline.vcf"),
        snplohhc = temp("results/varscan2/vs2-processed/{pair_id}/{pair_id}.snp.LOH.hc.vcf"),
        snploh = temp("results/varscan2/vs2-processed/{pair_id}/{pair_id}.snp.LOH.vcf"),
        snpsomatichc = temp("results/varscan2/vs2-processed/{pair_id}/{pair_id}.snp.Somatic.hc.vcf"),
        snpsomatic = temp("results/varscan2/vs2-processed/{pair_id}/{pair_id}.snp.Somatic.vcf"),
        indelgermlinehc = temp("results/varscan2/vs2-processed/{pair_id}/{pair_id}.indel.Germline.hc.vcf"),
        indelgermline = temp("results/varscan2/vs2-processed/{pair_id}/{pair_id}.indel.Germline.vcf"),
        indellohhc = temp("results/varscan2/vs2-processed/{pair_id}/{pair_id}.indel.LOH.hc.vcf"),
        indelloh = temp("results/varscan2/vs2-processed/{pair_id}/{pair_id}.indel.LOH.vcf"),
        indelsomatichc = temp("results/varscan2/vs2-processed/{pair_id}/{pair_id}.indel.Somatic.hc.vcf"),
        indelsomatic = temp("results/varscan2/vs2-processed/{pair_id}/{pair_id}.indel.Somatic.vcf")
    params:
        mem = CLUSTER_META["processsomatic"]["mem"]
    threads:
        CLUSTER_META["processsomatic"]["ppn"]
    log:
        "logs/processsomatic/{pair_id}.log"
    benchmark:
        "benchmarks/processsomatic/{pair_id}.txt"
    message:
        "Isolate Germline/LOH/Somatic calls from output and identifies high confidence calls (Varscan2)\n"
        "Pair: {wildcards.pair_id}"
    shell:
    	"zcat {input.snp} > {output.tmpsnp}; "
    	"java -Xmx{params.mem}g -Djava.io.tmpdir={config[tempdir]} \
            -jar jar/VarScan.v2.4.3.jar processSomatic {output.tmpsnp} \
            > {log} 2>&1; "
        "zcat {input.indel} > {output.tmpindel}; "
    	"java -Xmx{params.mem}g -Djava.io.tmpdir={config[tempdir]} \
            -jar jar/VarScan.v2.4.3.jar processSomatic {output.tmpindel} \
            >> {log} 2>&1; "

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Somatic Filter (Varscan)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Filter somatic variants for clusters/indels
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule somaticfilter:
    input:
        snphc = "results/varscan2/vs2-processed/{pair_id}/{pair_id}.snp.Somatic.hc.vcf",
        indel = "results/varscan2/vs2-processed/{pair_id}/{pair_id}.indel.vcf"
    output:
    	"results/varscan2/vs2-filter/{pair_id}.snp.Somatic.hc.filter.vcf"
    params:
        mem = CLUSTER_META["somaticfilter"]["mem"]
    threads:
        CLUSTER_META["somaticfilter"]["ppn"]
    log:
        "logs/somaticfilter/{pair_id}.log"
    benchmark:
        "benchmarks/somaticfilter/{pair_id}.txt"
    message:
        "Filter somatic variants for clusters/indels (Varscan2)\n"
        "Pair: {wildcards.pair_id}"
    shell:
    	"java -Xmx{params.mem}g -Djava.io.tmpdir={config[tempdir]} \
            -jar jar/VarScan.v2.4.3.jar somaticFilter {input.snphc} \
            --indel-file {input.indel} \
            --output-file {output} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## BAM readcount (Varscan2)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Obtain BAM readcounts for a given list of variants (Varscan2)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule bamreadcount:
    input:
        vcf = "results/varscan2/vs2-filter/{pair_id}.snp.Somatic.hc.filter.vcf",
        bam = lambda wildcards: ancient("results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["tumor_aliquot_id"]))
    output:
    	temp("results/varscan2/bam-readcount/{pair_id}.variants.readcounts")
    params:
        mem = CLUSTER_META["bamreadcount"]["mem"]
    threads:
        CLUSTER_META["bamreadcount"]["ppn"]
    log:
        "logs/bamreadcount/{pair_id}.log"
    benchmark:
        "benchmarks/bamreadcount/{pair_id}.txt"
    message:
        "Obtain BAM readcounts for a given list of variants (Varscan2)\n"
        "Pair: {wildcards.pair_id}"
    shell:
    	"awk '{{OFS=\"\\t\"; if (!/^#/){{print $1,$2,$2}}}}' {input.vcf} | \
    		bam-readcount \
    		-q 1 \
    		-b 20 \
    		-f {config[reference_fasta]} \
    		-l /dev/stdin \
    		{input.bam} \
    		2> {log} \
    		1> {output}"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## False positive filter (Varscan2)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Apply the false-positive filter (Varscan2)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule fpfilter:
    input:
        readcounts = "results/varscan2/bam-readcount/{pair_id}.variants.readcounts",
        vcf = "results/varscan2/vs2-filter/{pair_id}.snp.Somatic.hc.filter.vcf"
    output:
    	protected("results/varscan2/fpfilter/{pair_id}.snp.Somatic.hc.final.vcf.gz")
    params:
    	vcfout = "results/varscan2/fpfilter/{pair_id}.snp.Somatic.hc.final.vcf",
        mem = CLUSTER_META["fpfilter"]["mem"]
    threads:
        CLUSTER_META["fpfilter"]["ppn"]
    log:
        "logs/fpfilter/{pair_id}.log"
    benchmark:
        "benchmarks/fpfilter/{pair_id}.txt"
    message:
        "Apply the false-positive filter (Varscan2)\n"
        "Pair: {wildcards.pair_id}"
    shell:
    	"java -Xmx{params.mem}g -Djava.io.tmpdir={config[tempdir]} \
            -jar jar/VarScan.v2.4.3.jar fpfilter \
            {input.vcf} {input.readcounts} \
            --output-file {params.vcfout} \
            > {log} 2>&1; "
        "bgzip -i {params.vcfout} && bcftools index -t {output}"

## END ##