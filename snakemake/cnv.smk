## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Collect read counts
## See: https://software.broadinstitute.org/gatk/documentation/article?id=11682
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

def selectIntervals(is_exome, chromosomes):
    if is_exome:
        if chromosomes == "autosomal":
            return config["cnv"]["intervals_wes_autosomal"]
        elif chromosomes == "allosomal":
            return config["cnv"]["intervals_wes_allosomal"]
    else:
        if chromosomes == "autosomal":
            return config["cnv"]["intervals_wgs_autosomal"]
        elif chromosomes == "allosomal":
            return config["cnv"]["intervals_wgs_allosomal"]

rule collectreadcounts:
    input:
        bam = ancient("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam")
    output:
        "results/cnv/readcounts/{aliquot_barcode}.{chromosomes}.counts.hdf5"
    params:
        mem = CLUSTER_META["collectreadcounts"]["mem"],
        intervals = lambda wildcards: selectIntervals(manifest.isExome(wildcards.aliquot_barcode), wildcards["chromosomes"])
    threads:
        CLUSTER_META["collectreadcounts"]["ppn"]
    log:
        "logs/cnv/readcounts/{aliquot_barcode}.{chromosomes}.log"
    conda:
        "../envs/gatk4.yaml"
    benchmark:
        "benchmarks/cnv/readcounts/{aliquot_barcode}.{chromosomes}.txt"
    message:
        "Collecting read counts\n"
        "Sample: {wildcards.aliquot_barcode}\n"
        "Chromosomes: {wildcards.chromosomes}"
    shell:
        "gatk --java-options -Xmx{params.mem}g CollectReadCounts \
            -I {input.bam} \
            -L {params.intervals} \
            --interval-merging-rule OVERLAPPING_ONLY \
            -O {output} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Collect Allelic Counts
## Counts reference and alternative alleles at common germline variant sites
## See: https://gatkforums.broadinstitute.org/dsde/discussion/11683/
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule collecthets:
    input:
        ancient("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam")
    output:
        "results/cnv/hets/{aliquot_barcode}.allelicCounts.tsv"
    params:
        mem = CLUSTER_META["collecthets"]["mem"]
    threads:
        CLUSTER_META["collecthets"]["ppn"]
    log:
        "logs/cnv/hets/{aliquot_barcode}.log"
    conda:
        "../envs/gatk4.yaml"
    benchmark:
        "benchmarks/cnv/hets/{aliquot_barcode}.txt"
    message:
        "Collect allelic counts\n"
        "Aliquot: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g CollectAllelicCounts \
            -I {input} \
            -L {config[cnv][hetsites]} \
            -R {config[reference_fasta]} \
            -O {output} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Create CNV panel of normals
## See: https://software.broadinstitute.org/gatk/documentation/article?id=11682
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

def selectGCAnnotation(is_exome, chromosomes):
    if is_exome:
        if chromosomes == "autosomal":
            return config["cnv"]["intervals_wes_autosomal_gc"]
        elif chromosomes == "allosomal":
            return config["cnv"]["intervals_wes_allosomal_gc"]
    else:
        if chromosomes == "autosomal":
            return config["cnv"]["intervals_wgs_autosomal_gc"]
        elif chromosomes == "allosomal":
            return config["cnv"]["intervals_wgs_allosomal_gc"]

rule createcnvpon:
    input:
        lambda wildcards: ancient(expand("results/cnv/readcounts/{aliquot_barcode}.{chromosomes}.counts.hdf5", aliquot_barcode = manifest.getPONAliquotsByBatchAndSex(wildcards.aliquot_batch, wildcards.sex if wildcards.sex != "NA" else None), chromosomes = wildcards.chromosomes))
    output:
        "results/cnv/createcnvpon/{aliquot_batch}.{chromosomes}.{sex}.pon.hdf5"
    params:
        mem = CLUSTER_META["createcnvpon"]["mem"],
        input_files = lambda _, input: " ".join(["-I " + s for s in input]),
        gcanno = lambda wildcards: selectGCAnnotation(wildcards.aliquot_batch.endswith("WXS"), wildcards["chromosomes"])
    threads:
        CLUSTER_META["createcnvpon"]["ppn"]
    log:
        "logs/cnv/createcnvpon/{aliquot_batch}.{chromosomes}.{sex}.log"
    conda:
        "../envs/gatk4.yaml"
    benchmark:
        "benchmarks/cnv/createcnvpon/{aliquot_batch}.{chromosomes}.{sex}.txt"
    message:
        "Creating CNV panel of normals\n"
        "Batch: {wildcards.aliquot_batch}\n"
        "Chromosomes: {wildcards.chromosomes}\n"
        "Sex: {wildcards.sex}"
    shell:
        "gatk --java-options -Xmx{params.mem}g CreateReadCountPanelOfNormals \
            {params.input_files} \
            --annotated-intervals {params.gcanno} \
            -O {output} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Denoise read counts
## See: https://software.broadinstitute.org/gatk/documentation/article?id=11682
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

def getSexString(manifest, aliquot_barcode, chromosomes):
    if chromosomes == "autosomal":
        return "all"
    elif chromosomes == "allosomal":
        sex = manifest.getSex(aliquot_barcode)
        return sex if sex is not None else "NA"

rule denoisereadcounts:
    input:
        sample = lambda wildcards: ancient("results/cnv/readcounts/{aliquot_barcode}.{chromosomes}.counts.hdf5".format(aliquot_barcode = wildcards.aliquot_barcode, chromosomes = wildcards.chromosomes)),
        pon =  lambda wildcards: ancient("results/cnv/createcnvpon/{aliquot_batch}.{chromosomes}.{sex}.pon.hdf5".format(aliquot_batch = manifest.getBatch(wildcards.aliquot_barcode), chromosomes = wildcards.chromosomes, sex = getSexString(manifest, wildcards.aliquot_barcode, wildcards.chromosomes)))
    output:
        standardized = "results/cnv/denoisereadcounts/{aliquot_barcode}.{chromosomes}.standardizedCR.tsv",
        denoised = "results/cnv/denoisereadcounts/{aliquot_barcode}.{chromosomes}.denoisedCR.tsv"
    params:
        mem = CLUSTER_META["denoisereadcounts"]["mem"]
    threads:
        CLUSTER_META["denoisereadcounts"]["ppn"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/cnv/denoisereadcounts/{aliquot_barcode}.{chromosomes}.log"
    benchmark:
        "benchmarks/cnv/denoisereadcounts/{aliquot_barcode}.{chromosomes}.txt"
    message:
        "Denoising and standardizing read counts\n"
        "Aliquot: {wildcards.aliquot_barcode}\n"
        "Chromosomes: {wildcards.chromosomes}"
    shell:
        "gatk --java-options -Xmx{params.mem}g DenoiseReadCounts \
            -I {input.sample} \
            --count-panel-of-normals {input.pon} \
            --standardized-copy-ratios {output.standardized} \
            --denoised-copy-ratios {output.denoised} \
            > {log} 2>&1"


rule mergedenoisedreadcounts:
    input:
        standardized_autosomal =    "results/cnv/denoisereadcounts/{aliquot_barcode}.autosomal.standardizedCR.tsv",
        denoised_autosomal =        "results/cnv/denoisereadcounts/{aliquot_barcode}.autosomal.denoisedCR.tsv",
        standardized_allosomal =    "results/cnv/denoisereadcounts/{aliquot_barcode}.allosomal.standardizedCR.tsv",
        denoised_allosomal =        "results/cnv/denoisereadcounts/{aliquot_barcode}.allosomal.denoisedCR.tsv"
    output:
        standardized =  "results/cnv/mergedenoisedreadcounts/{aliquot_barcode}.standardizedCR.tsv",
        denoised =      "results/cnv/mergedenoisedreadcounts/{aliquot_barcode}.denoisedCR.tsv"
    params:
        mem = CLUSTER_META["mergedenoisedreadcounts"]["mem"]
    threads:
        CLUSTER_META["mergedenoisedreadcounts"]["ppn"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/cnv/mergedenoisedreadcounts/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/cnv/mergedenoisedreadcounts/{aliquot_barcode}.txt"
    message:
        "Merging denoising and standardizing read counts\n"
        "Aliquot: {wildcards.aliquot_barcode}"
    shell:
        "awk 'FNR==1 && NR!=1 {{ while (/^@|^CONTIG/) getline; }} 1 {{print}}' \
            {input.standardized_autosomal} {input.standardized_allosomal} \
            1> {output.standardized} \
            2> {log};"
        "awk 'FNR==1 && NR!=1 {{ while (/^@|^CONTIG/) getline; }} 1 {{print}}' \
            {input.denoised_autosomal} {input.denoised_allosomal} \
            1> {output.denoised} \
            2>> {log};"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Plot denoised and standardized copy ratios
## See: https://software.broadinstitute.org/gatk/documentation/article?id=11682
## NB need to optimize minimum contig length parameter for b37
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule plotcr:
    input:
        standardized = "results/cnv/mergedenoisedreadcounts/{aliquot_barcode}.standardizedCR.tsv",
        denoised = "results/cnv/mergedenoisedreadcounts/{aliquot_barcode}.denoisedCR.tsv"
    output:
        "results/cnv/plotcr/{aliquot_barcode}/{aliquot_barcode}.denoised.png",
        "results/cnv/plotcr/{aliquot_barcode}/{aliquot_barcode}.denoisedLimit4.png",
        "results/cnv/plotcr/{aliquot_barcode}/{aliquot_barcode}.standardizedMAD.txt",
        "results/cnv/plotcr/{aliquot_barcode}/{aliquot_barcode}.denoisedMAD.txt",
        "results/cnv/plotcr/{aliquot_barcode}/{aliquot_barcode}.deltaMAD.txt",
        "results/cnv/plotcr/{aliquot_barcode}/{aliquot_barcode}.scaledDeltaMAD.txt"
    params:
        mem = CLUSTER_META["plotcr"]["mem"],
        outputdir = "results/cnv/plotcr/{aliquot_barcode}",
        outputprefix = "{aliquot_barcode}"
    threads:
        CLUSTER_META["plotcr"]["ppn"]
    log:
        "logs/cnv/plotcr/{aliquot_barcode}.log"
    conda:
        "../envs/gatk4.yaml"
    benchmark:
        "benchmarks/cnv/plotcr/{aliquot_barcode}.txt"
    message:
        "Plot denoised and standardized read counts\n"
        "Aliquot: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g PlotDenoisedCopyRatios \
            --standardized-copy-ratios {input.standardized} \
            --denoised-copy-ratios {input.denoised} \
            --sequence-dictionary {config[reference_dict]} \
            --minimum-contig-length {config[cnv][min_contig_len]} \
            --output {params.outputdir} \
            --output-prefix {params.outputprefix} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Model segments
## Use a guassian-kernel binary-segmentation algorithm to group contiguouis copy ratios into segments
## See: https://gatkforums.broadinstitute.org/dsde/discussion/11683/
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule modelsegments:
    input:
        denoised = "results/cnv/mergedenoisedreadcounts/{aliquot_barcode}.denoisedCR.tsv",
        counts = ancient("results/cnv/hets/{aliquot_barcode}.allelicCounts.tsv")
    output:
        "results/cnv/modelsegments/{aliquot_barcode}/{aliquot_barcode}.modelBegin.seg",
        "results/cnv/modelsegments/{aliquot_barcode}/{aliquot_barcode}.modelFinal.seg",
        "results/cnv/modelsegments/{aliquot_barcode}/{aliquot_barcode}.cr.seg",
        "results/cnv/modelsegments/{aliquot_barcode}/{aliquot_barcode}.modelBegin.af.param",
        "results/cnv/modelsegments/{aliquot_barcode}/{aliquot_barcode}.modelBegin.cr.param",
        "results/cnv/modelsegments/{aliquot_barcode}/{aliquot_barcode}.modelFinal.af.param",
        "results/cnv/modelsegments/{aliquot_barcode}/{aliquot_barcode}.modelFinal.cr.param", #"results/cnv/modelsegments/{aliquot_barcode}/{aliquot_barcode}.hets.normal.tsv",
        "results/cnv/modelsegments/{aliquot_barcode}/{aliquot_barcode}.hets.tsv"
    params:
        mem = CLUSTER_META["modelsegments"]["mem"],
        outputdir = "results/cnv/modelsegments/{aliquot_barcode}",
        outputprefix = "{aliquot_barcode}",
        normal_counts_cmd = lambda wildcards: "--normal-allelic-counts results/cnv/hets/{}.allelicCounts.tsv".format(manifest.getMatchedNormal(wildcards.aliquot_barcode)) if manifest.getMatchedNormal(wildcards.aliquot_barcode) is not None else ""
    threads:
        CLUSTER_META["modelsegments"]["ppn"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/cnv/modelsegments/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/cnv/modelsegments/{aliquot_barcode}.txt"
    message:
        "Model segments\n"
        "Aliquot barcode: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g ModelSegments \
            --denoised-copy-ratios {input.denoised} \
            --allelic-counts {input.counts} \
            {params.normal_counts_cmd} \
            --output {params.outputdir} \
            --output-prefix {params.outputprefix} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Call Segments
## Systematic calling of copy-neutral, aplified and deleted segments
## See: https://gatkforums.broadinstitute.org/dsde/discussion/11683/
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule callsegments:
    input:
        "results/cnv/modelsegments/{aliquot_barcode}/{aliquot_barcode}.cr.seg"
    output:
        "results/cnv/callsegments/{aliquot_barcode}.called.seg"
    params:
        mem = CLUSTER_META["callsegments"]["mem"]
    threads:
        CLUSTER_META["callsegments"]["ppn"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/cnv/callsegments/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/cnv/callsegments/{aliquot_barcode}.txt"
    message:
        "Call segments\n"
        "Aliquot barcode: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g CallCopyRatioSegments \
            --input {input} \
            --output {output} \
            > {log} 2>&1"


rule seg2db:
    input:
        seg = lambda wildcards: expand("results/cnv/callsegments/{aliquot_barcode}.called.seg", aliquot_barcode = [aliquot for aliquot in manifest.getSelectedAliquots() if manifest.getBatch(aliquot) is not None])
    output:
        tsv = "results/cnv/callsegments.merged.tsv"
    params:
        mem = CLUSTER_META["seg2db"]["mem"]
    threads:
        CLUSTER_META["seg2db"]["ppn"]
    log:
        "logs/cnv/seg2db/seg2db.log"
    benchmark:
        "benchmarks/cnv/seg2db/seg2db.txt"
    message:
        "Merge seg file and convert to TSV for easy database upload"
    script:
        "../R/snakemake/seg2db.R"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Plot modeled segments
## Use a guassian-kernel binary-segmentation algorithm to group contiguouis copy ratios into segments
## See: https://gatkforums.broadinstitute.org/dsde/discussion/11683/
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule plotmodeledsegments:
    input:
        denoisedCR = "results/cnv/mergedenoisedreadcounts/{aliquot_barcode}.denoisedCR.tsv",
        hets = "results/cnv/modelsegments/{aliquot_barcode}/{aliquot_barcode}.hets.tsv",
        segments = "results/cnv/modelsegments/{aliquot_barcode}/{aliquot_barcode}.modelFinal.seg"
    output:
        "results/cnv/plotmodeledsegments/{aliquot_barcode}.modeled.png"
    params:
        mem = CLUSTER_META["plotmodeledsegments"]["mem"],
        outputdir = "results/cnv/plotmodeledsegments",
        outputprefix = "{aliquot_barcode}"
    threads:
        CLUSTER_META["plotmodeledsegments"]["ppn"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/cnv/plotmodeledsegments/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/cnv/plotmodeledsegments/{aliquot_barcode}.txt"
    message:
        "Plot modelled segments\n"
        "Aliquot barcode: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g PlotModeledSegments \
            --denoised-copy-ratios {input.denoisedCR} \
            --allelic-counts {input.hets} \
            --segments {input.segments} \
            --sequence-dictionary {config[reference_dict]} \
            --minimum-contig-length 46709983 \
            --output {params.outputdir} \
            --output-prefix {params.outputprefix} \
            > {log} 2>&1"

rule combinecnvplots:
    input:
        ms      = "results/cnv/plotmodeledsegments/{aliquot_barcode}.modeled.png",
        cr      = "results/cnv/plotcr/{aliquot_barcode}/{aliquot_barcode}.denoisedLimit4.png"
    output:
        "results/cnv/plots/{aliquot_barcode}.pdf"
    threads:
        CLUSTER_META["combinecnvplots"]["ppn"]
    log:
        "logs/cnv/combinecnvplots/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/cnv/combinecnvplots/{aliquot_barcode}.txt"
    message:
        "Plot GATK CNV results\n"
        "Aliquot: {wildcards.aliquot_barcode}"
    shell:
        "/opt/software/helix/ImageMagick/7.0.7-26/bin/montage {input.cr} {input.ms} -tile 1x2 -geometry +0+0 {output}"

## END ##
