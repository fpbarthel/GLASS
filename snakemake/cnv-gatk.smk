## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Collect read counts
## See: https://software.broadinstitute.org/gatk/documentation/article?id=11682
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule collectreadcounts:
    input:
        "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam"
    output:
        "results/cnv/readcounts/{aliquot_barcode}.counts.hdf5"
    params:
        mem = CLUSTER_META["collectreadcounts"]["mem"]
    threads:
        CLUSTER_META["collectreadcounts"]["ppn"]
    log:
        "logs/cnv/readcounts/{aliquot_barcode}.log"
    conda:
        "../envs/gatk4.yaml"
    benchmark:
        "benchmarks/cnv/readcounts/{aliquot_barcode}.txt"
    message:
        "Collecting read counts\n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g CollectReadCounts \
            -I {input} \
            -L {config[cnv_intervals]} \
            --interval-merging-rule OVERLAPPING_ONLY \
            -O {output} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Create CNV panel of normals
## See: https://software.broadinstitute.org/gatk/documentation/article?id=11682
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule createcnvpon:
    input:
        lambda wildcards: expand("results/cnv/readcounts/{aliquot_barcode}.counts.hdf5", aliquot_barcode = manifest.getPONAliquots()) #wildcards.analysis_type))
    output:
        "results/cnv/createcnvpon/{analysis_type}.pon.hdf5"
    params:
        mem = CLUSTER_META["createcnvpon"]["mem"]
    threads:
        CLUSTER_META["createcnvpon"]["ppn"]
    log:
        "logs/cnv/createcnvpon/{analysis_type}.log"
    benchmark:
        "benchmarks/cnv/createcnvpon/{analysis_type}.txt"
    message:
        "Creating CNV panel of normals\n"
        "Analysis type: {wildcards.analysis_type}\n"
        "Batch: {wildcards.analysis_type}"
    run:
        vcfs = " ".join(["-I " + s for s in input])
        shell("gatk --java-options -Xmx{params.mem}g CreateReadCountPanelOfNormals \
                {vcfs} \
                --minimum-interval-median-percentile 5.0 \
                -O {output} \
                > {log} 2>&1")

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Denoise read counts
## See: https://software.broadinstitute.org/gatk/documentation/article?id=11682
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule denoisereadcounts:
    input:
        sample = lambda wildcards: "results/cnv/readcounts/{analysis_type}/{aliquot_barcode}.counts.hdf5".format(analysis_type = "WGS", aliquot_barcode = wildcards.aliquot_barcode),
        pon =  lambda wildcards: "results/cnv/createcnvpon/{analysis_type}.pon.hdf5".format(analysis_type = "WGS")
    output:
        standardized = "results/cnv/denoisereadcounts/{aliquot_barcode}.standardizedCR.tsv",
        denoised = "results/cnv/denoisereadcounts/{aliquot_barcode}.denoisedCR.tsv"
    params:
        mem = CLUSTER_META["denoisereadcounts"]["mem"]
    threads:
        CLUSTER_META["denoisereadcounts"]["ppn"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/cnv/denoisereadcounts/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/cnv/denoisereadcounts/{aliquot_barcode}.txt"
    message:
        "Denoising and standardizing read counts\n"
        "Aliquot: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g DenoiseReadCounts \
            -I {input.sample} \
            --count-panel-of-normals {input.pon} \
            --standardized-copy-ratios {output.standardized} \
            --denoised-copy-ratios {output.denoised} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Plot denoised and standardized copy ratios
## See: https://software.broadinstitute.org/gatk/documentation/article?id=11682
## NB need to optimize minimum contig length parameter for b37
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule plotcr:
    input:
        standardized = "results/cnv/denoisereadcounts/{aliquot_barcode}.standardizedCR.tsv",
        denoised = "results/cnv/denoisereadcounts/{aliquot_barcode}.denoisedCR.tsv"
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
            --minimum-contig-length 46709983 \
            --output {params.outputdir} \
            --output-prefix {params.outputprefix} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Collect Allelic Counts
## Counts reference and alternative alleles at common germline variant sites
## See: https://gatkforums.broadinstitute.org/dsde/discussion/11683/
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule collectalleliccounts:
    input:
        "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam"
    output:
        "results/cnv/alleliccounts/{aliquot_barcode}.allelicCounts.tsv"
    params:
        mem = CLUSTER_META["collectalleliccounts"]["mem"]
    threads:
        CLUSTER_META["collectalleliccounts"]["ppn"]
    log:
        "logs/cnv/alleliccounts/{aliquot_barcode}.log"
    conda:
        "../envs/gatk4.yaml"
    benchmark:
        "benchmarks/cnv/alleliccounts/{aliquot_barcode}.txt"
    message:
        "Collect allelic counts\n"
        "Aliquot: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g CollectAllelicCounts \
            -I {input} \
            -L {config[cnv_common_snp]} \
            -R {config[reference_fasta]} \
            -O {output} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Model segments
## Use a guassian-kernel binary-segmentation algorithm to group contiguouis copy ratios into segments
## See: https://gatkforums.broadinstitute.org/dsde/discussion/11683/
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule modelsegments:
    input:
        tumor_denoised = lambda wildcards: "results/cnv/denoisereadcounts/{aliquot_barcode}.denoisedCR.tsv".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode)),
        tumor_counts = lambda wildcards: "results/cnv/alleliccounts/{aliquot_barcode}.allelicCounts.tsv".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode)),
        normal_counts = lambda wildcards: "results/cnv/alleliccounts/{aliquot_barcode}.allelicCounts.tsv".format(aliquot_barcode = manifest.getNormal(wildcards.pair_barcode))
    output:
        "results/cnv/modelsegments/{pair_barcode}/{pair_barcode}.modelBegin.seg",
        "results/cnv/modelsegments/{pair_barcode}/{pair_barcode}.modelFinal.seg",
        "results/cnv/modelsegments/{pair_barcode}/{pair_barcode}.cr.seg",
        "results/cnv/modelsegments/{pair_barcode}/{pair_barcode}.modelBegin.af.param",
        "results/cnv/modelsegments/{pair_barcode}/{pair_barcode}.modelBegin.cr.param",
        "results/cnv/modelsegments/{pair_barcode}/{pair_barcode}.modelFinal.af.param",
        "results/cnv/modelsegments/{pair_barcode}/{pair_barcode}.modelFinal.cr.param",
        "results/cnv/modelsegments/{pair_barcode}/{pair_barcode}.hets.normal.tsv",
        "results/cnv/modelsegments/{pair_barcode}/{pair_barcode}.hets.tsv"
    params:
        mem = CLUSTER_META["modelsegments"]["mem"],
        outputdir = "results/cnv/modelsegments/{pair_barcode}",
        outputprefix = "{pair_barcode}"
    threads:
        CLUSTER_META["modelsegments"]["ppn"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/cnv/modelsegments/{pair_barcode}.log"
    benchmark:
        "benchmarks/cnv/modelsegments/{pair_barcode}.txt"
    message:
        "Model segments\n"
        "Pair ID: {wildcards.pair_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g ModelSegments \
            --denoised-copy-ratios {input.tumor_denoised} \
            --allelic-counts {input.tumor_counts} \
            --normal-allelic-counts {input.normal_counts} \
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
        "results/cnv/modelsegments/{pair_barcode}/{pair_barcode}.cr.seg"
    output:
        "results/cnv/callsegments/{pair_barcode}.called.seg"
    params:
        mem = CLUSTER_META["callsegments"]["mem"]
    threads:
        CLUSTER_META["callsegments"]["ppn"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/cnv/callsegments/{pair_barcode}.log"
    benchmark:
        "benchmarks/cnv/callsegments/{pair_barcode}.txt"
    message:
        "Call segments\n"
        "Pair ID: {wildcards.pair_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g CallCopyRatioSegments \
            --input {input} \
            --output {output} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Plot modeled segments
## Use a guassian-kernel binary-segmentation algorithm to group contiguouis copy ratios into segments
## See: https://gatkforums.broadinstitute.org/dsde/discussion/11683/
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule plotmodeledsegments:
    input:
        tumor_denoised = lambda wildcards: "results/cnv/denoisereadcounts/{aliquot_barcode}.denoisedCR.tsv".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode)),
        tumor_counts = "results/cnv/modelsegments/{pair_barcode}/{pair_barcode}.hets.tsv",
        tumor_segments = "results/cnv/modelsegments/{pair_barcode}/{pair_barcode}.modelFinal.seg"
    output:
        "results/cnv/plotmodeledsegments/{pair_barcode}/{pair_barcode}.modeled.png"
    params:
        mem = CLUSTER_META["plotmodeledsegments"]["mem"],
        outputdir = "results/cnv/plotmodeledsegments/{pair_barcode}",
        outputprefix = "{pair_barcode}"
    threads:
        CLUSTER_META["plotmodeledsegments"]["ppn"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/cnv/plotmodeledsegments/{pair_barcode}.log"
    benchmark:
        "benchmarks/cnv/plotmodeledsegments/{pair_barcode}.txt"
    message:
        "Plot modelled segments\n"
        "Pair ID: {wildcards.pair_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g PlotModeledSegments \
            --denoised-copy-ratios {input.tumor_denoised} \
            --allelic-counts {input.tumor_counts} \
            --segments {input.tumor_segments} \
            --sequence-dictionary {config[reference_dict]} \
            --minimum-contig-length 46709983 \
            --output {params.outputdir} \
            --output-prefix {params.outputprefix} \
            > {log} 2>&1"

## END ##