## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Collect read counts
## See: https://software.broadinstitute.org/gatk/documentation/article?id=11682
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule collectreadcounts:
    input:
        "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam"
    output:
        "results/collectreadcounts/{batch}/{aliquot_id}.counts.hdf5"
    params:
        mem = CLUSTER_META["collectreadcounts"]["mem"]
    threads:
        CLUSTER_META["collectreadcounts"]["ppn"]
    log:
        "logs/collectreadcounts/{batch}/{aliquot_id}.log"
    benchmark:
        "benchmarks/collectreadcounts/{batch}/{aliquot_id}.txt"
    message:
        "Collecting read counts\n"
        "Batch: {wildcards.batch}\n"
        "Sample: {wildcards.aliquot_id}"
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
        lambda wildcards: expand("results/collectreadcounts/{batch}/{aliquot_id}.pon.vcf", batch=wildcards.batch, aliquot_id=BATCH_TO_NORMAL[wildcards.batch]) 
    output:
        "results/createcnvpon/{batch}.pon.hdf5"
    params:
        mem = CLUSTER_META["createcnvpon"]["mem"]
    threads:
        CLUSTER_META["createcnvpon"]["ppn"]
    log:
        "logs/createcnvpon/{batch}.log"
    benchmark:
        "benchmarks/createcnvpon/{batch}.txt"
    message:
        "Creating CNV panel of normals\n"
        "Batch: {wildcards.batch}"
    run:
        vcfs = " ".join(["-I " + s for s in input])
        shell("gatk --java-options -Xmx{params.mem}g CreateReadCountPanelOfNormals \
                {vcfs} \
                --minimum-interval-median-percentile 5.0 \
                -O {output} \
                > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Denoise read counts
## See: https://software.broadinstitute.org/gatk/documentation/article?id=11682
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule denoisereadcounts:
    input:
        sample = "results/collectreadcounts/{batch}/{aliquot_id}.counts.hdf5",
        pon = "results/createcnvpon/{batch}.pon.hdf5"
    output:
        standardized = "results/denoisereadcounts/{aliquot_id}.standardizedCR.tsv",
        denoised = "results/denoisereadcounts/{aliquot_id}.denoisedCR.tsv"
    params:
        mem = CLUSTER_META["denoisereadcounts"]["mem"]
    threads:
        CLUSTER_META["denoisereadcounts"]["ppn"]
    log:
        "logs/denoisereadcounts/{aliquot_id}.log"
    benchmark:
        "benchmarks/denoisereadcounts/{aliquot_id}.txt"
    message:
        "Denoising and standardizing read counts\n"
        "Aliquot: {wildcards.aliquot_id}"
    shell:
        "gatk --java-options -Xmx{params.mem}g DenoiseReadCounts \
            -I {input.sample} \
            --count-panel-of-normals {input.pon} \
            --interval-merging-rule OVERLAPPING_ONLY \
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
        standardized = "results/denoisereadcounts/{aliquot_id}.standardizedCR.tsv",
        denoised = "results/denoisereadcounts/{aliquot_id}.denoisedCR.tsv"
    output:
        "results/plotcr/{aliquot_id}/{aliquot_id}.denoised.png",
        "results/plotcr/{aliquot_id}/{aliquot_id}.denoisedLimit4.png",
        "results/plotcr/{aliquot_id}/{aliquot_id}.standardizedMAD.txt",
        "results/plotcr/{aliquot_id}/{aliquot_id}.denoisedMAD.txt",
        "results/plotcr/{aliquot_id}/{aliquot_id}.deltaMAD.txt",
        "results/plotcr/{aliquot_id}/{aliquot_id}.scaledDeltaMAD.txt"
    params:
        mem = CLUSTER_META["plotcr"]["mem"],
        outputdir = "results/plotcr/{aliquot_id}",
        outputprefix = "{aliquot_id}"
    threads:
        CLUSTER_META["plotcr"]["ppn"]
    log:
        "logs/plotcr/{aliquot_id}.log"
    benchmark:
        "benchmarks/plotcr/{aliquot_id}.txt"
    message:
        "Plot denoised and standardized read counts\n"
        "Aliquot: {wildcards.aliquot_id}"
    shell:
        "gatk --java-options -Xmx{params.mem}g PlotDenoisedCopyRatios \
            --standardized-copy-ratios {input.standardized} \
            --denoised-copy-ratios {input.denoised} \
            --sequence-dictionary {config[reference_dict]} \
            --minimum-contig-length 46709983 \
            --standardized-copy-ratios {output.standardized} \
            --denoised-copy-ratios {output.denoised} \
            > {log} 2>&1"

## END ##