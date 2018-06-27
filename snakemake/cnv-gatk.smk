## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Collect read counts
## See: https://software.broadinstitute.org/gatk/documentation/article?id=11682
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule collectreadcounts:
    input:
        "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam"
    output:
        "results/cnv/readcounts/{batch}/{aliquot_id}.counts.hdf5"
    params:
        mem = CLUSTER_META["collectreadcounts"]["mem"]
    threads:
        CLUSTER_META["collectreadcounts"]["ppn"]
    log:
        "logs/cnv/readcounts/{batch}/{aliquot_id}.log"
    benchmark:
        "benchmarks/cnv/readcounts/{batch}/{aliquot_id}.txt"
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
        lambda wildcards: expand("results/cnv/readcounts/{batch}/{aliquot_id}.counts.hdf5", batch=wildcards.batch, aliquot_id=BATCH_TO_NORMAL[wildcards.batch]) 
    output:
        "results/cnv/createcnvpon/{batch}.pon.hdf5"
    params:
        mem = CLUSTER_META["createcnvpon"]["mem"]
    threads:
        CLUSTER_META["createcnvpon"]["ppn"]
    log:
        "logs/cnv/createcnvpon/{batch}.log"
    benchmark:
        "benchmarks/cnv/createcnvpon/{batch}.txt"
    message:
        "Creating CNV panel of normals\n"
        "Batch: {wildcards.batch}"
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
        sample = lambda wildcards: "results/cnv/readcounts/{batch}/{aliquot_id}.counts.hdf5".format(batch = CASES_DICT[SAMPLES_DICT[ALIQUOTS_DICT[wildcards.aliquot_id]["sample_id"]]["case_id"]]["project_id"], aliquot_id=wildcards.aliquot_id),
        pon =  lambda wildcards: "results/cnv/createcnvpon/{batch}.pon.hdf5".format(batch = CASES_DICT[SAMPLES_DICT[ALIQUOTS_DICT[wildcards.aliquot_id]["sample_id"]]["case_id"]]["project_id"])
    output:
        standardized = "results/cnv/denoisereadcounts/{aliquot_id}.standardizedCR.tsv",
        denoised = "results/cnv/denoisereadcounts/{aliquot_id}.denoisedCR.tsv"
    params:
        mem = CLUSTER_META["denoisereadcounts"]["mem"]
    threads:
        CLUSTER_META["denoisereadcounts"]["ppn"]
    log:
        "logs/cnv/denoisereadcounts/{aliquot_id}.log"
    benchmark:
        "benchmarks/cnv/denoisereadcounts/{aliquot_id}.txt"
    message:
        "Denoising and standardizing read counts\n"
        "Aliquot: {wildcards.aliquot_id}"
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
        standardized = "results/cnv/denoisereadcounts/{aliquot_id}.standardizedCR.tsv",
        denoised = "results/cnv/denoisereadcounts/{aliquot_id}.denoisedCR.tsv"
    output:
        "results/cnv/plotcr/{aliquot_id}/{aliquot_id}.denoised.png",
        "results/cnv/plotcr/{aliquot_id}/{aliquot_id}.denoisedLimit4.png",
        "results/cnv/plotcr/{aliquot_id}/{aliquot_id}.standardizedMAD.txt",
        "results/cnv/plotcr/{aliquot_id}/{aliquot_id}.denoisedMAD.txt",
        "results/cnv/plotcr/{aliquot_id}/{aliquot_id}.deltaMAD.txt",
        "results/cnv/plotcr/{aliquot_id}/{aliquot_id}.scaledDeltaMAD.txt"
    params:
        mem = CLUSTER_META["plotcr"]["mem"],
        outputdir = "results/cnv/plotcr/{aliquot_id}",
        outputprefix = "{aliquot_id}"
    threads:
        CLUSTER_META["plotcr"]["ppn"]
    log:
        "logs/cnv/plotcr/{aliquot_id}.log"
    benchmark:
        "benchmarks/cnv/plotcr/{aliquot_id}.txt"
    message:
        "Plot denoised and standardized read counts\n"
        "Aliquot: {wildcards.aliquot_id}"
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
        "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam"
    output:
        "results/cnv/alleliccounts/{aliquot_id}.allelicCounts.tsv"
    params:
        mem = CLUSTER_META["collectalleliccounts"]["mem"]
    threads:
        CLUSTER_META["collectalleliccounts"]["ppn"]
    log:
        "logs/cnv/alleliccounts/{aliquot_id}.log"
    benchmark:
        "benchmarks/cnv/alleliccounts/{aliquot_id}.txt"
    message:
        "Collect allelic counts\n"
        "Aliquot: {wildcards.aliquot_id}"
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
        tumor_denoised = lambda wildcards: "results/cnv/denoisereadcounts/{aliquot_id}.denoisedCR.tsv".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["tumor_aliquot_id"]),
        tumor_counts = lambda wildcards: "results/cnv/alleliccounts/{aliquot_id}.allelicCounts.tsv".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["tumor_aliquot_id"]),
        normal_counts = lambda wildcards: "results/cnv/alleliccounts/{aliquot_id}.allelicCounts.tsv".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["normal_aliquot_id"])
    output:
        "results/cnv/modelsegments/{pair_id}/{pair_id}.modelBegin.seg",
        "results/cnv/modelsegments/{pair_id}/{pair_id}.modelFinal.seg",
        "results/cnv/modelsegments/{pair_id}/{pair_id}.cr.seg",
        "results/cnv/modelsegments/{pair_id}/{pair_id}.modelBegin.af.param",
        "results/cnv/modelsegments/{pair_id}/{pair_id}.modelBegin.cr.param",
        "results/cnv/modelsegments/{pair_id}/{pair_id}.modelFinal.af.param",
        "results/cnv/modelsegments/{pair_id}/{pair_id}.modelFinal.cr.param",
        "results/cnv/modelsegments/{pair_id}/{pair_id}.hets.normal.tsv",
        "results/cnv/modelsegments/{pair_id}/{pair_id}.hets.tsv"
    params:
        mem = CLUSTER_META["modelsegments"]["mem"],
        outputdir = "results/cnv/modelsegments/{pair_id}",
        outputprefix = "{pair_id}"
    threads:
        CLUSTER_META["modelsegments"]["ppn"]
    log:
        "logs/cnv/modelsegments/{pair_id}.log"
    benchmark:
        "benchmarks/cnv/modelsegments/{pair_id}.txt"
    message:
        "Model segments\n"
        "Pair ID: {wildcards.pair_id}"
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
        "results/cnv/modelsegments/{pair_id}/{pair_id}.cr.seg"
    output:
        "results/cnv/callsegments/{pair_id}.called.seg"
    params:
        mem = CLUSTER_META["callsegments"]["mem"]
    threads:
        CLUSTER_META["callsegments"]["ppn"]
    log:
        "logs/cnv/callsegments/{pair_id}.log"
    benchmark:
        "benchmarks/cnv/callsegments/{pair_id}.txt"
    message:
        "Call segments\n"
        "Pair ID: {wildcards.pair_id}"
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
        tumor_denoised = lambda wildcards: "results/cnv/denoisereadcounts/{aliquot_id}.denoisedCR.tsv".format(aliquot_id=PAIRS_DICT[wildcards.pair_id]["tumor_aliquot_id"]),
        tumor_counts = "results/cnv/modelsegments/{pair_id}/{pair_id}.hets.tsv",
        tumor_segments = "results/cnv/modelsegments/{pair_id}/{pair_id}.modelFinal.seg"
    output:
        "results/cnv/plotmodeledsegments/{pair_id}/{pair_id}.modeled.png"
    params:
        mem = CLUSTER_META["plotmodeledsegments"]["mem"],
        outputdir = "results/cnv/plotmodeledsegments/{pair_id}",
        outputprefix = "{pair_id}"
    threads:
        CLUSTER_META["plotmodeledsegments"]["ppn"]
    log:
        "logs/cnv/plotmodeledsegments/{pair_id}.log"
    benchmark:
        "benchmarks/cnv/plotmodeledsegments/{pair_id}.txt"
    message:
        "Plot modelled segments\n"
        "Pair ID: {wildcards.pair_id}"
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