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
            -L {config[cnv][intervals]} \
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
        "Analysis type: {wildcards.analysis_type}"
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
        sample = lambda wildcards: "results/cnv/readcounts/{aliquot_barcode}.counts.hdf5".format(aliquot_barcode = wildcards.aliquot_barcode),
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
            --minimum-contig-length {config[cnv][min_contig_len]} \
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
        denoised = "results/cnv/denoisereadcounts/{aliquot_barcode}.denoisedCR.tsv",
        counts = "results/cnv/alleliccounts/{aliquot_barcode}.allelicCounts.tsv"
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
        normal_counts_cmd = lambda wildcards: "--normal-allelic-counts results/cnv/alleliccounts/{}.allelicCounts.tsv".format(manifest.getMatchedNormal(wildcards.aliquot_barcode)) if manifest.getMatchedNormal(wildcards.aliquot_barcode) is not None else ""
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

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Plot modeled segments
## Use a guassian-kernel binary-segmentation algorithm to group contiguouis copy ratios into segments
## See: https://gatkforums.broadinstitute.org/dsde/discussion/11683/
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule plotmodeledsegments:
    input:
        tumor_denoised = lambda wildcards: "results/cnv/denoisereadcounts/{aliquot_barcode}.denoisedCR.tsv",
        tumor_counts = "results/cnv/modelsegments/{aliquot_barcode}/{aliquot_barcode}.hets.tsv",
        tumor_segments = "results/cnv/modelsegments/{aliquot_barcode}/{aliquot_barcode}.modelFinal.seg"
    output:
        "results/cnv/plotmodeledsegments/{aliquot_barcode}/{aliquot_barcode}.modeled.png"
    params:
        mem = CLUSTER_META["plotmodeledsegments"]["mem"],
        outputdir = "results/cnv/plotmodeledsegments/{aliquot_barcode}",
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
            --denoised-copy-ratios {input.tumor_denoised} \
            --allelic-counts {input.tumor_counts} \
            --segments {input.tumor_segments} \
            --sequence-dictionary {config[reference_dict]} \
            --minimum-contig-length 46709983 \
            --output {params.outputdir} \
            --output-prefix {params.outputprefix} \
            > {log} 2>&1"

rule combinetracks:
    input:
        tumor_called_seg = lambda wildcards: "results/cnv/callsegments/{aliquot_barcode}.called.seg".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode)),
        normal_called_seg = lambda wildcards: "results/cnv/callsegments/{aliquot_barcode}.called.seg".format(aliquot_barcode = manifest.getNormal(wildcards.pair_barcode))
    output:
        germline_tagged_seg = temp("results/cnv/combinetracks/{pair_barcode}.germline_tagged.seg"),
        centromere_tagged_seg = temp("results/cnv/combinetracks/{pair_barcode}.centromere_tagged.seg"),
        final_seg = "results/cnv/combinetracks/{pair_barcode}.final.seg"
    params:
        mem = CLUSTER_META["combinetracks"]["mem"]
    threads:
        CLUSTER_META["combinetracks"]["ppn"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/cnv/combinetracks/{pair_barcode}.log"
    benchmark:
        "benchmarks/cnv/combinetracks/{pair_barcode}.txt"
    message:
        "Combine tumor and normal segmentation\n"
        "Pair ID: {wildcards.pair_barcode}"
    shell:
        "echo \"======= Germline Tagging\";"
        "gatk --java-options -Xmx{params.mem}g TagGermlineEvents \
            --segments {input.tumor_called_seg} \
            --called-matched-normal-seg-file {input.normal_called_seg} \
            -O {output.germline_tagged_seg} \
            -R {config[reference_fasta]} \
            > {log} 2>&1;"
        
        "echo \"======= Centromeres\";"
        "gatk --java-options -Xmx{params.mem}g CombineSegmentBreakpoints \
            --segments {output.germline_tagged_seg} \
            --segments {config[cnv][centromere]} \
            --columns-of-interest NUM_POINTS_COPY_RATIO \
            --columns-of-interest MEAN_LOG2_COPY_RATIO \
            --columns-of-interest CALL \
            --columns-of-interest POSSIBLE_GERMLINE \
            --columns-of-interest type \
            -O {output.centromere_tagged_seg} \
            -R {config[reference_fasta]} \
            >> {log} 2>&1;"
        
        "echo \"======= GISTIC blacklist\";"
        "gatk --java-options -Xmx{params.mem}g CombineSegmentBreakpoints \
            --segments {output.centromere_tagged_seg} \
            --segments {config[cnv][gistic]} \
            --columns-of-interest NUM_POINTS_COPY_RATIO \
            --columns-of-interest MEAN_LOG2_COPY_RATIO \
            --columns-of-interest CALL \
            --columns-of-interest POSSIBLE_GERMLINE \
            --columns-of-interest type \
            --columns-of-interest ID \
            -O {output.final_seg} \
            -R {config[reference_fasta]} \
            >> {log} 2>&1;"

# rule prepare_acs:
#     input:
#         called_seg = "results/cnv/combinetracks/{pair_barcode}.final.seg",
#         modeled_seg = "results/cnv/modelsegments/{pair_barcode}/{pair_barcode}.modelFinal.seg" # "results/cnv/callsegments/{aliquot_barcode}.called.seg".format(aliquot_barcode = manifest.getNormal(wildcards.pair_barcode))
#     output:
#         temp("results/cnv/prepare_acs/{pair_barcode}.prepare_acs.seg")
#     params:
#         mem = CLUSTER_META["prepare_acs"]["mem"],
#     threads:
#         CLUSTER_META["prepare_acs"]["ppn"]
#     conda:
#         "../envs/gatk4.yaml"
#     log:
#         "logs/cnv/prepare_acs/{pair_barcode}.log"
#     benchmark:
#         "benchmarks/cnv/prepare_acs/{pair_barcode}.txt"
#     message:
#         "Prepare ACS conversion by Merging GATK Model Seg and GATK Segment caller file\n"
#         "Pair ID: {wildcards.pair_barcode}"
#     shell:
#         "echo \"======= Merging GATK Model Seg and GATK Segment caller file\""
#         "gatk --java-options -Xmx{params.mem}g \
#             CombineSegmentBreakpoints \
#             --segments {input.called_seg} \
#             --segments {input.modeled_seg} \
#             --columns-of-interest NUM_POINTS_COPY_RATIO \
#             --columns-of-interest NUM_POINTS_ALLELE_FRACTION \
#             --columns-of-interest LOG2_COPY_RATIO_POSTERIOR_10 \
#             --columns-of-interest LOG2_COPY_RATIO_POSTERIOR_50 \
#             --columns-of-interest LOG2_COPY_RATIO_POSTERIOR_90 \
#             --columns-of-interest MINOR_ALLELE_FRACTION_POSTERIOR_10 \
#             --columns-of-interest MINOR_ALLELE_FRACTION_POSTERIOR_50 \
#             --columns-of-interest MINOR_ALLELE_FRACTION_POSTERIOR_90 \
#             --columns-of-interest CALL \
#             --columns-of-interest NUM_POINTS_COPY_RATIO \
#             --columns-of-interest MEAN_LOG2_COPY_RATIO \
#             --columns-of-interest POSSIBLE_GERMLINE \
#             --columns-of-interest type \
#             --columns-of-interest ID \
#             -O {output} \
#             -R {config[reference_fasta]} \
#             > {log} 2>&1"

# rule filter_tagged:
#     input:
#         "results/cnv/prepare_acs/{pair_barcode}.prepare_acs.seg"
#     output:
#         temp("results/cnv/filter_tagged/{pair_barcode}.pruned.seg")
#     params:
#         mem = CLUSTER_META["filter_tagged"]["mem"],
#     threads:
#         CLUSTER_META["filter_tagged"]["ppn"]
#     conda:
#         "../envs/gatk4.yaml"
#     log:
#         "logs/cnv/filter_tagged/{pair_barcode}.log"
#     benchmark:
#         "benchmarks/cnv/filter_tagged/{pair_barcode}.txt"
#     message:
#         "Prune tagged variants\n"
#         "Pair ID: {wildcards.pair_barcode}"
#     run:
#         import pandas
#         import os.path
#         tumor_tagged_df = pandas.read_csv(input, delimiter="\t", comment="@")
#         tumor_tagged_pruned_df = tumor_tagged_df[(tumor_tagged_df["POSSIBLE_GERMLINE"] == "0") & (tumor_tagged_df["type"] != "centromere") & (tumor_tagged_df["ID"].isna())]
#         output_filename = output
#         print(output_filename)
#         tumor_tagged_pruned_df.to_csv(output_filename, sep="\t", index=False)

# rule merge_annotation:
#     input:
#         "results/cnv/filter_tagged/{pair_barcode}.pruned.seg"
#     output:
#         ...
#     params:
#         mem = CLUSTER_META["merge_annotation"]["mem"],
#     threads:
#         CLUSTER_META["merge_annotation"]["ppn"]
#     conda:
#         "../envs/gatk4.yaml"
#     log:
#         "logs/cnv/merge_annotation/{pair_barcode}.log"
#     benchmark:
#         "benchmarks/cnv/merge_annotation/{pair_barcode}.txt"
#     message:
#         "Prune tagged variants\n"
#         "Pair ID: {wildcards.pair_barcode}"

## END ##