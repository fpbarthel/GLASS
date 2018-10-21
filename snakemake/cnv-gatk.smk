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
        mem = CLUSTER_META["collectreadcounts"]["mem"],
        intervals = lambda wildcards: config["cnv"]["intervals_wes"] if manifest.isExome(wildcards.aliquot_barcode) else config["cnv"]["intervals_wgs"]
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
            -L {params.intervals} \
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

rule prepare_acs:
    input:
        called_seg = "results/cnv/combinetracks/{pair_barcode}.final.seg",
        modeled_seg = lambda wildcards: "results/cnv/modelsegments/{aliquot_barcode}/{aliquot_barcode}.modelFinal.seg".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode)) # "results/cnv/modelsegments/{pair_barcode}/{pair_barcode}.modelFinal.seg" # 
    output:
        temp("results/cnv/prepare_acs/{pair_barcode}.prepare_acs.seg")
    params:
        mem = CLUSTER_META["prepare_acs"]["mem"]
    threads:
        CLUSTER_META["prepare_acs"]["ppn"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/cnv/prepare_acs/{pair_barcode}.log"
    benchmark:
        "benchmarks/cnv/prepare_acs/{pair_barcode}.txt"
    message:
        "Prepare ACS conversion by Merging GATK Model Seg and GATK Segment caller file\n"
        "Pair ID: {wildcards.pair_barcode}"
    shell:
        "echo \"======= Merging GATK Model Seg and GATK Segment caller file\";"
        "gatk --java-options -Xmx{params.mem}g \
            CombineSegmentBreakpoints \
            --segments {input.called_seg} \
            --segments {input.modeled_seg} \
            --columns-of-interest NUM_POINTS_COPY_RATIO \
            --columns-of-interest NUM_POINTS_ALLELE_FRACTION \
            --columns-of-interest LOG2_COPY_RATIO_POSTERIOR_10 \
            --columns-of-interest LOG2_COPY_RATIO_POSTERIOR_50 \
            --columns-of-interest LOG2_COPY_RATIO_POSTERIOR_90 \
            --columns-of-interest MINOR_ALLELE_FRACTION_POSTERIOR_10 \
            --columns-of-interest MINOR_ALLELE_FRACTION_POSTERIOR_50 \
            --columns-of-interest MINOR_ALLELE_FRACTION_POSTERIOR_90 \
            --columns-of-interest CALL \
            --columns-of-interest NUM_POINTS_COPY_RATIO \
            --columns-of-interest MEAN_LOG2_COPY_RATIO \
            --columns-of-interest POSSIBLE_GERMLINE \
            --columns-of-interest type \
            --columns-of-interest ID \
            -O {output} \
            -R {config[reference_fasta]} \
            > {log} 2>&1"

rule filter_tagged:
    input:
        "results/cnv/prepare_acs/{pair_barcode}.prepare_acs.seg"
    output:
        temp("results/cnv/filter_tagged/{pair_barcode}.pruned.seg")
    params:
        mem = CLUSTER_META["filter_tagged"]["mem"]
    threads:
        CLUSTER_META["filter_tagged"]["ppn"]
    #conda:
    #    "../envs/gatk4.yaml"
    log:
        "logs/cnv/filter_tagged/{pair_barcode}.log"
    benchmark:
        "benchmarks/cnv/filter_tagged/{pair_barcode}.txt"
    message:
        "Prune tagged variants\n"
        "Pair ID: {wildcards.pair_barcode}"
    run:
        import pandas
        import os.path
        tumor_tagged_df = pandas.read_csv(input[0], delimiter="\t", comment="@")
        tumor_tagged_pruned_df = tumor_tagged_df[(tumor_tagged_df["POSSIBLE_GERMLINE"] == "0") & (tumor_tagged_df["type"] != "centromere") & (tumor_tagged_df["ID"].isnull())]
        output_filename = output[0]
        print(output_filename)
        tumor_tagged_pruned_df.to_csv(output_filename, sep="\t", index=False)

rule merge_annotation:
    input:
        "results/cnv/filter_tagged/{pair_barcode}.pruned.seg"
    output:
        "results/cnv/merge_annotation/{pair_barcode}.merged.seg"
    params:
        mem = CLUSTER_META["merge_annotation"]["mem"]
    threads:
        CLUSTER_META["merge_annotation"]["ppn"]
    #conda: ## COMMENTED OUT BECAUSE USES GATK 4.0.10.1 which is not on conda
    #    "../envs/gatk4.yaml"
    log:
        "logs/cnv/merge_annotation/{pair_barcode}.log"
    benchmark:
        "benchmarks/cnv/merge_annotation/{pair_barcode}.txt"
    message:
        "Merge adjacent genomic regions when annotations match\n"
        "Pair ID: {wildcards.pair_barcode}"
    shell:
        "echo \"======= Merging\";"
        "gatk --java-options -Xmx{params.mem}g MergeAnnotatedRegionsByAnnotation \
            --segments {input} \
            --annotations-to-match MEAN_LOG2_COPY_RATIO \
            --annotations-to-match LOG2_COPY_RATIO_POSTERIOR_10 \
            --annotations-to-match LOG2_COPY_RATIO_POSTERIOR_50 \
            --annotations-to-match LOG2_COPY_RATIO_POSTERIOR_90 \
            --annotations-to-match MINOR_ALLELE_FRACTION_POSTERIOR_10 \
            --annotations-to-match MINOR_ALLELE_FRACTION_POSTERIOR_50 \
            --annotations-to-match MINOR_ALLELE_FRACTION_POSTERIOR_90 \
            -O {output} \
            -R {config[reference_fasta]} \
            > {log} 2>&1"

rule acs_convert:
    input:
        af_param = lambda wildcards: "results/cnv/modelsegments/{aliquot_barcode}/{aliquot_barcode}.modelFinal.af.param".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode)),
        merged_seg = "results/cnv/merge_annotation/{pair_barcode}.merged.seg"
    output:
        seg = "results/cnv/acs_convert/{pair_barcode}.acs.seg",
        skew = "results/cnv/acs_convert/{pair_barcode}.skew"
    params:
        mem = CLUSTER_META["acs_convert"]["mem"]
    threads:
        CLUSTER_META["acs_convert"]["ppn"]
    #conda:
    #    "../envs/gatk4.yaml"
    log:
        "logs/cnv/acs_convert/{pair_barcode}.log"
    benchmark:
        "benchmarks/cnv/acs_convert/{pair_barcode}.txt"
    message:
        "Produces a seg file of identical format to AllelicCapSeg\n"
        "Pair ID: {wildcards.pair_barcode}"
    run:
        import sys
        import re
        import pandas as pd
        import numpy as np
        from collections import defaultdict
        import scipy
        from scipy import special as sp
        import os.path

        model_segments_seg_input_file = input["merged_seg"]
        model_segments_af_param_input_file = input["af_param"]
        alleliccapseg_seg_output_file = output["seg"]
        alleliccapseg_skew_output_file = output["skew"]

        HAM_FIST_THRESHOLD = config["cnv"]["maf90_threshold"]

        # regular expression for matching sample name from header comment line
        sample_name_header_regexp = "^@RG.*SM:(.*)[\t]*.*$"

        #define AllelicCapSeg columns
        alleliccapseg_seg_columns = [
            'Chromosome',
            'Start.bp',
            'End.bp',
            'n_probes',
            'length',
            'n_hets',
            'f',
            'tau',
            'sigma.tau',
            'mu.minor',
            'sigma.minor',
            'mu.major',
            'sigma.major',
            'SegLabelCNLOH']

        def read_sample_name(input_file, max_scan_lines=10000):
            with open(input_file, 'r') as f:
                for _ in range(max_scan_lines):
                    line = f.readline()
                    match = re.search(sample_name_header_regexp, line, re.M)
                    if match is None:
                        continue
                    groups = match.groups()
                    return groups[0]
            raise Exception("Sample name could not be found in \"{0}\"".format(input_file))

        #read GATK ModelSegments files and perform some basic checks
        model_segments_seg_pd = pd.read_csv(model_segments_seg_input_file,
                                            sep='\t', comment='@', na_values='NA')
        model_segments_af_param_pd = pd.read_csv(model_segments_af_param_input_file, sep='\t', comment='@')

        def simple_determine_allelic_fraction(model_segments_seg_pd):
            result = model_segments_seg_pd['MINOR_ALLELE_FRACTION_POSTERIOR_50']
            result[model_segments_seg_pd['MINOR_ALLELE_FRACTION_POSTERIOR_90'] > HAM_FIST_THRESHOLD] = 0.5
            return result

        def convert_model_segments_to_alleliccapseg(model_segments_seg_pd,
                                                    model_segments_af_param_pd):
            alleliccapseg_seg_pd = pd.DataFrame(columns=alleliccapseg_seg_columns)

            #The following conversions are trivial.
            alleliccapseg_seg_pd['Chromosome'] = model_segments_seg_pd['CONTIG']
            alleliccapseg_seg_pd['Start.bp'] = model_segments_seg_pd['START']
            alleliccapseg_seg_pd['End.bp'] = model_segments_seg_pd['END']
            alleliccapseg_seg_pd['n_probes'] = model_segments_seg_pd['NUM_POINTS_COPY_RATIO_1']
            alleliccapseg_seg_pd['length'] = alleliccapseg_seg_pd['End.bp'] - alleliccapseg_seg_pd['Start.bp']
            alleliccapseg_seg_pd['n_hets'] = model_segments_seg_pd['NUM_POINTS_ALLELE_FRACTION']

            #ModelSegments estimates posterior credible intervals, while AllelicCapSeg performs maximum a posteriori (MAP) estimation.
            #The copy-ratio and allele-fraction models fit by both also differ.
            # We will attempt a rough translation of the model fits here.

            alleliccapseg_seg_pd['f'] = simple_determine_allelic_fraction(model_segments_seg_pd)

            alleliccapseg_seg_pd['tau'] = 2. * 2**model_segments_seg_pd['LOG2_COPY_RATIO_POSTERIOR_50']
            alleliccapseg_seg_pd['sigma.tau'] = 2**model_segments_seg_pd['LOG2_COPY_RATIO_POSTERIOR_90'] - 2**model_segments_seg_pd['LOG2_COPY_RATIO_POSTERIOR_10']
            sigma_f = (model_segments_seg_pd['MINOR_ALLELE_FRACTION_POSTERIOR_90'].values - model_segments_seg_pd['MINOR_ALLELE_FRACTION_POSTERIOR_10'].values) / 2.
            sigma_mu = np.sqrt(sigma_f**2 + alleliccapseg_seg_pd['sigma.tau']**2) #we propagate errors in the products f * tau and (1 - f) * tau in the usual way
            alleliccapseg_seg_pd['mu.minor'] = alleliccapseg_seg_pd['f'] * alleliccapseg_seg_pd['tau']
            alleliccapseg_seg_pd['sigma.minor'] = sigma_mu
            alleliccapseg_seg_pd['mu.major'] = (1. - alleliccapseg_seg_pd['f']) * alleliccapseg_seg_pd['tau']
            alleliccapseg_seg_pd['sigma.major'] = sigma_mu

            #For whatever reason, AllelicCapSeg attempts to call CNLOH.  Documentation is spotty, but it seems like it attempts
            # to distinguish between three states ("0 is flanked on both sides, 1 is one side, 2 is no cn.loh").
            # Let's just set everything to 2 for now.
            # Hopefully, ABSOLUTE is robust to this ...
            alleliccapseg_seg_pd['SegLabelCNLOH'] = 2

            #One important caveat: for segments with less than 10 hets, AllelicCapSeg also tries to call whether a segment is "split" or not.
            #  This script will attempt to call "split" on all segments.
            # ACS performs a simple hypothesis test on the alternate-allele fractions to see if
            # a unimodal distribution peaked at 0.5 is supported over a bimodal distribution peaked at f and 1 - f.
            # If the former is supported, then AllelicCapSeg ignores the MAP estimate of f and simply sets it to be 0.5.
            # ABSOLUTE may actually be rather sensitive to this.  Again, let's ignore for now, and we can later port this
            # statistical test if necessary.

            #Finally, I believe that ABSOLUTE requires the value of the "skew" parameter from the AllelicCapSeg
            #allele-fraction model.  This parameter is supposed to allow the model to account for reference bias,
            #  but the model likelihood that AllelicCapSeg uses is not valid over the entire range of the skew parameter.
            # We corrected this during the development of AllelicCNV and retain the same corrected model in ModelSegments.
            # We will try to transform the relevant parameter in the corrected model back to a "skew",
            # but this operation is ill defined.  Luckily, for WGS, the reference bias is typically negligible.
            model_segments_reference_bias = model_segments_af_param_pd[
                model_segments_af_param_pd['PARAMETER_NAME'] == 'MEAN_BIAS']['POSTERIOR_50']
            alleliccapseg_skew = 2. / (1. + model_segments_reference_bias)

            return alleliccapseg_seg_pd, alleliccapseg_skew


        #do the conversion
        alleliccapseg_seg_pd, alleliccapseg_skew = convert_model_segments_to_alleliccapseg(model_segments_seg_pd,
                                                                                           model_segments_af_param_pd)

        #write the results
        alleliccapseg_seg_pd.to_csv(alleliccapseg_seg_output_file, sep='\t', index=False, na_rep='NaN')
        np.savetxt(alleliccapseg_skew_output_file, alleliccapseg_skew)

rule igv_convert:
    input:
        "results/cnv/merge_annotation/{pair_barcode}.merged.seg"
    output:
        tmp = temp("results/cnv/igv_convert/{pair_barcode}.igv.tmp"),
        tmp2 = temp("results/cnv/igv_convert/{pair_barcode}.igv.tmp2"),
        seg = "results/cnv/igv_convert/{pair_barcode}.igv.seg"
    params:
        comment_char = "@",
        field = "SAMPLE",
        segment_mean_col = "MEAN_LOG2_COPY_RATIO",
        value = lambda wildcards: "{aliquot_barcode}.called.seg".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode)),
        mem = CLUSTER_META["igv_convert"]["mem"]
    threads:
        CLUSTER_META["igv_convert"]["ppn"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/cnv/igv_convert/{pair_barcode}.log"
    benchmark:
        "benchmarks/cnv/igv_convert/{pair_barcode}.txt"
    message:
        "Convert merged segments for IGV use\n"
        "Pair ID: {wildcards.pair_barcode}"
    shell:"""
        # Modified from original task written by Chip Stewart
        grep -v "^"{params.comment_char} {input} > tmp.tsv
        head -1 tmp.tsv > tmp.header.txt
        cat tmp.header.txt | while IFS=$'\n\r' read -r line
        do
                echo {params.field} > tmp.x.tsv
        done
        sed 1,1d tmp.tsv > tmp.rest.txt

        cat tmp.rest.txt | while IFS=$'\n\r' read -r line
        do
            echo {params.value} >> tmp.x.tsv
        done

        paste tmp.x.tsv tmp.tsv > {output.tmp}
        head -1 {output.tmp} > tmp_header2.txt

        tr "\t" "\n" < tmp_header2.txt | grep -n {params.segment_mean_col} | cut -f1 -d: > tmp_col_num
        COL_NUM=`cat tmp_col_num`
        echo $COL_NUM

        cut -f$COL_NUM  {output.tmp} > col_data
        cut --complement -f $COL_NUM {output.tmp} > {output.tmp2}
        paste {output.tmp2} col_data > {output.seg}
        """

rule gistic_convert:
    input:
        "results/cnv/igv_convert/{pair_barcode}.igv.seg"
    output:
        "results/cnv/gistic_convert/{pair_barcode}.gistic2.seg"
    params:
        mem = CLUSTER_META["gistic_convert"]["mem"]
    threads:
        CLUSTER_META["gistic_convert"]["ppn"]
    #conda:
    #    "../envs/gatk4.yaml"
    log:
        "logs/cnv/gistic_convert/{pair_barcode}.log"
    benchmark:
        "benchmarks/cnv/gistic_convert/{pair_barcode}.txt"
    message:
        "Convert for GISTIC 2\n"
        "Pair ID: {wildcards.pair_barcode}"
    run:
        import csv
        input_file = input[0]
        output_file = output[0]

        """
        The column headers are:
        (1)  Sample           (sample name)
        (2)  Chromosome  (chromosome number)
        (3)  Start Position  (segment start position, in bases)
        (4)  End Position   (segment end position, in bases)
        (5)  Num markers      (number of markers in segment)
        (6)  Seg.CN       (log2() -1 of copy number)
        """
            
        with open(input_file, 'r') as tsvinfp, open(output_file, 'w') as tsvoutfp:
            tsvin = csv.DictReader(tsvinfp, delimiter='\t')
            tsvout = csv.writer(tsvoutfp, delimiter="\t")
            for r in tsvin:
                int_ify_num_points = r["NUM_POINTS_COPY_RATIO_1"].replace(".0", "")
                outrow = [r["SAMPLE"], r["CONTIG"], r["START"], r["END"], int_ify_num_points, r["MEAN_LOG2_COPY_RATIO"]]
                print(outrow)
                tsvout.writerow(outrow)

rule runabsolute:
    input:
        seg = "results/cnv/acs_convert/{pair_barcode}.acs.seg",
        skew = "results/cnv/acs_convert/{pair_barcode}.skew",
        maf = "results/mutect2/vcf2maf/{pair_barcode}.final.maf"
    output:
        res = "results/cnv/acs_convert/{pair_barcode}/{pair_barcode}.ABSOLUTE_mode.res.Rds",
        tab = "results/cnv/acs_convert/{pair_barcode}/{pair_barcode}.ABSOLUTE_mode.tab.Rds",
        rda = "results/cnv/acs_convert/{pair_barcode}/{pair_barcode}.ABSOLUTE.RData",
        pdf = "results/cnv/acs_convert/{pair_barcode}/{pair_barcode}.ABSOLUTE_plot.pdf"
    params:
        outdir = "results/cnv/acs_convert/{pair_barcode}/"
    threads:
        CLUSTER_META["runabsolute"]["ppn"]
    conda:
        "../envs/absolute.yaml"
    log:
        "logs/cnv/runabsolute/{pair_barcode}.log"
    benchmark:
        "benchmarks/cnv/runabsolute/{pair_barcode}.txt"
    message:
        "Run ABSOLUTE 1.2 (patched)\n"
        "Pair ID: {wildcards.pair_barcode}"
    shell:"""
        Rscript -e \
            "if (!require('ABSOLUTE')) devtools::install_github('TheJacksonLaboratory/Broad-ABSOLUTE');
             RunAbsolute({input.seg}, 
             genome_build = 'hg19',
             platform = 'Illumina_WES',
             copy_num_type = 'allelic',
             results.dir = '{params.outdir}',
             sample.name = '{wildcards.pair_barcode}',
             gender = NA,
             output.fn.base = '{wildcards.pair_barcode}',
             maf.fn = '{input.maf}',
             min.mut.af = 0.05,
             SSNV_skew = as.numeric(readLines('{input.skew}',warn=F)),
             verbose = T,
             max.as.seg.count = 10E10,
             primary.disease='Glioma')
    """
    
## END ##
