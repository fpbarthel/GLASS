## CNV
## postproessing

# rule combinetracks:
#     input:
#         tumor_called_seg = lambda wildcards: "results/cnv/callsegments/{aliquot_barcode}.called.seg".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode)),
#         normal_called_seg = lambda wildcards: "results/cnv/callsegments/{aliquot_barcode}.called.seg".format(aliquot_barcode = manifest.getNormal(wildcards.pair_barcode))
#     output:
#         germline_tagged_seg = temp("results/cnv/combinetracks/{pair_barcode}.germline_tagged.seg"),
#         centromere_tagged_seg = temp("results/cnv/combinetracks/{pair_barcode}.centromere_tagged.seg"),
#         final_seg = "results/cnv/combinetracks/{pair_barcode}.final.seg"
#     params:
#         mem = CLUSTER_META["combinetracks"]["mem"]
#     threads:
#         CLUSTER_META["combinetracks"]["ppn"]
#     conda:
#         "../envs/gatk4.yaml"
#     log:
#         "logs/cnv/combinetracks/{pair_barcode}.log"
#     benchmark:
#         "benchmarks/cnv/combinetracks/{pair_barcode}.txt"
#     message:
#         "Combine tumor and normal segmentation\n"
#         "Pair ID: {wildcards.pair_barcode}"
#     shell:
#         "echo \"======= Germline Tagging\";"
#         "gatk --java-options -Xmx{params.mem}g TagGermlineEvents \
#             --segments {input.tumor_called_seg} \
#             --called-matched-normal-seg-file {input.normal_called_seg} \
#             -O {output.germline_tagged_seg} \
#             -R {config[reference_fasta]} \
#             > {log} 2>&1;"
        
#         "echo \"======= Centromeres\";"
#         "gatk --java-options -Xmx{params.mem}g CombineSegmentBreakpoints \
#             --segments {output.germline_tagged_seg} \
#             --segments {config[cnv][centromere]} \
#             --columns-of-interest NUM_POINTS_COPY_RATIO \
#             --columns-of-interest MEAN_LOG2_COPY_RATIO \
#             --columns-of-interest CALL \
#             --columns-of-interest POSSIBLE_GERMLINE \
#             --columns-of-interest type \
#             -O {output.centromere_tagged_seg} \
#             -R {config[reference_fasta]} \
#             >> {log} 2>&1;"
        
#         "echo \"======= GISTIC blacklist\";"
#         "gatk --java-options -Xmx{params.mem}g CombineSegmentBreakpoints \
#             --segments {output.centromere_tagged_seg} \
#             --segments {config[cnv][gistic]} \
#             --columns-of-interest NUM_POINTS_COPY_RATIO \
#             --columns-of-interest MEAN_LOG2_COPY_RATIO \
#             --columns-of-interest CALL \
#             --columns-of-interest POSSIBLE_GERMLINE \
#             --columns-of-interest type \
#             --columns-of-interest ID \
#             -O {output.final_seg} \
#             -R {config[reference_fasta]} \
#             >> {log} 2>&1;"

# rule prepare_acs:
#     input:
#         called_seg = "results/cnv/combinetracks/{pair_barcode}.final.seg",
#         modeled_seg = lambda wildcards: "results/cnv/modelsegments/{aliquot_barcode}/{aliquot_barcode}.modelFinal.seg".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode)) # "results/cnv/modelsegments/{pair_barcode}/{pair_barcode}.modelFinal.seg" # 
#     output:
#         temp("results/cnv/prepare_acs/{pair_barcode}.prepare_acs.seg")
#     params:
#         mem = CLUSTER_META["prepare_acs"]["mem"]
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
#         "echo \"======= Merging GATK Model Seg and GATK Segment caller file\";"
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
#         mem = CLUSTER_META["filter_tagged"]["mem"]
#     threads:
#         CLUSTER_META["filter_tagged"]["ppn"]
#     #conda:
#     #    "../envs/gatk4.yaml"
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
#         tumor_tagged_df = pandas.read_csv(input[0], delimiter="\t", comment="@")
#         tumor_tagged_pruned_df = tumor_tagged_df[(tumor_tagged_df["POSSIBLE_GERMLINE"] == 0) & (tumor_tagged_df["type"] != "centromere") & (tumor_tagged_df["ID"].isnull())]
#         output_filename = output[0]
#         print(output_filename)
#         tumor_tagged_pruned_df.to_csv(output_filename, sep="\t", index=False)

# rule merge_annotation:
#     input:
#         "results/cnv/filter_tagged/{pair_barcode}.pruned.seg"
#     output:
#         "results/cnv/merge_annotation/{pair_barcode}.merged.seg"
#     params:
#         mem = CLUSTER_META["merge_annotation"]["mem"]
#     threads:
#         CLUSTER_META["merge_annotation"]["ppn"]
#     #conda: ## COMMENTED OUT BECAUSE USES GATK 4.0.10.1 which is not on conda
#     #    "../envs/gatk4.yaml"
#     log:
#         "logs/cnv/merge_annotation/{pair_barcode}.log"
#     benchmark:
#         "benchmarks/cnv/merge_annotation/{pair_barcode}.txt"
#     message:
#         "Merge adjacent genomic regions when annotations match\n"
#         "Pair ID: {wildcards.pair_barcode}"
#     shell:
#         "echo \"======= Merging\";"
#         "gatk --java-options -Xmx{params.mem}g MergeAnnotatedRegionsByAnnotation \
#             --segments {input} \
#             --annotations-to-match MEAN_LOG2_COPY_RATIO \
#             --annotations-to-match LOG2_COPY_RATIO_POSTERIOR_10 \
#             --annotations-to-match LOG2_COPY_RATIO_POSTERIOR_50 \
#             --annotations-to-match LOG2_COPY_RATIO_POSTERIOR_90 \
#             --annotations-to-match MINOR_ALLELE_FRACTION_POSTERIOR_10 \
#             --annotations-to-match MINOR_ALLELE_FRACTION_POSTERIOR_50 \
#             --annotations-to-match MINOR_ALLELE_FRACTION_POSTERIOR_90 \
#             -O {output} \
#             -R {config[reference_fasta]} \
#             > {log} 2>&1"

# rule acs_convert:
#     input:
#         af_param = lambda wildcards: "results/cnv/modelsegments/{aliquot_barcode}/{aliquot_barcode}.modelFinal.af.param".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode)),
#         merged_seg = "results/cnv/merge_annotation/{pair_barcode}.merged.seg"
#     output:
#         seg = "results/cnv/acs_convert/{pair_barcode}.acs.seg",
#         skew = "results/cnv/acs_convert/{pair_barcode}.skew"
#     params:
#         mem = CLUSTER_META["acs_convert"]["mem"]
#     threads:
#         CLUSTER_META["acs_convert"]["ppn"]
#     #conda:
#     #    "../envs/gatk4.yaml"
#     log:
#         "logs/cnv/acs_convert/{pair_barcode}.log"
#     benchmark:
#         "benchmarks/cnv/acs_convert/{pair_barcode}.txt"
#     message:
#         "Produces a seg file of identical format to AllelicCapSeg\n"
#         "Pair ID: {wildcards.pair_barcode}"
#     run:
#         import sys
#         import re
#         import pandas as pd
#         import numpy as np
#         from collections import defaultdict
#         import scipy
#         from scipy import special as sp
#         import os.path

#         model_segments_seg_input_file = input["merged_seg"]
#         model_segments_af_param_input_file = input["af_param"]
#         alleliccapseg_seg_output_file = output["seg"]
#         alleliccapseg_skew_output_file = output["skew"]

#         HAM_FIST_THRESHOLD = config["cnv"]["maf90_threshold"]

#         # regular expression for matching sample name from header comment line
#         sample_name_header_regexp = "^@RG.*SM:(.*)[\t]*.*$"

#         #define AllelicCapSeg columns
#         alleliccapseg_seg_columns = [
#             'Chromosome',
#             'Start.bp',
#             'End.bp',
#             'n_probes',
#             'length',
#             'n_hets',
#             'f',
#             'tau',
#             'sigma.tau',
#             'mu.minor',
#             'sigma.minor',
#             'mu.major',
#             'sigma.major',
#             'SegLabelCNLOH']

#         def read_sample_name(input_file, max_scan_lines=10000):
#             with open(input_file, 'r') as f:
#                 for _ in range(max_scan_lines):
#                     line = f.readline()
#                     match = re.search(sample_name_header_regexp, line, re.M)
#                     if match is None:
#                         continue
#                     groups = match.groups()
#                     return groups[0]
#             raise Exception("Sample name could not be found in \"{0}\"".format(input_file))

#         #read GATK ModelSegments files and perform some basic checks
#         model_segments_seg_pd = pd.read_csv(model_segments_seg_input_file,
#                                             sep='\t', comment='@', na_values='NA')
#         model_segments_af_param_pd = pd.read_csv(model_segments_af_param_input_file, sep='\t', comment='@')

#         def simple_determine_allelic_fraction(model_segments_seg_pd):
#             result = model_segments_seg_pd['MINOR_ALLELE_FRACTION_POSTERIOR_50']
#             result[model_segments_seg_pd['MINOR_ALLELE_FRACTION_POSTERIOR_90'] > HAM_FIST_THRESHOLD] = 0.5
#             return result

#         def convert_model_segments_to_alleliccapseg(model_segments_seg_pd,
#                                                     model_segments_af_param_pd):
#             alleliccapseg_seg_pd = pd.DataFrame(columns=alleliccapseg_seg_columns)

#             #The following conversions are trivial.
#             alleliccapseg_seg_pd['Chromosome'] = model_segments_seg_pd['CONTIG']
#             alleliccapseg_seg_pd['Start.bp'] = model_segments_seg_pd['START']
#             alleliccapseg_seg_pd['End.bp'] = model_segments_seg_pd['END']
#             alleliccapseg_seg_pd['n_probes'] = model_segments_seg_pd['NUM_POINTS_COPY_RATIO_1']
#             alleliccapseg_seg_pd['length'] = alleliccapseg_seg_pd['End.bp'] - alleliccapseg_seg_pd['Start.bp']
#             alleliccapseg_seg_pd['n_hets'] = model_segments_seg_pd['NUM_POINTS_ALLELE_FRACTION']

#             #ModelSegments estimates posterior credible intervals, while AllelicCapSeg performs maximum a posteriori (MAP) estimation.
#             #The copy-ratio and allele-fraction models fit by both also differ.
#             # We will attempt a rough translation of the model fits here.

#             alleliccapseg_seg_pd['f'] = simple_determine_allelic_fraction(model_segments_seg_pd)

#             alleliccapseg_seg_pd['tau'] = 2. * 2**model_segments_seg_pd['LOG2_COPY_RATIO_POSTERIOR_50']
#             alleliccapseg_seg_pd['sigma.tau'] = 2**model_segments_seg_pd['LOG2_COPY_RATIO_POSTERIOR_90'] - 2**model_segments_seg_pd['LOG2_COPY_RATIO_POSTERIOR_10']
#             sigma_f = (model_segments_seg_pd['MINOR_ALLELE_FRACTION_POSTERIOR_90'].values - model_segments_seg_pd['MINOR_ALLELE_FRACTION_POSTERIOR_10'].values) / 2.
#             sigma_mu = np.sqrt(sigma_f**2 + alleliccapseg_seg_pd['sigma.tau']**2) #we propagate errors in the products f * tau and (1 - f) * tau in the usual way
#             alleliccapseg_seg_pd['mu.minor'] = alleliccapseg_seg_pd['f'] * alleliccapseg_seg_pd['tau']
#             alleliccapseg_seg_pd['sigma.minor'] = sigma_mu
#             alleliccapseg_seg_pd['mu.major'] = (1. - alleliccapseg_seg_pd['f']) * alleliccapseg_seg_pd['tau']
#             alleliccapseg_seg_pd['sigma.major'] = sigma_mu

#             #For whatever reason, AllelicCapSeg attempts to call CNLOH.  Documentation is spotty, but it seems like it attempts
#             # to distinguish between three states ("0 is flanked on both sides, 1 is one side, 2 is no cn.loh").
#             # Let's just set everything to 2 for now.
#             # Hopefully, ABSOLUTE is robust to this ...
#             alleliccapseg_seg_pd['SegLabelCNLOH'] = 2

#             #One important caveat: for segments with less than 10 hets, AllelicCapSeg also tries to call whether a segment is "split" or not.
#             #  This script will attempt to call "split" on all segments.
#             # ACS performs a simple hypothesis test on the alternate-allele fractions to see if
#             # a unimodal distribution peaked at 0.5 is supported over a bimodal distribution peaked at f and 1 - f.
#             # If the former is supported, then AllelicCapSeg ignores the MAP estimate of f and simply sets it to be 0.5.
#             # ABSOLUTE may actually be rather sensitive to this.  Again, let's ignore for now, and we can later port this
#             # statistical test if necessary.

#             #Finally, I believe that ABSOLUTE requires the value of the "skew" parameter from the AllelicCapSeg
#             #allele-fraction model.  This parameter is supposed to allow the model to account for reference bias,
#             #  but the model likelihood that AllelicCapSeg uses is not valid over the entire range of the skew parameter.
#             # We corrected this during the development of AllelicCNV and retain the same corrected model in ModelSegments.
#             # We will try to transform the relevant parameter in the corrected model back to a "skew",
#             # but this operation is ill defined.  Luckily, for WGS, the reference bias is typically negligible.
#             model_segments_reference_bias = model_segments_af_param_pd[
#                 model_segments_af_param_pd['PARAMETER_NAME'] == 'MEAN_BIAS']['POSTERIOR_50']
#             alleliccapseg_skew = 2. / (1. + model_segments_reference_bias)

#             return alleliccapseg_seg_pd, alleliccapseg_skew


#         #do the conversion
#         alleliccapseg_seg_pd, alleliccapseg_skew = convert_model_segments_to_alleliccapseg(model_segments_seg_pd,
#                                                                                            model_segments_af_param_pd)

#         #write the results
#         alleliccapseg_seg_pd.to_csv(alleliccapseg_seg_output_file, sep='\t', index=False, na_rep='NaN')
#         np.savetxt(alleliccapseg_skew_output_file, alleliccapseg_skew)

# rule igv_convert:
#     input:
#         "results/cnv/merge_annotation/{pair_barcode}.merged.seg"
#     output:
#         tmp1 = temp("results/cnv/igv_convert/{pair_barcode}.igv.tmp1"),
#         tmp2 = temp("results/cnv/igv_convert/{pair_barcode}.igv.tmp2"),
#         tmp3 = temp("results/cnv/igv_convert/{pair_barcode}.igv.tmp3"),
#         tmp4 = temp("results/cnv/igv_convert/{pair_barcode}.igv.tmp4"),
#         tmp5 = temp("results/cnv/igv_convert/{pair_barcode}.igv.tmp5"),
#         tmp6 = temp("results/cnv/igv_convert/{pair_barcode}.igv.tmp6"),
#         tmp7 = temp("results/cnv/igv_convert/{pair_barcode}.igv.tmp7"),
#         tmp8 = temp("results/cnv/igv_convert/{pair_barcode}.igv.tmp8"),
#         tmp9 = temp("results/cnv/igv_convert/{pair_barcode}.igv.tmp9"),
#         seg = "results/cnv/igv_convert/{pair_barcode}.igv.seg"
#     params:
#         comment_char = "@",
#         field = "SAMPLE",
#         segment_mean_col = "MEAN_LOG2_COPY_RATIO",
#         value = lambda wildcards: "{aliquot_barcode}.called.seg".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode)),
#         mem = CLUSTER_META["igv_convert"]["mem"]
#     threads:
#         CLUSTER_META["igv_convert"]["ppn"]
#     conda:
#         "../envs/gatk4.yaml"
#     log:
#         "logs/cnv/igv_convert/{pair_barcode}.log"
#     benchmark:
#         "benchmarks/cnv/igv_convert/{pair_barcode}.txt"
#     message:
#         "Convert merged segments for IGV use\n"
#         "Pair ID: {wildcards.pair_barcode}"
#     shell:"""
#         # Modified from original task written by Chip Stewart
#         grep -v "^"{params.comment_char} {input} > {output.tmp1}
#         head -1 {output.tmp1} > {output.tmp2}
#         cat {output.tmp2} | while IFS=$'\n\r' read -r line
#         do
#                 echo {params.field} > {output.tmp3}
#         done
#         sed 1,1d {output.tmp1} > {output.tmp4}

#         cat {output.tmp4} | while IFS=$'\n\r' read -r line
#         do
#             echo {params.value} >> {output.tmp3}
#         done

#         paste {output.tmp3} {output.tmp1} > {output.tmp5}
#         head -1 {output.tmp5} > {output.tmp6}

#         tr "\t" "\n" < {output.tmp6} | grep -n {params.segment_mean_col} | cut -f1 -d: > {output.tmp7}
#         COL_NUM=`cat {output.tmp7}`
#         echo $COL_NUM

#         cut -f$COL_NUM  {output.tmp5} > {output.tmp8}
#         cut --complement -f $COL_NUM {output.tmp5} > {output.tmp9}
#         paste {output.tmp9} {output.tmp8} > {output.seg}
#         """

# rule gistic_convert:
#     input:
#         "results/cnv/igv_convert/{pair_barcode}.igv.seg"
#     output:
#         "results/cnv/gistic_convert/{pair_barcode}.gistic2.seg"
#     params:
#         mem = CLUSTER_META["gistic_convert"]["mem"]
#     threads:
#         CLUSTER_META["gistic_convert"]["ppn"]
#     #conda:
#     #    "../envs/gatk4.yaml"
#     log:
#         "logs/cnv/gistic_convert/{pair_barcode}.log"
#     benchmark:
#         "benchmarks/cnv/gistic_convert/{pair_barcode}.txt"
#     message:
#         "Convert for GISTIC 2\n"
#         "Pair ID: {wildcards.pair_barcode}"
#     run:
#         import csv
#         input_file = input[0]
#         output_file = output[0]

#         """
#         The column headers are:
#         (1)  Sample           (sample name)
#         (2)  Chromosome  (chromosome number)
#         (3)  Start Position  (segment start position, in bases)
#         (4)  End Position   (segment end position, in bases)
#         (5)  Num markers      (number of markers in segment)
#         (6)  Seg.CN       (log2() -1 of copy number)
#         """
            
#         with open(input_file, 'r') as tsvinfp, open(output_file, 'w') as tsvoutfp:
#             tsvin = csv.DictReader(tsvinfp, delimiter='\t')
#             tsvout = csv.writer(tsvoutfp, delimiter="\t")
#             for r in tsvin:
#                 int_ify_num_points = r["NUM_POINTS_COPY_RATIO_1"].replace(".0", "")
#                 outrow = [r["SAMPLE"], r["CONTIG"], r["START"], r["END"], int_ify_num_points, r["MEAN_LOG2_COPY_RATIO"]]
#                 print(outrow)
#                 tsvout.writerow(outrow)

# rule runabsolute:
#     input:
#         seg = "results/cnv/acs_convert/{pair_barcode}.acs.seg",
#         skew = "results/cnv/acs_convert/{pair_barcode}.skew",
#         maf = ancient("results/mutect2/vcf2maf/{pair_barcode}.final.maf")
#     output:
#         res = "results/cnv/absolute/{pair_barcode}/{pair_barcode}.ABSOLUTE_mode.res.Rds",
#         tab = "results/cnv/absolute/{pair_barcode}/{pair_barcode}.ABSOLUTE_mode.tab.Rds",
#         rda = "results/cnv/absolute/{pair_barcode}/{pair_barcode}.ABSOLUTE.RData",
#         pdf = "results/cnv/absolute/{pair_barcode}/{pair_barcode}.ABSOLUTE_plot.pdf"
#     params:
#         outdir = "results/cnv/absolute/{pair_barcode}/"
#     threads:
#         CLUSTER_META["runabsolute"]["ppn"]
#     conda:
#         "../envs/absolute.yaml"
#     log:
#         "logs/cnv/runabsolute/{pair_barcode}.log"
#     benchmark:
#         "benchmarks/cnv/runabsolute/{pair_barcode}.txt"
#     message:
#         "Run ABSOLUTE 1.2 (patched)\n"
#         "Pair ID: {wildcards.pair_barcode}"
#     shell:"""
#         conda activate R341
#         Rscript -e 'ABSOLUTE::RunAbsolute("{input.seg}", genome_build = "hg19", platform = "Illumina_WES", copy_num_type = "allelic", results.dir = "{params.outdir}", sample.name = "{wildcards.pair_barcode}", gender = NA, output.fn.base = "{wildcards.pair_barcode}", maf.fn = "{input.maf}", min.mut.af = 0.05, SSNV_skew = as.numeric(readLines("{input.skew}",warn=F)), verbose = T, max.as.seg.count = 10E10, primary.disease="Glioma")' > {log} 2>&1
#     """


# rule titanhets:
#     input:
#         tumor = lambda wildcards: "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode)),
#         normal = lambda wildcards: "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam".format(aliquot_barcode = manifest.getNormal(wildcards.pair_barcode))
#     output:
#         sites = temp("results/cnv/titanhets/{pair_barcode}.{chr}.sites.tsv"),
#         hets = temp("results/cnv/titanhets/{pair_barcode}.{chr}.hets.tsv")
#     params:
#         basequal = 10,
#         mapqual = 10,
#         vcfqual = 100
#     threads:
#         CLUSTER_META["titanhets"]["ppn"]
#     log:
#         "logs/cnv/titanhets/{pair_barcode}.{chr}.log"
#     benchmark:
#         "benchmarks/cnv/titanhets/{pair_barcode}.{chr}.txt"
#     message:
#         "Call hets for titan\n"
#         "Pair ID: {wildcards.pair_barcode}\n"
#         "Chr: {wildcards.chr}" 
#     shell:"""
#         source activate samtools
#         samtools mpileup -uv -I -f {config[reference_fasta]} -r {wildcards.chr} -l {config[cnv][hetsites]} {input.normal} | bcftools call -v -c - | grep -e '0/1' -e '#' > {output.sites} 2> {log};
#         python python/countPysam.py {wildcards.chr} {output.sites} {input.tumor} {params.basequal} {params.mapqual} {params.vcfqual} > {output.hets} 2>> {log}
#     """

# rule preparetitan:
#     input:
#         seg     = lambda wildcards: "results/cnv/mergedenoisedreadcounts/{aliquot_barcode}.denoisedCR.tsv".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode)),
#         hets    = lambda wildcards: expand("results/cnv/titanhets/{pair_barcode}.{chr}.hets.tsv", pair_barcode = wildcards.pair_barcode, chr = list(range(1,23)) + ['X'])
#     output:
#         hets    = "results/cnv/titan/{pair_barcode}/{pair_barcode}.hets.tsv",
#         seg     = "results/cnv/titan/{pair_barcode}/{pair_barcode}.seg"
#     threads:
#         CLUSTER_META["preparetitan"]["ppn"]
#     log:
#         "logs/cnv/preparetitan/{pair_barcode}.log"
#     benchmark:
#         "benchmarks/cnv/preparetitan/{pair_barcode}.txt"
#     message:
#         "Prepare TitanCNA\n"
#         "Pair ID: {wildcards.pair_barcode}"
#     shell:"""
#         cat {input.hets} | grep -v Chr > {output.hets} 2> {log};
#         cat {input.seg} | awk -F\\t 'BEGIN{{print"chr\\tstart\\tend\\tlog2_TNratio_corrected"}} /^[^@C]/ {{ print $0 }}' > {output.seg}
#     """

rule preparetitan:
    input:
        seg     = lambda wildcards: "results/cnv/mergedenoisedreadcounts/{aliquot_barcode}.denoisedCR.tsv".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode)),
        hets    = lambda wildcards: "results/cnv/modelsegments/{aliquot_barcode}/{aliquot_barcode}.hets.tsv".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode))
    output:
        hets    = "results/cnv/titan/{pair_barcode}/{pair_barcode}.hets.tsv",
        seg     = "results/cnv/titan/{pair_barcode}/{pair_barcode}.seg"
    threads:
        CLUSTER_META["preparetitan"]["ppn"]
    log:
        "logs/cnv/preparetitan/{pair_barcode}.log"
    benchmark:
        "benchmarks/cnv/preparetitan/{pair_barcode}.txt"
    message:
        "Prepare TitanCNA\n"
        "Pair ID: {wildcards.pair_barcode}"
    shell:"""
        cat {input.hets} | awk '/^[^@]/ {{ print $1,$2,$5,$3,$6,$4 }}' | tr ' ' '\\t' > {output.hets}
        cat {input.seg} | awk -F\\t 'BEGIN{{print"chr\\tstart\\tend\\tlog2_TNratio_corrected"}} /^[^@C]/ {{ print $0 }}' > {output.seg}
    """

rule titan:
    input:
        hets    = "results/cnv/titan/{pair_barcode}/{pair_barcode}.hets.tsv",
        seg     = "results/cnv/titan/{pair_barcode}/{pair_barcode}.seg"
    output:
        seg     = "results/cnv/titan/{pair_barcode}/ploidy{ploidy}/{pair_barcode}_cluster0{cluster}.seg",
        segs    = "results/cnv/titan/{pair_barcode}/ploidy{ploidy}/{pair_barcode}_cluster0{cluster}.segs.txt",
        titan   = "results/cnv/titan/{pair_barcode}/ploidy{ploidy}/{pair_barcode}_cluster0{cluster}.titan.txt",
        params  = "results/cnv/titan/{pair_barcode}/ploidy{ploidy}/{pair_barcode}_cluster0{cluster}.params.txt",
        cf      = "results/cnv/titan/{pair_barcode}/ploidy{ploidy}/{pair_barcode}_cluster0{cluster}/{pair_barcode}_cluster0{cluster}_CF.pdf",
        cna     = "results/cnv/titan/{pair_barcode}/ploidy{ploidy}/{pair_barcode}_cluster0{cluster}/{pair_barcode}_cluster0{cluster}_CNA.pdf",
        cnaseg  = "results/cnv/titan/{pair_barcode}/ploidy{ploidy}/{pair_barcode}_cluster0{cluster}/{pair_barcode}_cluster0{cluster}_CNASEG.pdf",
        loh     = "results/cnv/titan/{pair_barcode}/ploidy{ploidy}/{pair_barcode}_cluster0{cluster}/{pair_barcode}_cluster0{cluster}_LOH.pdf",
        lohseg  = "results/cnv/titan/{pair_barcode}/ploidy{ploidy}/{pair_barcode}_cluster0{cluster}/{pair_barcode}_cluster0{cluster}_LOHSEG.pdf"
    params:
        outdir =        "results/cnv/titan/{pair_barcode}/ploidy{ploidy}/",
        alphaK =        lambda wildcards: 2500 if manifest.isExome(manifest.getTumor(wildcards.pair_barcode)) else 10000,
        alphaKHigh =    lambda wildcards: 2500 if manifest.isExome(manifest.getTumor(wildcards.pair_barcode)) else 10000,
        gender = 		lambda wildcards: manifest.getSex(manifest.getTumor(wildcards.pair_barcode)) if manifest.getSex(manifest.getTumor(wildcards.pair_barcode)) is not None else "NA"
    threads:
        CLUSTER_META["titan"]["ppn"]
    #conda:
    #    "../envs/titan.yaml"
    log:
        "logs/cnv/titan/{pair_barcode}_ploidy{ploidy}_cluster{cluster}.log"
    benchmark:
        "benchmarks/cnv/titan/{pair_barcode}_ploidy{ploidy}_cluster{cluster}.txt"
    message:
        "Run TitanCNA\n"
        "Pair ID: {wildcards.pair_barcode}\n"
        "Cluster: {wildcards.cluster}\n"
        "Ploidy: {wildcards.ploidy}"
    shell:"""
        source activate R341
        Rscript /projects/barthf/opt/TitanCNA/scripts/R_scripts/titanCNA.R \
            --id {wildcards.pair_barcode} \
            --hetFile {input.hets} \
            --cnFile {input.seg} \
            --numClusters {wildcards.cluster} \
            --numCores {threads} \
            --normal_0 0.5 \
            --ploidy_0 {wildcards.ploidy} \
            --alphaK {params.alphaK} \
            --alphaKHigh {params.alphaKHigh} \
            --minDepth 5 \
            --maxDepth 50000 \
            --gender {params.gender} \
            --estimatePloidy TRUE \
            --outDir {params.outdir} \
            --libdir /projects/barthf/opt/TitanCNA \
            > {log} 2>&1
    """

rule selecttitan:
    input:
        seg     = lambda wildcards: expand("results/cnv/titan/{pair_barcode}/ploidy{ploidy}/{pair_barcode}_cluster0{cluster}.seg", ploidy = [2,3], cluster = [1,3], pair_barcode = wildcards.pair_barcode),
        segs    = lambda wildcards: expand("results/cnv/titan/{pair_barcode}/ploidy{ploidy}/{pair_barcode}_cluster0{cluster}.segs.txt", ploidy = [2,3], cluster = [1,3], pair_barcode = wildcards.pair_barcode),
        titan   = lambda wildcards: expand("results/cnv/titan/{pair_barcode}/ploidy{ploidy}/{pair_barcode}_cluster0{cluster}.titan.txt", ploidy = [2,3], cluster = [1,3], pair_barcode = wildcards.pair_barcode),
        params  = lambda wildcards: expand("results/cnv/titan/{pair_barcode}/ploidy{ploidy}/{pair_barcode}_cluster0{cluster}.params.txt", ploidy = [2,3], cluster = [1,3], pair_barcode = wildcards.pair_barcode)
    output:
        txt     = "results/cnv/titan/{pair_barcode}/{pair_barcode}.optimalClusters.txt"
    threads:
        CLUSTER_META["selecttitan"]["ppn"]
    log:
        "logs/cnv/selecttitan/{pair_barcode}.log"
    benchmark:
        "benchmarks/cnv/selecttitan/{pair_barcode}.txt"
    message:
        "Select optimal TitanCNA cluster and ploidy\n"
        "Pair ID: {wildcards.pair_barcode}"
    shell:"""
        source activate R341
        Rscript /projects/barthf/opt/TitanCNA/scripts/R_scripts/selectSolution.R \
            --ploidyRun2=results/cnv/titan/{wildcards.pair_barcode}/ploidy2 \
            --ploidyRun3=results/cnv/titan/{wildcards.pair_barcode}/ploidy3 \
            --threshold=0.15 \
            --outFile {output.txt} \
            > {log} 2>&1
    """


rule finaltitan:
    input:
        txt     = "results/cnv/titan/{pair_barcode}/{pair_barcode}.optimalClusters.txt"
    output:
        segs    = "results/cnv/titanfinal/seg/{pair_barcode}.seg.txt",
        params  = "results/cnv/titanfinal/params/{pair_barcode}.params.txt",
        pdf     = "results/cnv/titanfinal/pdf/{pair_barcode}.pdf",
        igv		= "results/cnv/titanfinal/igv/{pair_barcode}.igv.seg"
    threads:
        CLUSTER_META["finaltitan"]["ppn"]
    log:
        "logs/cnv/finaltitan/{pair_barcode}.log"
    benchmark:
        "benchmarks/cnv/finaltitan/{pair_barcode}.txt"
    message:
        "Copy selected (final) TitanCNA results\n"
        "Pair ID: {wildcards.pair_barcode}"
    shell:"""
        TITANDIR=$(cat {input} | awk -F'\\t' '{{print $11}}' | sed -n 2p)
        cp "${{TITANDIR}}.segs.txt" {output.segs}
        cp "${{TITANDIR}}.seg" {output.igv}
        cp "${{TITANDIR}}.params.txt" {output.params}
        /opt/software/helix/ImageMagick/7.0.7-26/bin/montage $TITANDIR/*_CNA.pdf $TITANDIR/*_LOH.pdf $TITANDIR/*_CF.pdf -tile 1x3 -geometry +0+0 {output.pdf}
    """

## END ##
