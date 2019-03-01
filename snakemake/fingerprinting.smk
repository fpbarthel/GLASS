## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## FingerprintSample
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## This rule uses CrosscheckFingerprints to check that all readgroups in a singke sample
## come from the same individual
## See: 
## https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/picard_fingerprint_CrosscheckFingerprints.php
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule fingerprintsample:
    input:
        "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam"
    output:
        "results/fingerprinting/sample/{aliquot_barcode}.crosscheck_metrics"
    params:
        mem = CLUSTER_META["fingerprintsample"]["mem"]
    threads:
        CLUSTER_META["fingerprintsample"]["ppn"]
    log:
        "logs/fingerprinting/{aliquot_barcode}.fingerprintsample.log"
    benchmark:
        "benchmarks/fingerprinting/{aliquot_barcode}.fingerprintsample.txt"
    message:
        "Running Picard CrosscheckFingerprints to check that all readgroups in a sample come from the same individual\n"
        "Aliquot: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g CrosscheckFingerprints \
            --HAPLOTYPE_MAP {config[haplotype_map][file]} \
            --INPUT {input} \
            --OUTPUT {output} \
            > {log} 2>&1 \
            || true"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## FingerprintCase
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Running Picard CrosscheckFingerprints across multiple samples from the sample individual 
## to check for mismatches
## See: 
## https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/picard_fingerprint_CrosscheckFingerprints.php
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule fingerprintcase:
    input:
        lambda wildcards: expand("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam", aliquot_barcode = manifest.getAliquotsByCase(wildcards.case_barcode))
    output:
        "results/fingerprinting/case/{case_barcode}.crosscheck_metrics"
    params:
        mem = CLUSTER_META["fingerprintcase"]["mem"],
        samples = lambda _, input: " ".join(["--INPUT " + s for s in input])
    threads:
        CLUSTER_META["fingerprintcase"]["ppn"]
    log:
        "logs/fingerprinting/{case_barcode}.fingerprintcase.log"
    benchmark:
        "benchmarks/fingerprinting/{case_barcode}.fingerprintcase.txt"
    message:
        "Running Picard CrosscheckFingerprints across multiple samples from the sample individual to check for mismatches\n"
        "Case: {wildcards.case_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem}g CrosscheckFingerprints \
            --HAPLOTYPE_MAP {config[haplotype_map][file]} \
            --CROSSCHECK_BY SAMPLE \
            {params.samples} \
            --OUTPUT {output} \
            > {log} 2>&1 \
            || true"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## FingerprintProject
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Running Picard CrosscheckFingerprints across all samples from a single project
## to check for mismatches
## See: 
## https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/picard_fingerprint_CrosscheckFingerprints.php
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# rule fingerprintproject:
#     input:
#         lambda wildcards: expand("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam", aliquot_barcode = manifest.getAliquotsByProject(wildcards.project))
#     output:
#         "results/fingerprinting/project/{project}.crosscheck_metrics"
#     params:
#         mem = CLUSTER_META["fingerprintproject"]["mem"]
#     threads:
#         CLUSTER_META["fingerprintproject"]["ppn"]
#     log:
#         "logs/fingerprinting/{project}.fingerprintproject.log"
#     benchmark:
#         "benchmarks/fingerprinting/{project}.fingerprintproject.txt"
#     message:
#         "Running Picard CrosscheckFingerprints across all samples from a single project to check for mismatches\n"
#         "Project: {wildcards.project}"
#     run:
#         input_samples = " ".join(["--INPUT " + s for s in input])
#         shell("gatk --java-options -Xmx{params.mem}g CrosscheckFingerprints \
#             --HAPLOTYPE_MAP {config[haplotype_map]} \
#             --LOD_THRESHOLD -5 \
#             --CROSSCHECK_BY SAMPLE \
#             {input_samples} \
#             --OUTPUT {output} \
#             > {log} 2>&1 \
#             || true")

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## FingerprintAll
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Running Picard CrosscheckFingerprints across all samples in the cohort
## See: 
## https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/picard_fingerprint_CrosscheckFingerprints.php
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# rule fingerprintall:
#     input:
#         lambda wildcards: expand("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam", aliquot_barcode = manifest.getSelectedAliquots())
#     output:
#         "results/fingerprinting/GLASS.crosscheck_metrics"
#     params:
#         mem = CLUSTER_META["fingerprintall"]["mem"],
#         samples = lambda _, input: " ".join(["--INPUT " + s for s in input])
#     threads:
#         CLUSTER_META["fingerprintall"]["ppn"]
#     log:
#         "logs/fingerprinting/GLASS.fingerprintall.log"
#     benchmark:
#         "benchmarks/fingerprinting/GLASS.fingerprintall.txt"
#     message:
#         "Running Picard CrosscheckFingerprints across entire cohort"
#     shell:
#         "gatk --java-options -Xmx{params.mem}g CrosscheckFingerprints \
#             --HAPLOTYPE_MAP {config[haplotype_map][file]} \
#             --LOD_THRESHOLD -5 \
#             --CROSSCHECK_BY SAMPLE \
#             {params.samples} \
#             --OUTPUT {output} \
#             > {log} 2>&1 \
#             || true"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ClusterFingerprintBatch
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Running Picard ClusterCrosscheckMetrics on samples from one batch
## See: 
## https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/picard_fingerprint_CrosscheckFingerprints.php
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# rule clusterfingerprintproject:
#     input:
#         "results/fingerprinting/project/{project}.crosscheck_metrics"
#     output:
#         "results/fingerprinting/project/{project}.clustered.crosscheck_metrics"
#     params:
#         mem = CLUSTER_META["clusterfingerprintproject"]["mem"]
#     threads:
#         CLUSTER_META["clusterfingerprintproject"]["ppn"]
#     log:
#         "logs/fingerprinting/{project}.clusterfingerprintproject.log"
#     benchmark:
#         "benchmarks/fingerprinting/{project}.clusterfingerprintproject.txt"
#     message:
#         "Running Picard ClusterCrosscheckMetrics on project\n"
#         "Project: {wildcards.project}"
#     shell:
#         "gatk --java-options -Xmx{params.mem}g ClusterCrosscheckMetrics \
#             --INPUT {input} \
#             --LOD_THRESHOLD 5 \
#             --OUTPUT {output} \
#             > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ClusterFingerprintAll
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Running Picard ClusterCrosscheckMetrics on all samples
## See: 
## https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/picard_fingerprint_CrosscheckFingerprints.php
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# rule clusterfingerprintall:
#     input:
#         "results/fingerprinting/GLASS-.crosscheck_metrics"
#     output:
#         "results/fingerprinting/GLASS-WG.clustered.crosscheck_metrics"
#     params:
#         mem = CLUSTER_META["clusterfingerprintall"]["mem"]
#     threads:
#         CLUSTER_META["clusterfingerprintall"]["ppn"]
#     log:
#         "logs/fingerprinting/GLASS-WG.clusterfingerprintall.log"
#     benchmark:
#         "benchmarks/fingerprinting/GLASS-WG.clusterfingerprintall.txt"
#     message:
#         "Running Picard ClusterCrosscheckMetrics on entire GLASS-WG cohort"
#     shell:
#         "gatk --java-options -Xmx{params.mem}g ClusterCrosscheckMetrics \
#             --INPUT {input} \
#             --LOD_THRESHOLD 5 \
#             --OUTPUT {output} \
#             > {log} 2>&1"

## END ##