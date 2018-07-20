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
        "results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam"
    output:
        "results/fingerprinting/sample/{aliquot_id}.crosscheck_metrics"
    params:
        mem = CLUSTER_META["fingerprintsample"]["mem"]
    threads:
        CLUSTER_META["fingerprintsample"]["ppn"]
    log:
        "logs/fingerprinting/{aliquot_id}.fingerprintsample.log"
    benchmark:
        "benchmarks/fingerprinting/{aliquot_id}.fingerprintsample.txt"
    message:
        "Running Picard CrosscheckFingerprints to check that all readgroups in a sample come from the same individual\n"
        "Aliquot: {wildcards.aliquot_id}"
    shell:
        "gatk --java-options -Xmx{params.mem}g CrosscheckFingerprints \
            --HAPLOTYPE_MAP {config[haplotype_map]} \
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
        lambda wildcards: expand("results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam", aliquot_id = CASE_TO_ALIQUOT[wildcards.case_id])
    output:
        "results/fingerprinting/case/{case_id}.crosscheck_metrics"
    params:
        mem = CLUSTER_META["fingerprintcase"]["mem"]
    threads:
        CLUSTER_META["fingerprintcase"]["ppn"]
    log:
        "logs/fingerprinting/{case_id}.fingerprintcase.log"
    benchmark:
        "benchmarks/fingerprinting/{case_id}.fingerprintcase.txt"
    message:
        "Running Picard CrosscheckFingerprints across multiple samples from the sample individual to check for mismatches\n"
        "Case: {wildcards.case_id}"
    run:
        input_samples = " ".join(["--INPUT " + s for s in input])
        shell("gatk --java-options -Xmx{params.mem}g CrosscheckFingerprints \
            --HAPLOTYPE_MAP {config[haplotype_map]} \
            {input_samples} \
            --OUTPUT {output} \
            > {log} 2>&1 \
            || true")

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## FingerprintBatch
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Running Picard CrosscheckFingerprints across all samples from a single batch
## to check for mismatches
## See: 
## https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/picard_fingerprint_CrosscheckFingerprints.php
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule fingerprintbatch:
    input:
        lambda wildcards: expand("results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam", aliquot_id = BATCH_TO_ALIQUOT[wildcards.batch])
    output:
        "results/fingerprinting/batch/{batch}.crosscheck_metrics"
    params:
        mem = CLUSTER_META["fingerprintbatch"]["mem"]
    threads:
        CLUSTER_META["fingerprintbatch"]["ppn"]
    log:
        "logs/fingerprinting/{batch}.fingerprintbatch.log"
    benchmark:
        "benchmarks/fingerprinting/{batch}.fingerprintbatch.txt"
    message:
        "Running Picard CrosscheckFingerprints across all samples from a single batch to check for mismatches\n"
        "Batch: {wildcards.batch}"
    run:
        input_samples = " ".join(["--INPUT " + s for s in input])
        shell("gatk --java-options -Xmx{params.mem}g CrosscheckFingerprints \
            --HAPLOTYPE_MAP {config[haplotype_map]} \
            --LOD_THRESHOLD -5 \
            --CROSSCHECK_BY SAMPLE \
            {input_samples} \
            --OUTPUT {output} \
            > {log} 2>&1 \
            || true")

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## FingerprintCohort
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Running Picard CrosscheckFingerprints across all samples in the cohort
## See: 
## https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/picard_fingerprint_CrosscheckFingerprints.php
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule fingerprintcohort:
    input:
        lambda wildcards: expand("results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam", aliquot_id = ALIQUOT_TO_READGROUP.keys())
    output:
        "results/fingerprinting/GLASS-WG.crosscheck_metrics"
    params:
        mem = CLUSTER_META["fingerprintcohort"]["mem"]
    threads:
        CLUSTER_META["fingerprintcohort"]["ppn"]
    log:
        "logs/fingerprinting/GLASS-WG.fingerprintcohort.log"
    benchmark:
        "benchmarks/fingerprinting/GLASS-WG.fingerprintcohort.txt"
    message:
        "Running Picard CrosscheckFingerprints across entire cohort"
    run:
        input_samples = " ".join(["--INPUT " + s for s in input])
        shell("gatk --java-options -Xmx{params.mem}g CrosscheckFingerprints \
            --HAPLOTYPE_MAP {config[haplotype_map]} \
            --LOD_THRESHOLD -5 \
            --CROSSCHECK_BY SAMPLE \
            {input_samples} \
            --OUTPUT {output} \
            > {log} 2>&1 \
            || true")

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ClusterFingerprintBatch
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Running Picard ClusterCrosscheckMetrics on samples from one batch
## See: 
## https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/picard_fingerprint_CrosscheckFingerprints.php
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule clusterfingerprintbatch:
    input:
        "results/fingerprinting/batch/{batch}.crosscheck_metrics"
    output:
        "results/fingerprinting/batch/{batch}.clustered.crosscheck_metrics"
    params:
        mem = CLUSTER_META["clusterfingerprintbatch"]["mem"]
    threads:
        CLUSTER_META["clusterfingerprintbatch"]["ppn"]
    log:
        "logs/fingerprinting/{batch}.clusterfingerprintbatch.log"
    benchmark:
        "benchmarks/fingerprinting/{batch}.clusterfingerprintbatch.txt"
    message:
        "Running Picard ClusterCrosscheckMetrics on batch\n"
        "Batch: {wildcards.batch}"
    shell:
        "gatk --java-options -Xmx{params.mem}g ClusterCrosscheckMetrics \
            --INPUT {input} \
            --LOD_THRESHOLD 5 \
            --OUTPUT {output} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ClusterFingerprintCohort
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Running Picard ClusterCrosscheckMetrics on all samples
## See: 
## https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/picard_fingerprint_CrosscheckFingerprints.php
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule clusterfingerprintcohort:
    input:
        "results/fingerprinting/GLASS-WG.crosscheck_metrics"
    output:
        "results/fingerprinting/GLASS-WG.clustered.crosscheck_metrics"
    params:
        mem = CLUSTER_META["clusterfingerprintcohort"]["mem"]
    threads:
        CLUSTER_META["clusterfingerprintcohort"]["ppn"]
    log:
        "logs/fingerprinting/GLASS-WG.clusterfingerprintcohort.log"
    benchmark:
        "benchmarks/fingerprinting/GLASS-WG.clusterfingerprintcohort.txt"
    message:
        "Running Picard ClusterCrosscheckMetrics on entire GLASS-WG cohort"
    shell:
        "gatk --java-options -Xmx{params.mem}g ClusterCrosscheckMetrics \
            --INPUT {input} \
            --LOD_THRESHOLD 5 \
            --OUTPUT {output} \
            > {log} 2>&1"

## END ##