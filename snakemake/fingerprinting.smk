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
        lambda wildcards: expand("results/bqsr/{aliquot_id}.realn.mdup.bqsr.bam", aliquot_id = CASE_TO_ALIQUOTS[wildcards.case_id])
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

## END ##