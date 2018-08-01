## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Download BAM file from GDC
## GDC key needs to be re-downloaded and updated from time to time
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule download:
    output:
        "data/download/{uuid}/{filename}.bam"
    threads:
        CLUSTER_META["download"]["ppn"]
    message:
        "Downloading from GDC\n"
        "UUID {wildcards.uuid}\n"
        "File {wildcards.filename}"
    conda:
        "../envs/gdc-client.yaml"
    log:
        "logs/download/{uuid}.{filename}.log"
    benchmark:
        "benchmarks/download/{uuid}.{filename}.txt"
    shell:
        "gdc-client download \
            -d download \
            -n {threads} \
            -t {config[gdc_token]} \
            {wildcards.uuid} \
            > {log} 2>&1"