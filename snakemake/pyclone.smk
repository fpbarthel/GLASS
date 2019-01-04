## BEGIN

rule pyclone_tsv:
    """
    PyClone TSV rule
    """
    output:
        tsv = "results/pyclone/mytsv/{aliquot_barcode}.tsv"
    params:
        mem = CLUSTER_META["pyclone_tsv"]["mem"]
    threads:
        CLUSTER_META["pyclone_tsv"]["ppn"]
    conda:
        "../envs/r.yaml"
    log:
        "logs/pyclone/tsv/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/pyclone/tsv/{aliquot_barcode}.txt"
    message:
        "Create PyClone input TSV file\n"
        "Aliquot: {wildcards.aliquot_barcode}"
    script:
        "../R/snv/pyclone_create_tsv.R"

rule pyclone_setup:
    """
    PyClone setup
    """
    input:
        lambda wildcards: expand("results/pyclone/mytsv/{barcode}.tsv", barcode = manifest.getPyCloneAliquots(wildcards.case_barcode, wildcards.analysis_type))
    output:
        "results/pyclone/run/{case_barcode}-{analysis_type}/config.yaml"
    params:
        mem = CLUSTER_META["pyclone_setup"]["mem"],
        samples = lambda wildcards: " ".join(manifest.getPyCloneAliquots(wildcards.case_barcode, wildcards.analysis_type)),
        purity = lambda wildcards: " ".join(manifest.getPyClonePurity(wildcards.case_barcode, wildcards.analysis_type)),
        workdir = "results/pyclone/run/{case_barcode}-{analysis_type}"
    threads:
        CLUSTER_META["pyclone_setup"]["ppn"]
    conda:
        "../envs/pyclone.yaml"
    log:
        "logs/pyclone/setup/{case_barcode}.{analysis_type}.log"
    benchmark:
        "benchmarks/pyclone/setup/{case_barcode}.{analysis_type}.txt"
    message:
        "Setup PyClone\n"
        "Analysis type: {wildcards.analysis_type}\n"
        "Case: {wildcards.case_barcode}"
    shell:
        "PyClone setup_analysis \
            --in_files {input} \
            --working_dir {params.workdir} \
            --tumour_contents {params.purity} \
            --samples {params.samples} \
            --density pyclone_binomial \
            --init_method connected \
            --num_iters 10000 \
            --prior parental_copy_number \
            > {log} 2>&1"

rule pyclone_run:
    """
    PyClone run
    """
    input:
        "results/pyclone/run/{case_barcode}-{analysis_type}/config.yaml"
    output:
        "results/pyclone/run/{case_barcode}-{analysis_type}/trace/alpha.tsv.bz2",
        "results/pyclone/run/{case_barcode}-{analysis_type}/trace/labels.tsv.bz2"
    params:
        mem = CLUSTER_META["pyclone_run"]["mem"]
    threads:
        CLUSTER_META["pyclone_run"]["ppn"]
    conda:
        "../envs/pyclone.yaml"
    log:
        "logs/pyclone/run/{case_barcode}.{analysis_type}.log"
    benchmark:
        "benchmarks/pyclone/run/{case_barcode}.{analysis_type}.txt"
    message:
        "Run PyClone\n"
        "Analysis type: {wildcards.analysis_type}\n"
        "Case: {wildcards.case_barcode}"
    shell:
        "PyClone run_analysis \
            --config_file {input} \
            --seed 24081987 \
            > {log} 2>&1"

rule pyclone_plot_loci:
    """
    PyClone plot loci
    """
    input:
        conf = "results/pyclone/run/{case_barcode}-{analysis_type}/config.yaml",
        alpha = "results/pyclone/run/{case_barcode}-{analysis_type}/trace/alpha.tsv.bz2",
        levels = "results/pyclone/run/{case_barcode}-{analysis_type}/trace/labels.tsv.bz2"
    output:
        "results/pyclone/run/{case_barcode}-{analysis_type}/plots/loci/{plot_type}.pdf"
    params:
        mem = CLUSTER_META["pyclone_plot_loci"]["mem"]
    threads:
        CLUSTER_META["pyclone_plot_loci"]["ppn"]
    conda:
        "../envs/pyclone.yaml"
    log:
        "logs/pyclone/plot_loci/{case_barcode}.{analysis_type}.{plot_type}.log"
    benchmark:
        "benchmarks/pyclone/plot_loci/{case_barcode}.{analysis_type}.{plot_type}.txt"
    message:
        "Plotting loci (PyClone)\n"
        "Analysis type: {wildcards.analysis_type}\n"
        "Case: {wildcards.case_barcode}\n"
        "Plot type: {wildcards.plot_type}"
    shell:
        "PyClone plot_loci \
            --config_file {input.conf} \
            --plot_file {output} \
            --plot_type {wildcards.plot_type} \
            --burnin 1000 \
            --thin 1 \
            --min_cluster_size 2 \
            --max_clusters 6 \
            > {log} 2>&1"

rule pyclone_plot_clusters:
    """
    PyClone plot clusters
    """
    input:
        conf = "results/pyclone/run/{case_barcode}-{analysis_type}/config.yaml",
        alpha = "results/pyclone/run/{case_barcode}-{analysis_type}/trace/alpha.tsv.bz2",
        levels = "results/pyclone/run/{case_barcode}-{analysis_type}/trace/labels.tsv.bz2"
    output:
        "results/pyclone/run/{case_barcode}-{analysis_type}/plots/clusters/{plot_type}.pdf"
    params:
        mem = CLUSTER_META["pyclone_plot_clusters"]["mem"]
    threads:
        CLUSTER_META["pyclone_plot_clusters"]["ppn"]
    conda:
        "../envs/pyclone.yaml"
    log:
        "logs/pyclone/plot_clusters/{case_barcode}.{analysis_type}.{plot_type}.log"
    benchmark:
        "benchmarks/pyclone/plot_clusters/{case_barcode}.{analysis_type}.{plot_type}.txt"
    message:
        "Plotting clusters (PyClone)\n"
        "Analysis type: {wildcards.analysis_type}\n"
        "Case: {wildcards.case_barcode}\n"
        "Plot type: {wildcards.plot_type}"
    shell:
        "PyClone plot_clusters \
            --config_file {input.conf} \
            --plot_file {output} \
            --plot_type {wildcards.plot_type} \
            --burnin 1000 \
            --thin 1 \
            --min_cluster_size 2 \
            --max_clusters 6 \
            > {log} 2>&1"

rule pyclone_build_table:
    """
    PyClone build table
    """
    input:
        conf = "results/pyclone/run/{case_barcode}-{analysis_type}/config.yaml",
        alpha = "results/pyclone/run/{case_barcode}-{analysis_type}/trace/alpha.tsv.bz2",
        levels = "results/pyclone/run/{case_barcode}-{analysis_type}/trace/labels.tsv.bz2"
    output:
        "results/pyclone/run/{case_barcode}-{analysis_type}/tables/{table_type}.tsv"
    params:
        mem = CLUSTER_META["pyclone_build_table"]["mem"]
    threads:
        CLUSTER_META["pyclone_build_table"]["ppn"]
    conda:
        "../envs/pyclone.yaml"
    log:
        "logs/pyclone/build_table/{case_barcode}.{analysis_type}.{table_type}.log"
    benchmark:
        "benchmarks/pyclone/build_table/{case_barcode}.{analysis_type}.{table_type}.txt"
    message:
        "Building table (PyClone)\n"
        "Analysis type: {wildcards.analysis_type}\n"
        "Case: {wildcards.case_barcode}\n"
        "Table type: {wildcards.table_type}"
    shell:
        "PyClone build_table \
            --config_file {input.conf} \
            --out_file {output} \
            --table_type {wildcards.table_type} \
            --max_clusters 6 \
            --burnin 1000 \
            --thin 1 \
            > {log} 2>&1"

## END ##