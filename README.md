# GLASS-WG

This Github is home to alignment, variant calling and analysis pipelines for the GLASS whole genomes (GLASS-WG) analysis project.

## Version 0.1a

Alignment pipeline working with most standard settings. It can take FASTQ files or TCGA UUIDs as input that are specified in `samples.json` file (see `data/ref`). TCGA UUIDS are used to download BAM files which are then reverted back to uBAM (preserving original base qualities) and split across readgroups. Multiple FASTQ files from the same sample should reside in the same directory.

Pipeline output can be found in `bqsr/`

To-Do for alignment pipeline:
- Implement rigorous QC metrics (eg. FASTQ, samtools, MultiQC)
- Enforce QC standards, eg. create reports that indicate PASS or FAIL
- Accept gzipped FASTQ files as input

General To-Do:
- SNV calling pipeline
- CNV calling pipeline
- etc..

The rule diagram is as follows:

![rulegraph](https://user-images.githubusercontent.com/9220167/39890983-7795e5a0-546a-11e8-8605-cec1f15d5f50.png)

This was applied to test data (see `download/` for BAM files and `fastq/` for FASTQ files) which led to the following flow:

![dag](https://user-images.githubusercontent.com/9220167/39890999-8092bffc-546a-11e8-8e9d-dcb58aa4e78a.png)
