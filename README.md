# GLASS-WG

This Github is home to alignment, variant calling and analysis pipelines for the GLASS whole genomes (GLASS-WG) analysis project.

The "home" directory for the project on Helix was moved to `/fastscratch/verhaak-lab/GLASS-WG`. This will always contain a mirror of the `master` branch. Keep your own testing/development branches (eg. `devel`) under a personal directory, eg. I keep mine under `/projects/barthf/GLASS-WG`.

## Changelog

### Version 0.3
- Finalized JSON format, see `data/manifest` for details
- Alignment pipeline works successfully for BAM and FQ input, both TCGA and HK dataset processed succesfully
- Updated reference genome to b37-decoy
- SNV pipeline works (Mutect2)
- To-Do: CNV, SV, alternative SNV
- Cleaned up repository

### Version 0.2

- Improved JSON format and parsing
- JSON now indicates patients, samples, files and readgroups
- Added QC steps
- Streamlined workflow
- Support for BAM and FASTQ as input
- No longer relies on dynamic() statement in snakemake and instead statically reads readgroups from JSON (this is difficult for TCGA as BAM files need to be downloaded before readgroups can be determined, thus download of TCGA files should now be considered a seperate step. Regardless, determination of readgroups in TCGA BAMs has been completed and are stored in `data/` directory).
- Fixed a lot of errors/bugsb

### Version 0.1

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
