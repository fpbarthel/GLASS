### Generating manifest files for Snakemake pipeline

There were a large number of differences in how metadata was stored for each of the GLASS cohorts. We sought to standardize metadata and sequencing information. Each R script in this directory represents our attempt to wrangle the data into a structure that works with the Snakemake pipeline. We later migrated to a PostgreSQL format, but these scripts may be helpful to others trying to implement the GLASS workflow.
