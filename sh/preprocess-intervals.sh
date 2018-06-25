#!/bin/bash
## GATK pre-proccess intervals for CNV calling
## See https://software.broadinstitute.org/gatk/documentation/article?id=11682

gatk PreprocessIntervals \
    -R /projects/verhaak-lab/verhaak_ref/gatk-legacy-bundles/b37/human_g1k_v37_decoy.fasta \
    --bin-length 1000 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O human_g1k_v37_decoy_CNV.interval_list