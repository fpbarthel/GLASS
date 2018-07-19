#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=72:00:00
#PBS -q batch
#PBS -N manta_cancer_1
#PBS -m bea
#PBS -M lucas.lochovsky@jax.org
#PBS -j oe
#PBS -V

tumor_bam=$1
normal_bam=$2
output_directory=$3

cd /home/lochol/projects/code/manta-build/bin
python ./configManta.py --normalBam $normal_bam --tumorBam $tumor_bam --referenceFasta /home/lochol/projects/data/g1k_v37/human_g1k_v37_decoy.fasta --runDir $output_directory
python $output_directory/runWorkflow.py -m local -j 8
