#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=144:00:00
#PBS -q long
#PBS -N delly_cancer_helix
#PBS -m bea
#PBS -M lucas.lochovsky@jax.org
#PBS -j oe
#PBS -V

tumor_bam=$1
normal_bam=$2
output_file=$3

cd /home/lochol/code/SVE/src/delly/src
./delly call -x /home/lochol/code/SVE/scripts/FusorSV/data/human_g1k_v37_decoy_svmask.json -o $output_file -g /home/lochol/projects/data/g1k_v37/human_g1k_v37_decoy.fasta $tumor_bam $normal_bam
