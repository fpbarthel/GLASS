#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=72:00:00
#PBS -q batch
#PBS -m ea
#PBS -M lucas.lochovsky@jax.org
#PBS -j oe
#PBS -V
#PBS -N svelter_helix

input_bam=$1
output_directory=$2

mkdir -p $output_directory
svelter.py Setup --reference /home/lochol/projects/data/g1k_v37/human_g1k_v37_decoy.fasta --workdir $output_directory --support /home/lochol/projects/code/svelter/Support/hg19/
svelter.py NullModel --sample $input_bam --workdir $output_directory
svelter.py BPSearch --sample $input_bam --workdir $output_directory
svelter.py BPIntegrate --sample $input_bam --workdir $output_directory
