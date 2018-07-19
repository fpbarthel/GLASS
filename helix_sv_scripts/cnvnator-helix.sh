#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=72:00:00
#PBS -q batch
#PBS -m ea
#PBS -M lucas.lochovsky@jax.org
#PBS -j oe
#PBS -V
#PBS -N cnvnator_helix

input_bam=$1
output_directory=$2

cd $output_directory
/home/lochol/code/SVE/bin/sve call -r /home/lochol/projects/data/g1k_v37/human_g1k_v37_decoy.fasta -g hg19 -a cnvnator $input_bam
