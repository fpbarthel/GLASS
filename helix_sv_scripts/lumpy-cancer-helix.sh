#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=72:00:00
#PBS -l mem=128GB
#PBS -q batch
#PBS -N lumpy_cancer_2
#PBS -m bea
#PBS -M lucas.lochovsky@jax.org
#PBS -j oe
#PBS -V

tumor_bam=$1
normal_bam=$2
tumor_splitter_bam=$3
normal_splitter_bam=$4
tumor_discordant_bam=$5
normal_discordant_bam=$6
output_file=$7

module load python/2.7.3
cd /home/lochol/projects/tools/lumpy-sv/bin
./lumpyexpress -B $tumor_bam,$normal_bam -S $tumor_splitter_bam,$normal_splitter_bam -D $tumor_discordant_bam,$normal_discordant_bam -o $output_file
