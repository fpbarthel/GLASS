#!/bin/bash

#PBS -N svaba-test
#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=15
#PBS -r n
#PBS -M kevin.c.johnson@jax.org
#PBS -m a
#PBS -k oe
#PBS -q batch
#PBS -o /projects/verhaak-lab/sulman_GSCs/GSC-BAM/svaba-test/logs/svaba-${PBS_JOBID}.log
#PBS -e /projects/verhaak-lab/sulman_GSCs/GSC-BAM/svaba-test/logs/svaba-${PBS_JOBID}.err
#PBS -V

# Example bam file from the Sulman GSC project using a matched normal.
TUM_BAM="/projects/verhaak-lab/sulman_GSCs/GSC-BAM/BAM/GS6-27-sample1_S1.aln.dup.realn.recal.rp.bam"
NORM_BAM="/projects/verhaak-lab/sulman_GSCs/GSC-BAM/BAM/N6-27-sample2_S2.aln.dup.realn.recal.rp.bam"
DBSNP="/projects/verhaak-lab/sulman_GSCs/GSC-BAM/svaba-test/test-reference/dbsnp_indel.vcf"
CORES=10
REF="/projects/verhaak-lab/glassdir/dockscratch/bundle/human_g1k_v37_decoy.fasta"

# Samir downloaded svaba and saved it to Verhaak_env.
module load rvsvaba

# Change our working directory to the SVABA test.
cd $PBS_O_WORKDIR

# Set the date for which a sample was run.
STARTTIME=`date`
echo $STARTTIME
echo "Processing $TUM_BAM"

# Run on all chromosomes, including scaffolds.
svaba run -t $TUM_BAM -n $NORM_BAM -D $DBSNP -a GS6-27 -G $REF --hp -p $CORES