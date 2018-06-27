#!/bin/bash
## Script to quickly turn uBAM into FASTQ with appropriate filenames
## Tested on "_test2" per readgroup "RevertSam" output uBAM files
for i in *.bam; do ID="_test2"; SM=`samtools view -H $i | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g"`; FC=`samtools view -H $i | grep '^@RG' | sed "s/.*PU:[^_]*_[^_]*_[^_]*_[^_]*_\([^_]*\).*/\1/g"`; LN=`samtools view -H $i | grep '^@RG' | sed "s/.*PU:[^_]*_[^_]*_[^_]*_[^_]*_[^_]*_\([^#]*\).*/\1/g"`; bedtools bamtofastq -i $i -fq ${ID}_${SM}_${FC}_L${LN}_R1.fq -fq2 ${ID}_${SM}_${FC}_L${LN}_R2.fq; done