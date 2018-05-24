#!/bin/sh
## Get readgroups
BAMDIR="/fastscratch/barthf/GLASS-WG/download"
find "${BAMDIR}" -maxdepth 2 -type f -name "*bam*" | xargs -I% sh -c "samtools view -H % | grep ^@RG | sed 's|^|%\t|' | grep -v '^\['"