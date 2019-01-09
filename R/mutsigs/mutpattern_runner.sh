#!/bin/bash
#PBS -l nodes=1:ppn=12,walltime=70:00:00,mem=72gb

#Rscript /projects/varnf/GLASS/GLASS/R/mutsigs/mutationalpatterns.r
#Rscript /projects/varnf/GLASS/GLASS/R/mutsigs/individual_mutsig_heatmap.r
Rscript /projects/varnf/GLASS/GLASS/R/mutsigs/primary_recur_sig_comparison.r
