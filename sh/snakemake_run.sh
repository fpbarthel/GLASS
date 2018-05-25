#!/bin/bash
#snakemake -k --jobs 999 --latency-wait 120 --max-jobs-per-second 8 --cluster-config cluster.json --cluster "qsub -N {cluster.name} -j {cluster.j} -M {cluster.M} -m {cluster.m} -l nodes={cluster.nodes}:ppn={cluster.ppn},walltime={cluster.walltime} -l mem={cluster.mem} -e {cluster.stderr} -o {cluster.stdout}"

## CLUSTETR
#snakemake --jobs 100 --latency-wait 120 --max-jobs-per-second 8 --cluster-config test_cluster.json --jobname "{jobid}.{cluster.name}" --cluster "qsub -S /bin/bash -j {cluster.j} -M {cluster.M} -m {cluster.m} -l nodes={cluster.nodes}:ppn={cluster.ppn},walltime={cluster.walltime} -l mem={cluster.mem} -e {cluster.stderr} -o {cluster.stdout}"

## DRMAA
snakemake qc/testT-A.ValidateSamFile.txt qc/testT-B.ValidateSamFile.txt qc/testT-A.WgsMetrics.txt qc/testT-B.WgsMetrics.txt --jobs 100 --latency-wait 120 --max-jobs-per-second 8 --cluster-config cluster.json --jobname "{jobid}.{cluster.name}" --drmaa " -S /bin/bash -j {cluster.j} -M {cluster.M} -m {cluster.m} -l nodes={cluster.nodes}:ppn={cluster.ppn},walltime={cluster.walltime} -l mem={cluster.mem} -e {cluster.stderr} -o {cluster.stdout}"