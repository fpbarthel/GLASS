#!/bin/bash
#snakemake -k --jobs 999 --latency-wait 120 --max-jobs-per-second 8 --cluster-config cluster.json --cluster "qsub -N {cluster.name} -j {cluster.j} -M {cluster.M} -m {cluster.m} -l nodes={cluster.nodes}:ppn={cluster.ppn},walltime={cluster.walltime} -l mem={cluster.mem} -e {cluster.stderr} -o {cluster.stdout}"

## CLUSTETR
#snakemake --jobs 100 --latency-wait 120 --max-jobs-per-second 8 --cluster-config test_cluster.json --jobname "{jobid}.{cluster.name}" --cluster "qsub -S /bin/bash -j {cluster.j} -M {cluster.M} -m {cluster.m} -l nodes={cluster.nodes}:ppn={cluster.ppn},walltime={cluster.walltime} -l mem={cluster.mem} -e {cluster.stderr} -o {cluster.stdout}"

## DRMAA
#snakemake -p markduplicates/testT-A.realn.dedup.bam markduplicates/testT-B.realn.dedup.bam --jobs 100 --latency-wait 120 --max-jobs-per-second 8 --cluster-config cluster.json --jobname "{jobid}.{cluster.name}" --drmaa " -S /bin/bash -j {cluster.j} -M {cluster.M} -m {cluster.m} -l nodes={cluster.nodes}:ppn={cluster.ppn},walltime={cluster.walltime} -l mem={cluster.mem} -e {cluster.stderr} -o {cluster.stdout}"

CONFIGFILE="config.yaml"
CLUSTRCONF="cluster.json"
OPTS="-p"

POSITIONAL=()
while [[ $# -gt 0 ]]
do
	key="$1"

	case $key in
	    -t|--target)
		[ -z "$2" ] && echo "No target specified" && exit 1
	    TARGET="$2"
	    shift # past argument
	    shift # past value
	    ;;
	    -d|--dry-run)
	    CONFIGFILE="test_config.yaml"
		CLUSTRCONF="test_cluster.json"
	    shift # past argument
	    ;;
	    -n|--not-real)
	    OPTS="-np"
	    shift # past argument
	    ;;*)    # unknown option
	    POSITIONAL+=("$1") # save it in an array for later
	    shift # past argument
	    ;;
	esac
done

if [ ! -f "${CONFIGFILE}" ]
then
	echo "CONFIG ${CONFIGFILE} does not exist"
	exit 1
fi

if [ ! -f "${CLUSTRCONF}" ]
then
	echo "CLUSTER-CONFIG ${CLUSTRCONF} does not exist"
	exit 1
fi

echo "Running snakemake with"
echo "Directory: " `pwd`
echo "TARGET: ${TARGET}"
echo "CONFIG: ${CONFIGFILE}"
echo "CLUSTER-CONFIG: ${CLUSTRCONF}"
echo "OPTS: ${OPTS}"

read -p "Continue? (y/n)" -n 1 -r
echo    # (optional) move to a new line
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    [[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1 # handle exits from shell or function but don't exit interactive shell
fi

echo "Starting snakemake..."
sleep 1
snakemake ${OPTS} --jobs 800 --latency-wait 120 --max-jobs-per-second 8 --configfile "${CONFIGFILE}" --cluster-config "${CLUSTRCONF}" --jobname "{jobid}.{cluster.name}" --drmaa " -S /bin/bash -j {cluster.j} -M {cluster.M} -m {cluster.m} -l nodes={cluster.nodes}:ppn={cluster.ppn},walltime={cluster.walltime} -l mem={cluster.mem}gb -e {cluster.stderr} -o {cluster.stdout}" ${TARGET}