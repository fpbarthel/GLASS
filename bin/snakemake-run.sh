#!/bin/bash
#snakemake -k --jobs 999 --latency-wait 120 --max-jobs-per-second 8 --cluster-config cluster.json --cluster "qsub -N {cluster.name} -j {cluster.j} -M {cluster.M} -m {cluster.m} -l nodes={cluster.nodes}:ppn={cluster.ppn},walltime={cluster.walltime} -l mem={cluster.mem} -e {cluster.stderr} -o {cluster.stdout}"

## CLUSTETR
#snakemake --jobs 100 --latency-wait 120 --max-jobs-per-second 8 --cluster-config test_cluster.json --jobname "{jobid}.{cluster.name}" --cluster "qsub -S /bin/bash -j {cluster.j} -M {cluster.M} -m {cluster.m} -l nodes={cluster.nodes}:ppn={cluster.ppn},walltime={cluster.walltime} -l mem={cluster.mem} -e {cluster.stderr} -o {cluster.stdout}"

## DRMAA
#snakemake -p markduplicates/testT-A.realn.dedup.bam markduplicates/testT-B.realn.dedup.bam --jobs 100 --latency-wait 120 --max-jobs-per-second 8 --cluster-config cluster.json --jobname "{jobid}.{cluster.name}" --drmaa " -S /bin/bash -j {cluster.j} -M {cluster.M} -m {cluster.m} -l nodes={cluster.nodes}:ppn={cluster.ppn},walltime={cluster.walltime} -l mem={cluster.mem} -e {cluster.stderr} -o {cluster.stdout}"

CONFIGFILE="conf/config.yaml"
CLUSTRCONF="conf/cluster.json"
OPTS="-p"
RULEGRAPH=0
DAG=0


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
	    --rulegraph)
	    RULEGRAPH=1
	    shift # past argument
	    ;;
	    --dag)
	    DAG=1
	    shift # past argument
	    ;;
	    -d|--dry-run)
	    CONFIGFILE="conf/test_config.yaml"
		CLUSTRCONF="conf/test_cluster.json"
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

EXTRA_OPTS=$(printf " %s" "${POSITIONAL[@]}")
EXTRA_OPTS=${EXTRA_OPTS:1}

echo "Running snakemake with"
echo "Directory: " `pwd`
echo "TARGET: ${TARGET}"
echo "DAG: ${DAG}"
echo "RULEGRAPH: ${RULEGRAPH}"
echo "CONFIG: ${CONFIGFILE}"
echo "CLUSTER-CONFIG: ${CLUSTRCONF}"
echo "OPTS: ${OPTS}"
echo "EXTRA OPTS: ${EXTRA_OPTS}"

read -p "Continue? (y/n)" -n 1 -r
echo    # (optional) move to a new line
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    [[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1 # handle exits from shell or function but don't exit interactive shell
fi

echo "Starting snakemake..."
sleep 1
if [ ${RULEGRAPH} == 1 ];
then
	snakemake --configfile "${CONFIGFILE}" --rulegraph | dot -Tpng > rulegraph.png
elif [ ${DAG} == 1 ];
then
	snakemake --configfile "${CONFIGFILE}" --dag | dot -Tpng > dag.png
else
	snakemake ${OPTS} --jobs 400 -k --use-conda --latency-wait 120 --max-jobs-per-second 8 --configfile "${CONFIGFILE}" --cluster-config "${CLUSTRCONF}" --jobname "{jobid}.{cluster.name}" --drmaa " -S /bin/bash -j {cluster.j} -M {cluster.M} -m {cluster.m} -q {cluster.queu} -l nodes={cluster.nodes}:ppn={cluster.ppn},walltime={cluster.walltime} -l mem={cluster.mem}gb -e {cluster.stderr} -o {cluster.stdout}" $EXTRA_OPTS $TARGET
fi

## END ##
