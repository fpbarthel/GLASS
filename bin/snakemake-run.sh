#!/bin/bash
#snakemake -k --jobs 999 --latency-wait 120 --max-jobs-per-second 8 --cluster-config cluster.json --cluster "qsub -N {cluster.name} -j {cluster.j} -M {cluster.M} -m {cluster.m} -l nodes={cluster.nodes}:ppn={cluster.ppn},walltime={cluster.walltime} -l mem={cluster.mem} -e {cluster.stderr} -o {cluster.stdout}"

## CLUSTETR
#snakemake --jobs 100 --latency-wait 120 --max-jobs-per-second 8 --cluster-config test_cluster.json --jobname "{jobid}.{cluster.name}" --cluster "qsub -S /bin/bash -j {cluster.j} -M {cluster.M} -m {cluster.m} -l nodes={cluster.nodes}:ppn={cluster.ppn},walltime={cluster.walltime} -l mem={cluster.mem} -e {cluster.stderr} -o {cluster.stdout}"

## DRMAA
#snakemake -p markduplicates/testT-A.realn.dedup.bam markduplicates/testT-B.realn.dedup.bam --jobs 100 --latency-wait 120 --max-jobs-per-second 8 --cluster-config cluster.json --jobname "{jobid}.{cluster.name}" --drmaa " -S /bin/bash -j {cluster.j} -M {cluster.M} -m {cluster.m} -l nodes={cluster.nodes}:ppn={cluster.ppn},walltime={cluster.walltime} -l mem={cluster.mem} -e {cluster.stderr} -o {cluster.stdout}"

CONFIGFILE="conf/config.yaml"
CLUSTRCONF="conf/cluster.json"
NUM_JOBS=800
OPTS="-p"
RULEGRAPH=0
DAG=0
WORKDIR=`pwd`
FROMSOURCE=0
COHORT=""

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
	    -nj|--num-jobs)
		[ -z "$2" ] && echo "No number of jobs specified" && exit 1
	    NUM_JOBS=$2
	    shift # past argument
	    shift # past value
	    ;;
	    -wd|--work-dir)
		[ -z "$2" ] && echo "No workdir specified" && exit 1
	    WORKDIR="$2"
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
	    -cc|--cluster-conf)
	    [ -z "$2" ] && echo "No cluster config specified" && exit 1
	    CLUSTRCONF="$2"
	    shift # past value
	    shift # past argument
	    ;;
	    -c|--cohort)
	    [ -z "$2" ] && echo "No cohort specified" && exit 1
	    COHORT="$2"
	    shift # past value
	    shift # past argument
	    ;;
	    -n|--not-real)
	    OPTS="-np"
	    shift # past argument
	    ;;
	    -s|--source)
	    FROMSOURCE=1
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
echo "Directory: ${WORKDIR}" 
echo "TARGET: ${TARGET}"
echo "DAG: ${DAG}"
echo "RULEGRAPH: ${RULEGRAPH}"
echo "CONFIG: ${CONFIGFILE}"
echo "CLUSTER-CONFIG: ${CLUSTRCONF}"
echo "OPTS: ${OPTS}"
echo "NUM_JOBS: ${NUM_JOBS}"
echo "FROMSOURCE: ${FROMSOURCE}"
echo "COHORT: ${COHORT}"
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
	mkdir -p "${WORKDIR}/logs/drmaa"
	snakemake ${OPTS} --jobs ${NUM_JOBS} -k --use-conda --latency-wait 120 --max-jobs-per-second 2 --config workdir="${WORKDIR}" cluster_json="${CLUSTRCONF}" from_source="${FROMSOURCE}" cohort="${COHORT}" --restart-times 0 --cluster-config "${CLUSTRCONF}" --configfile "${CONFIGFILE}" --jobname "{jobid}.{cluster.name}" --drmaa " -S /bin/bash -j {cluster.j} -M {cluster.M} -m {cluster.m} -q {cluster.queu} -l nodes={cluster.nodes}:ppn={cluster.ppn},walltime={cluster.walltime} -l mem={cluster.mem}gb -e ${WORKDIR}/{cluster.stderr} -o ${WORKDIR}/{cluster.stdout}" $EXTRA_OPTS $TARGET
fi

## END ##
