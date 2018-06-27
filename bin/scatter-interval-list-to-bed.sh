#!\bin\bash
cd /projects/verhaak-lab/verhaak_ref/gatk-legacy-bundles/b37/scattered_wgs_intervals/scatter-5
for f in `find . -type f -name "scattered.interval_list"`;
do
    cat $f | grep -vE "^@" | awk 'OFS="\t" {print $1, $2-1, $3, $5, 0, $4}' > ${f%.*}.bed;
done