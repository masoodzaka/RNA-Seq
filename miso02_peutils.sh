#!/bin/bash -e
# $ -cwd -V
#$ -pe smp 1
#$ -l h_vmem=40G
# $ -M louise.pease2@ncl.ac.uk
#$ -m a
#$ -N pe_utils
#$ -hold_jid exon_utils
#$ -o /user/logs/peutils.$JOB_ID.$TASK_ID.out
#$ -e /user/logs/peutils.$JOB_ID.$TASK_ID.err

#module load apps/python27/2.7.8
#module load libs/python/numpy/1.9.1-python27-2.7.8
#module load libs/python/matplotlib/1.3.1-python27
#module load libs/python/pysam/0.7.5
module load apps/miso/0.5.3
module load apps/bedtools

start_time=$(date)
start_secs=$(date +%s)
echo "### Logging Info ###"
echo "### Job: ${JOB_ID}: ${JOB_NAME} ###"
echo "### Array ID: ${SGE_TASK_ID} ###"
echo "### Job Occupying ${NSLOTS} slots ###"
echo "### Job running on ${HOSTNAME} ###"
echo "### Started at: ${start_time} ###"
echo

source /user/scripts/C_variables_hg19.sh

readarray -t FILENAME < ${scripts}/samples.txt

echo pe_utils --compute-insert-len ${twopass}/${FILENAME[$SGE_TASK_ID-1]}.sorted.coord.bam ${miso_out}_ensembl.min_1000.const_exons.gff --output-dir ${miso_out}.insert-dist/
pe_utils --compute-insert-len ${twopass}/${FILENAME[$SGE_TASK_ID-1]}.sorted.coord.bam ${miso_out}_ensembl.min_1000.const_exons.gff --output-dir ${miso_out}.insert-dist/

end_time=$(date)
end_secs=$(date +%s)
time_elapsed=$(echo "${end_secs} - ${start_secs}" | bc)
echo
echo "### Ended at: ${end_time} ###"
echo "### Time Elapsed: ${time_elapsed} ###"
