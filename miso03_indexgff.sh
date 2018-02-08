#!/bin/bash -e
#$ -cwd -V
#$ -pe smp 1
#$ -l h_vmem=40G
# $ -M name@domain.com
#$ -m a
#$ -N indexgff
#$ -hold_jid pe_utils
#$ -o /user/logs/indexgff.$JOB_ID.$TASK_ID.out
#$ -e /user/logs/indexgff.$JOB_ID.$TASK_ID.err

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

index_gff --index ${GFF} ${miso_out}.indexed/

end_time=$(date)
end_secs=$(date +%s)
time_elapsed=$(echo "${end_secs} - ${start_secs}" | bc)
echo
echo "### Ended at: ${end_time} ###"
echo "### Time Elapsed: ${time_elapsed} ###"
