#!/bin/bash -e
#$ -cwd -V
#$ -pe smp 1
#$ -l h_vmem=40G
# $ -M louise.pease2@ncl.ac.uk
#$ -m a
#$ -N exon_utils
#$ -o /user/logs/exonutils.$JOB_ID.$TASK_ID.out
#$ -e /user/logs/exonutils.$JOB_ID.$TASK_ID.err

#module load apps/python27/2.7.8
#module load libs/python/numpy/1.9.1-python27-2.7.8
#module load libs/python/matplotlib/1.3.1-python27
#module load libs/python/pysam/0.7.5
module load apps/miso/0.5.3

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

echo "exon_utils --get-const-exons ${GFF}  --min-exon-size 1000 --output-dir ${miso_out}"
exon_utils --get-const-exons ${GFF}  --min-exon-size 1000 --output-dir ${miso_out}

end_time=$(date)
end_secs=$(date +%s)
time_elapsed=$(echo "${end_secs} - ${start_secs}" | bc)
echo
echo "### Ended at: ${end_time} ###"
echo "### Time Elapsed: ${time_elapsed} ###"
