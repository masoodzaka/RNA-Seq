#!/bin/bash -e
#$ -cwd -V
#$ -pe smp 1
#$ -l h_vmem=40G
# $ -M louise.pease@ncl.ac.uk
#$ -m a
#$ -N RenameSam
#$ -hold_jid star
#$ -o /user/logs/RenameSam.$JOB_ID.$TASK_ID.out
#$ -e /user/logs/RenameSam.$JOB_ID.$TASK_ID.err


start_time=$(date)
start_secs=$(date +%s)
echo "### Logging Info ###"
echo "### Job: ${JOB_ID}: ${JOB_NAME} ###"
echo "### Array ID: ${SGE_TASK_ID} ###"
echo "### Job Occupying ${NSLOTS} slots ###"
echo "### Job running on ${HOSTNAME} ###"
echo "### Started at: ${start_time} ###"

source /user/scripts/C_variables_hg19.sh

readarray -t FILENAME < ${scripts}/samples.txt

mv ${twopass}/Aligned.out.sam ${twopass}/${FILENAME[$SGE_TASK_ID-1]}.Aligned.out.sam

end_time=$(date)
end_secs=$(date +%s)
time_elapsed=$(echo "${end_secs} - ${start_secs}" | bc)
echo
echo "### Ended at: ${end_time} ###"
echo "### Time Elapsed: ${time_elapsed} ###"
