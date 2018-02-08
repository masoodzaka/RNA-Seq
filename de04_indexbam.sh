
#!/bin/bash -e
#$ -cwd -V
#$ -pe smp 1
#$ -l h_vmem=40G
#$ -m a
#$ -N de_indexbam
#$ -hold_jid de_sortbam
#$ -o /user/logs/de_indexbam.$JOB_ID.$TASK_ID.out
#$ -e /user/logs/de_indexbam.$JOB_ID.$TASK_ID.err

module load apps/samtools/1.2

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
#source ${scripts}/C_variables_hg19.sh

readarray -t FILENAME < ${scripts}/samples.txt

#echo samtools index ${bamOuput}/${FILENAME[$SGE_TASK_ID-1]}.sorted.coord.bam
#samtools index ${bamOutput}/${FILENAME[$SGE_TASK_ID-1]}.sorted.coord.bam

echo samtools index ${twopass}/${FILENAME[$SGE_TASK_ID-1]}.sorted.coord.bam
samtools index ${twopass}/${FILENAME[$SGE_TASK_ID-1]}.sorted.coord.bam


end_time=$(date)
end_secs=$(date +%s)
time_elapsed=$(echo "${end_secs} - ${start_secs}" | bc)
echo
echo "### Ended at: ${end_time} ###"
echo "### Time Elapsed: ${time_elapsed} ###"

