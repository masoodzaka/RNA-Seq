#!/bin/bash -e
#$ -cwd -V
#$ -pe smp 1
#$ -l h_vmem=40G
# $ -M louise.pease2@newcastle.ac.uk
#$ -m a
#$ -N de_sortbamname
#$ -hold_jid de_mvbam
#$ -o /user/logs/de_sortbam.$JOB_ID.$TASK_ID.out
#$ -e /user/logs/de_sortbam.$JOB_ID.$TASK_ID.err

module load apps/samtools/1.3.1

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

echo samtools sort -o ${twopass}/${FILENAME[$SGE_TASK_ID-1]}.n.sorted.bam -O bam -n ${twopass}/${FILENAME[$SGE_TASK_ID-1]}.Aligned.out.sam -T ${TMPDIR}/samsorttmp
samtools sort -o ${twopass}/${FILENAME[$SGE_TASK_ID-1]}.n.sorted.bam -O bam -n ${twopass}/${FILENAME[$SGE_TASK_ID-1]}.Aligned.out.sam -T ${TMPDIR}/samsorttmp

end_time=$(date)
end_secs=$(date +%s)
time_elapsed=$(echo "${end_secs} - ${start_secs}" | bc)
echo
echo "### Ended at: ${end_time} ###"
echo "### Time Elapsed: ${time_elapsed} ###"
