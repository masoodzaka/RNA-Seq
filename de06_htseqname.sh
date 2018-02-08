#!/bin/bash -e
#$ -cwd -V
#$ -pe smp 1
#$ -l h_vmem=40G
# $ -M louise.pease2@ncl.ac.uk
#$ -m a
#$ -N de_htseq
#$ -hold_jid de_indexbam
#$ -o /user/logs/de_htseq.$JOB_ID.$TASK_ID.out
#$ -e /user/logs/de_htseq.$JOB_ID.$TASK_ID.err

module load apps/python27/2.7.8
module load libs/python/numpy/1.9.1-python27-2.7.8
module load libs/python/matplotlib/1.3.1-python27
module load libs/python/pysam/0.7.5
module load apps/HTSeq/0.6.1

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

#echo "htseq-count -f bam -r name -s reverse -a 4 -i gene_id -m union ${bamOutput}/${FILENAME[$SGE_TASK_ID-1]}.sorted.bam $GTF > ${bamOutput}/${FILENAME[$SGE_TASK_ID-1]}.counts"
#htseq-count -f bam -r name -s reverse -a 4 -i gene_id -m union \
#  "${bamOutput}/${FILENAME[$SGE_TASK_ID-1]}.sorted.bam" \
#  $GTF \
#  > ${htseqOutput}/${FILENAME[$SGE_TASK_ID-1]}.counts

echo "htseq-count -f bam -r name -s reverse -a 4 -i gene_id -m union ${twopass}/${FILENAME[$SGE_TASK_ID-1]}.sorted.bam $GTF > ${twopass}/${FILENAME[$SGE_TASK_ID-1]}.counts"
htseq-count -f bam -r name -s reverse -a 4 -i gene_id -m union \
  "${twopass}/${FILENAME[$SGE_TASK_ID-1]}.sorted.bam" \
  $GTF \
  > ${htseqOutput}/${FILENAME[$SGE_TASK_ID-1]}.counts


end_time=$(date)
end_secs=$(date +%s)
time_elapsed=$(echo "${end_secs} - ${start_secs}" | bc)
echo
echo "### Ended at: ${end_time} ###"
echo "### Time Elapsed: ${time_elapsed} ###"

