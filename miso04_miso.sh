#!/bin/bash -e
#$ -cwd -V
#$ -pe smp 1
#$ -l h_vmem=40G
# $ -M name@domain.com
#$ -m a
#$ -N miso
#$ -hold_jid indexgff
#$ -o /user/logs/miso.$JOB_ID.$TASK_ID.out
#$ -e /user/logs/miso.$JOB_ID.$TASK_ID.err

#module load apps/python27/2.7.8
#module load libs/python/numpy/1.9.1-python27-2.7.8
#module load libs/python/matplotlib/1.3.1-python27
#module load libs/python/pysam/0.7.5

module load apps/miso/0.5.3
module load apps/bedtools/2.25

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

INSERT_LEN=`grep mean ${miso_out}.insert-dist/${FILENAME[$SGE_TASK_ID-1]}.sorted.coord.bam.insert_len |cut -f 1 --delim=","|cut -f2 --delim="="`
SDEV=`grep sdev ${miso_out}.insert-dist/${FILENAME[$SGE_TASK_ID-1]}.sorted.coord.bam.insert_len |cut -f 2 --delim=","|cut -f2 --delim="="`

echo miso --run ${miso_out}.indexed/ ${twopass}/${FILENAME[$SGE_TASK_ID-1]}.sorted.coord.bam --output-dir ${miso_out}/${FILENAME[$SGE_TASK_ID-1]}_output/ --read-len 75 --paired-end ${INSERT_LEN} ${SDEV} --settings-filename=${scripts}/miso_settings.txt
miso \
  --run ${miso_out}.indexed/ ${twopass}/${FILENAME[$SGE_TASK_ID-1]}.sorted.coord.bam \
  --output-dir ${miso_out}/${FILENAME[$SGE_TASK_ID-1]}_output/ \
  --read-len 75 \
  --paired-end ${INSERT_LEN} ${SDEV} \
  --settings-filename=${scripts}/miso_settings.txt 

end_time=$(date)
end_secs=$(date +%s)
time_elapsed=$(echo "${end_secs} - ${start_secs}" | bc)
echo
echo "### Ended at: ${end_time} ###"
echo "### Time Elapsed: ${time_elapsed} ###"
