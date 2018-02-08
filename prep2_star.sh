#!/bin/bash
#$ -cwd -V
#$ -pe smp 5
#$ -l h_vmem=40G
#$ -m a
#$ -N starprep
#$ -o /user/logs/hg19_74.$JOB_ID.$TASK_ID.out
#$ -e /user/logs/hg19_74.$JOB_ID.$TASK_ID.err

module add apps/STAR/2.5.2b

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

# 1) STAR uses genome index files that must be saved in unique directories. The human genome index was built from the FASTA file hg19.fa as follows:

if [ ! -d ${genomeDirstar} ]; then
  mkdir ${genomeDirstar}
fi

STAR --runMode genomeGenerate --genomeDir ${genomeDirstar} \
  --genomeFastaFiles ${genomeFasta} \
  --sjdbGTFfile ${gtfFile} \
  --runThreadN ${cores} \
  --sjdbOverhang 74 

echo
echo "### Ended at: ${end_time} ###"
echo "### Time Elapsed: ${time_elapsed} ###"
