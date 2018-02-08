#!/bin/bash -e
#$ -cwd -V
#$ -pe smp 1-5
#$ -l h_vmem=100G
# $ -M name@domain.com
#$ -m a
#$ -N de_star
#$ -hold_jid genomeprep
#$ -o /user/logs/de_star.$JOB_ID.$TASK_ID.out
#$ -e /user/logs/de_star.$JOB_ID.$TASK_ID.err

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

#source /user/scripts/C_variables_hg19.sh
source ${scripts}/C_variables_hg19.sh

readarray -t FILENAME < ${scripts}/samples.txt
### Can change setings here to use the trimmed fastq file 
#FASTQ1=${trimmedFastqs}/${FILENAME[$SGE_TASK_ID-1]}*_R1*.fastq.gz.trimmed.fastq.gz
#FASTQ2=${trimmedFastqs}/${FILENAME[$SGE_TASK_ID-1]}*_R2*.fastq.gz.trimmed.fastq.gz
FASTQ1=${trimmedFastqs}/${FILENAME[$SGE_TASK_ID-1]}_1.fastq.gz.trimmed.fastq.gz 
FASTQ2=${trimmedFastqs}/${FILENAME[$SGE_TASK_ID-1]}_2.fastq.gz.trimmed.fastq.gz 


runDir=${TMPDIR}/${FILENAME[$SGE_TASK_ID-1]}
mkdir $runDir
cd $runDir

echo "STAR --genomeDir ${GENDIR}  --readFilesIn ${FASTQ1} ${FASTQ2}  --readFilesCommand gunzip -c  --runThreadN ${NSLOTS}  --outSAMtype BAM SortedByCoordinate"
STAR \
  --genomeDir ${GENDIR} \
  --readFilesIn ${FASTQ1} ${FASTQ2} \
  --readFilesCommand gunzip -c \
  --runThreadN ${NSLOTS} \
  --outSAMtype BAM Unsorted \
  --twopassMode Basic

mv $runDir ${bamOutput}

end_time=$(date)
end_secs=$(date +%s)
time_elapsed=$(echo "${end_secs} - ${start_secs}" | bc)
echo
echo "### Ended at: ${end_time} ###"
echo "### Time Elapsed: ${time_elapsed} ###"


