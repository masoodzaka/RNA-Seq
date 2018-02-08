# prep0_fastqc version 1
# Jannetta's bioinformatics scripts repuprosed for NICR by Louise Pease October 2017
# The purpose of this script is t run fastqc on all your sample files.
# The source file variables.sh, is included with values set for a few
# variables that will be used in the pipelines. Make sure these are set
# before you run the scripts. Of importance for this script is the value
# of WORK.
#
#!/bin/bash -e
#$ -cwd -V
#$ -pe smp 1
#$ -l h_vmem=40G
#$ -q all.q
#$ -m a
#$ -N fastqc
#$ -o /user/logs/fastqc.$JOB_ID.$TASK_ID.out
#$ -e /user/logs/fastqc.$JOB_ID.$TASK_ID.err

module add apps/fastqc/0.11.5
###module add apps/java/jre-1.8.0
##export DISPLAY=:0.0

start_time=$(date)
start_secs=$(date +%s)
echo "### Logging Info ###"
echo "###â| Job: ${JOB_ID}: ${JOB_NAME} ###"
echo "### Array ID: ${SGE_TASK_ID} ###"
echo "###â| Job Occupying ${NSLOTS} slots ###"
echo "### Job running on ${HOSTNAME} ###"
echo "### Started at: ${start_time} ###"
echo

source /user/scripts/C_variables_hg19.sh
source ${scripts}/fastqc_variables.sh


readarray -t FILENAME < ${scripts}/samples.txt

echo "fastqc -o $fastqcOutput "
fastqc $WORK/${FILENAME[$SGE_TASK_ID-1]}/${FILENAME[$SGE_TASK_ID-1]}_1.fastq.gz $WORK/${FILENAME[$SGE_TASK_ID-1]}/${FILENAME[$SGE_TASK_ID-1]}_2.fastq.gz -o $fastqcOutput 

end_time=$(date)
end_secs=$(date +%s)
time_elapsed=$(echo "${end_secs} - ${start_secs}" | bc)
echo
echo "### Ended at: ${end_time} ###"
echo "### Time Elapsed: ${time_elapsed} ###"
