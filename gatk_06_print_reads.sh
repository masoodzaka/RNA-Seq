#!/bin/bash
#$ -cwd -V
#$ -pe smp 10
#$ -l h_vmem=110G
# $ -M louise.pease2@ncl.ac.uk
#$ -m a
#$ -N print_reads
#$ -hold_jid_ad recalibrate 
#$ -o /user/logs/print_reads.$JOB_ID.$TASK_ID.out
#$ -e /user/logs/print_reads.$JOB_ID.$TASK_ID.err

module add apps/java/jre-1.7.0_75
module add apps/gatk/3.3-protected

start_time=$(date)
start_secs=$(date +%s)
echo "### Logging Info ###"
echo "### Job: ${JOB_ID}: ${JOB_NAME} ###"
echo "### Array ID: ${SGE_TASK_ID} ###"
echo "### Job Occupying ${NSLOTS} slots ###"
echo "### Job running on ${HOSTNAME} ###"
echo "### Started at: ${start_time} ###"
echo

source /user/scripts/C_variables_hg19.sh

#BAM=${bamOutput}/2pass/Aligned_realigned.bam
BAM=${bamOutput}/2pass/${FILENAME[$SGE_TASK_ID-1]}.Aligned_realigned.bam
#BAMFILE=Aligned_realigned.bam
BAMFILE=${FILENAME[$SGE_TASK_ID-1]}.Aligned_realigned.bam

# move bamfile and index to localscratch
cp ${BAM}* ${TMPDIR}

GATK_JAVA_ARGS='-Xmx96g -Xms96g -XX:-UseLargePages' 
java ${GATK_JAVA_ARGS} -jar ${GATK_ROOT}/GenomeAnalysisTK.jar \
   -T PrintReads \
   -nct 10 \
   -R ${HG19} \
   -I ${TMPDIR}/${BAMFILE} \
   -BQSR ${bamOutput}/2pass/${FILENAME[$SGE_TASK_ID-1]}.grp \
   -o ${bamOutput}/2pass/${FILENAME[$SGE_TASK_ID-1]}.Aligned_recalibrated.bam
                      
end_time=$(date)
end_secs=$(date +%s)
time_elapsed=$(echo "${end_secs} - ${start_secs}" | bc)
echo
echo "### Ended at: ${end_time} ###"
echo "### Time Elapsed: ${time_elapsed} ###"
