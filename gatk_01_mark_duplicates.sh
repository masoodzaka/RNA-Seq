#!/bin/bash
#$ -cwd -V
#$ -pe smp 10
#$ -l h_vmem=18G
# $ -M louise.pease2@ncl.ac.uk
#$ -m a
#$ -N markdups
#$ -hold_jid_ad star
#$ -o /user/logs/mark_dup.$JOB_ID.$TASK_ID.out
#$ -e /user/logs/mark_dup.$JOB_ID.$TASK_ID.err

module add apps/java/jre-1.7.0_75
module add apps/picard/1.130

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

readarray -t FILENAME < ${scripts}/samples.txt

BAM=$bamOutput/2pass/${FILENAME[$SGE_TASK_ID-1]}.Aligned.out.sam

PICARD_JAVA_ARGS='-XX:-UseLargePages -XX:ParallelGCThreads=2 -Xmx16g -Xms16g -Djava.io.tmpdir='${TMPDIR}

java ${PICARD_JAVA_ARGS} -jar ${PICARD_PATH}/picard.jar AddOrReplaceReadGroups I=${BAM} O=${BAM%.out.sam}_rg_sorted.bam SO=coordinate RGID=${FILENAME[$SGE_TASK_ID-1]} RGLB=library RGPL=illumina RGPU=barcode RGSM=${FILENAME[$SGE_TASK_ID-1]}

java $PICARD_JAVA_ARGS -jar ${PICARD_PATH}/picard.jar MarkDuplicates I=${BAM%.out.sam}_rg_sorted.bam O=${BAM%.out.sam}_dedup.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${bamOutput}/2pass/${FILENAME[$SGE_TASK_ID-1]}_output.metrics MAX_RECORDS_IN_RAM=4000000



#GATK_JAVA_ARGS='-Xmx4096m -Xms4096m -XX:-UseLargePages' 
#java $GATK_JAVA_ARGS -jar $GATK_ROOT/GenomeAnalysisTK.jar

end_time=$(date)
end_secs=$(date +%s)
time_elapsed=$(echo "${end_secs} - ${start_secs}" | bc)
echo
echo "### Ended at: ${end_time} ###"
echo "### Time Elapsed: ${time_elapsed} ###"
