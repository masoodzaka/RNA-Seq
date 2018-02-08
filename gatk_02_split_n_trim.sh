#!/bin/bash
#$ -cwd -V
#$ -pe smp 10
#$ -l h_vmem=110G
# $ -M louise.pease2@ncl.ac.uk
#$ -m a
#$ -N splittrim
#$ -hold_jid_ad markdups
#$ -o /user/logs/split_n_trim.$JOB_ID.$TASK_ID.out
#$ -e /user/logs/split_n_trim.$JOB_ID.$TASK_ID.err

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

readarray -t FILENAME < ${scripts}/samples.txt



BAM=$bamOutput/2pass/${FILENAME[$SGE_TASK_ID-1]}.Aligned_dedup.bam
BAMFILE=${FILENAME[$SGE_TASK_ID-1]}.Aligned_dedup.bam

# move bamfile and index to localscratch
cp ${BAM}* ${TMPDIR}
#cp ${BAM%.bam}.bai $TMPDIR

GATK_JAVA_ARGS='-Xmx96g -Xms96g -XX:-UseLargePages' 

java $GATK_JAVA_ARGS -jar $GATK_ROOT/GenomeAnalysisTK.jar \
    -T SplitNCigarReads -R ${HG19} \
    -I ${TMPDIR}/${BAMFILE} \
    -o ${TMPDIR}/${BAMFILE%_dedup.bam}_split.bam \
    -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 \
    -U ALLOW_N_CIGAR_READS

mv ${TMPDIR}/${BAMFILE%_dedup.bam}_split.bam ${bamOutput}/2pass

end_time=$(date)
end_secs=$(date +%s)
time_elapsed=$(echo "${end_secs} - ${start_secs}" | bc)
echo
echo "### Ended at: ${end_time} ###"
echo "### Time Elapsed: ${time_elapsed} ###"
