#!/bin/bash
#$ -cwd -V
#$ -pe smp 10
#$ -l h_vmem=110G
# $ -M louise.pease2@ncl.ac.uk
#$ -m a
#$ -N gatk_realign
#$ -hold_jid_ad splittrim
#$ -o /user/logs/realigner_target_creator.$JOB_ID.$TASK_ID.out
#$ -e /user/logs/realigner_target_creator.$JOB_ID.$TASK_ID.err

module add apps/java/jre-1.7.0_75
module add apps/gatk/3.3-protected
module add apps/samtools/1.3.1

start_time=$(date)
start_secs=$(date +%s)
echo "### Logging Info ###"
echo "### Job: ${JOB_ID}: ${JOB_NAME} ###"
echo "### Array ID: ${SGE_TASK_ID} ###"
echo "### Job Occupying ${NSLOTS} slots ###"
echo "### Job running on ${HOSTNAME} ###"
echo "### Started at: ${start_time} ###"

source /user/scripts/C_variables_hg19.sh

BAM=$bamOutput/2pass/${FILENAME[$SGE_TASK_ID-1]}.Aligned_split.bam
samtools index ${BAM}

BAMFILE=${FILENAME[$SGE_TASK_ID-1]}.Aligned_split.bam
#BAMFILE=Aligned_split.bam
# move bamfile and index to localscratch
cp ${BAM}* ${TMPDIR}
#cp ${BAM%.bam}.bai $TMPDIR
#cp ${BAM}.bai $TMPDIR


GATK_JAVA_ARGS='-Xmx96g -Xms96g -XX:-UseLargePages' 
java ${GATK_JAVA_ARGS} -jar ${GATK_ROOT}/GenomeAnalysisTK.jar \
   -T RealignerTargetCreator \
   -nt ${NSLOTS} \
   -R ${HG19} \
   -I ${TMPDIR}/${BAMFILE} \
   -o ${bamOutput}/2pass/${FILENAME[$SGE_TASK_ID-1]}.intervals \
   --known ${GATK_BUNDLE}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
   --known ${GATK_BUNDLE}/1000G_phase1.indels.hg19.sites.vcf

end_time=$(date)
end_secs=$(date +%s)
time_elapsed=$(echo "${end_secs} - ${start_secs}" | bc)
echo
echo "### Ended at: ${end_time} ###"
echo "### Time Elapsed: ${time_elapsed} ###"
