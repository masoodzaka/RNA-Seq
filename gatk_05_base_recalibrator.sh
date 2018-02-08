#!/bin/bash
#$ -cwd -V
#$ -pe smp 10
#$ -l h_vmem=110G
# $ -M louise.pease2@ncl.ac.uk
#$ -m a
#$ -N recalibrate
#$ -hold_jid_ad indel_realign
#$ -o /user/logs/base_recalibrator.$JOB_ID.$TASK_ID.out
#$ -e /user/logs/base_recalibrator.$JOB_ID.$TASK_ID.err

module add apps/java/jre-1.7.0_75
module add apps/gatk/3.3-protected
module add apps/samtools

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

BAM=${bamOutput}/2pass/${FILENAME[$SGE_TASK_ID-1]}.Aligned_realigned.bam
#BAMFILE=Aligned_realigned.bam
BAMFILE=${FILENAME[$SGE_TASK_ID-1]}.Aligned_realigned.bam
#BAMFILE=$(basename ${BAM})

samtools index ${BAM}

# move bamfile and index to localscratch
cp ${BAM}* ${TMPDIR}

GATK_JAVA_ARGS='-Xmx96g -Xms96g -XX:-UseLargePages' 
java ${GATK_JAVA_ARGS} -jar ${GATK_ROOT}/GenomeAnalysisTK.jar \
   -T BaseRecalibrator \
   -nct 10 \
   -R ${HG19} \
   -I ${TMPDIR}/${BAMFILE} \
   -knownSites ${GATK_BUNDLE}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
   -knownSites ${GATK_BUNDLE}/1000G_phase1.indels.hg19.sites.vcf \
   -knownSites ${GATK_BUNDLE}/dbsnp_138.hg19.vcf \
   -o ${bamOutput}/2pass/${FILENAME[$SGE_TASK_ID-1]}.grp
                      
end_time=$(date)
end_secs=$(date +%s)
time_elapsed=$(echo "${end_secs} - ${start_secs}" | bc)
echo
echo "### Ended at: ${end_time} ###"
echo "### Time Elapsed: ${time_elapsed} ###"
