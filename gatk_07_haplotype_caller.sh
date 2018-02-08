#!/bin/bash
#$ -cwd -V
#$ -pe smp 10
#$ -l h_vmem=110G
# $ -M louise.pease2@ncl.ac.uk
#$ -m a
#$ -N haplotypes
#$ -hold_jid print_reads
# $ -o /user/logs/haplotype_caller.$JOB_ID.$TASK_ID.out
# $ -e /user/logs/haplotype_caller.$JOB_ID.$TASK_ID.err

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
echo

source /user/scripts/C_variables_hg19.sh

BAM=${bamOutput}/2pass/${FILENAME[$SGE_TASK_ID-1]}.Aligned_recalibrated.bam
#BAMFILE=Aligned_recalibrated.bam
BAMFILE=${FILENAME[$SGE_TASK_ID-1]}.Aligned_recalibrated.bam

samtools index ${BAM}

# move bamfile and index to localscratch
cp ${BAM}* ${TMPDIR}

GATK_JAVA_ARGS='-Xmx96g -Xms96g -XX:-UseLargePages' 
java ${GATK_JAVA_ARGS} -jar ${GATK_ROOT}/GenomeAnalysisTK.jar \
   -T HaplotypeCaller \
   -nct 10 \
   -R ${HG19} \
   -I ${TMPDIR}/${BAMFILE} \
   --dbsnp ${GATK_BUNDLE}/dbsnp_138.hg19.vcf \
   -dontUseSoftClippedBases \
   -stand_call_conf 20.0 \
   -stand_emit_conf 20.0 \
   -o ${vcf}/${FILENAME[$SGE_TASK_ID-1]}.raw.vcf
                      
end_time=$(date)
end_secs=$(date +%s)
time_elapsed=$(echo "${end_secs} - ${start_secs}" | bc)
echo
echo "### Ended at: ${end_time} ###"
echo "### Time Elapsed: ${time_elapsed} ###"
