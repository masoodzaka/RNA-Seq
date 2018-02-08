#!/bin/bash
#$ -cwd -V
#$ -pe smp 1
#$ -l h_vmem=50G
# $ -M louise.pease2@ncl.ac.uk 
#$ -m a
#$ -N variant_filt
#$ -hold_jid_ad haplotypes
#$ -o /user/logs/variant_filtration.$JOB_ID.$TASK_ID.out
#$ -e /user/logs/variant_filtration.$JOB_ID.$TASK_ID.err

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

GATK_JAVA_ARGS='-Xmx42g -Xms42g -XX:-UseLargePages' 
java ${GATK_JAVA_ARGS} -jar ${GATK_ROOT}/GenomeAnalysisTK.jar \
   -T VariantFiltration \
   -R ${HG19} \
   -V ${vcf}/${FILENAME[$SGE_TASK_ID-1]}.raw.vcf \
   -window 35 \
   -cluster 3 \
   --filterName "RNASEQ_Variant_Filter" --filterExpression "FS > 30.0 || QD < 2.0" \
   -o ${vcf}/${FILENAME[$SGE_TASK_ID-1]}.filtered.vcf 
                      
end_time=$(date)
end_secs=$(date +%s)
time_elapsed=$(echo "${end_secs} - ${start_secs}" | bc)
echo
echo "### Ended at: ${end_time} ###"
echo "### Time Elapsed: ${time_elapsed} ###"
