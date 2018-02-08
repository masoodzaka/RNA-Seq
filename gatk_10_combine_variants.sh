#!/bin/bash
#$ -cwd -V
#$ -pe smp 1
#$ -l h_vmem=50G
# $ -M louise.pease2@ncl.ac.uk
#$ -m a
#$ -N combine_var
#$ -hold_jid select_variants
# $ -o /user/logs/combine_variants.$JOB_ID.$TASK_ID.out
# $ -e /user/logs/combine_variants.$JOB_ID.$TASK_ID.err

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
   -T CombineVariants \
   -R ${HG19} \
   --variant:Sample1 ${vcf}/${FILENAME[$SGE_TASK_ID-1]}.PASS.vcf \
   --variant:Sample2 ${vcf}/${FILENAME[$SGE_TASK_ID-1]}.PASS.vcf \
#   --variant:Sample3 ${vcf}/${FILENAME[$SGE_TASK_ID-1]}.PASS.vcf \
#   --variant:Sample4 ${vcf}/${FILENAME[$SGE_TASK_ID-1]}.PASS.vcf  \
#   --variant:Sample5 ${vcf}/${FILENAME[$SGE_TASK_ID-1]}.PASS.vcf  \
#   --variant:Sample6 ${vcf}/${FILENAME[$SGE_TASK_ID-1]}.PASS.vcf  \
#   --variant:Sample7 ${vcf}/${FILENAME[$SGE_TASK_ID-1]}.PASS.vcf  \
#   --variant:Sample8 ${vcf}/${FILENAME[$SGE_TASK_ID-1]}.PASS.vcf  \
#   --variant:Sample9 ${vcf}/${FILENAME[$SGE_TASK_ID-1]}.PASS.vcf  \
#   --variant:Sample10 ${vcf}/${FILENAME[$SGE_TASK_ID-1]}.PASS.vcf  \
   -genotypeMergeOptions PRIORITIZE \
   -priority Sample1,Sample2 \
#   -priority Sample1,Sample2,Sample3,Sample4,Sample5,Sample6,Sample7,Sample8,Sample9,Sample10 \
   -o ${WORK}/all_variants.vcf 
                      
end_time=$(date)
end_secs=$(date +%s)
time_elapsed=$(echo "${end_secs} - ${start_secs}" | bc)
echo
echo "### Ended at: ${end_time} ###"
echo "### Time Elapsed: ${time_elapsed} ###"
