#!/bin/bash
#$ -cwd -V
#$ -pe smp 1
#$ -l h_vmem=50G
# $ -M louise.pease2@ncl.ac.uk
#$ -m a
#$ -N select_variants
#$ -hold_jid_ad variant_filt 
#$ -o /user/logs/select_variants.$JOB_ID.$TASK_ID.out
#$ -e /user/logs/select_variants.$JOB_ID.$TASK_ID.err

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
   -T SelectVariants \
   -R ${HG19} \
   --variant ${vcf}/${FILENAME[$SGE_TASK_ID-1]}.filtered.vcf \
   -select "vc.isNotFiltered()" \
   -o ${vcf}/${FILENAME[$SGE_TASK_ID-1]}.PASS.vcf 
                      
end_time=$(date)
end_secs=$(date +%s)
time_elapsed=$(echo "${end_secs} - ${start_secs}" | bc)
echo
echo "### Ended at: ${end_time} ###"
echo "### Time Elapsed: ${time_elapsed} ###"
