#!/bin/bash
#$ -cwd
#$ -V
#$ -l h_vmem=10G # memory default is 6
#$ -pe smp 5 # CPUs default is 1
#$ -l h_rt=48:00:00
#$ -hold_jid star

source /user/scripts/C_variables_hg19.sh

${fusioncatcher} \
-d ${WORK}/ \
-i ${WORK}/${FILENAME[$SGE_TASK_ID-1]}/${FILENAME[$SGE_TASK_ID-1]}_1.fastq.gz,${WORK}/${FILENAME[$SGE_TASK_ID-1]}/${FILENAME[$SGE_TASK_ID-1]}_2.fastq.gz \
-o ${fusioncatcher_out}

#inDir=$1
#outDir=$( basename $1 )

#/sharedlustre/users/ndw31/fusioncatcher/fusioncatcher/bin/fusioncatcher \
#-d /sharedlustre/users/ndw31/fusioncatcher/fusioncatcher/data/current/ \
#-i $inDir \
#-o /sharedlustre/users/ndw31/output_fusioncatcher/$outDir/
