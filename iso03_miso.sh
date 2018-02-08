#!/bin/bash
#$ -cwd
#$ -V
#$ -l h_vmem=10G # memory default is 6
#$ -pe smp 5 # CPUs default is 1
#$ -l h_rt=48:00:00
#$ -q bigmem.q,all.q
# $ -hold_jid genomeprep
# $ -o /user/logs/miso.$JOB_ID.$TASK_ID.out
# $ -e /user/logs/miso.$JOB_ID.$TASK_ID.err

module add apps/miso/0.5.3
module add apps/bedtools/2.25

#source /user/scripts/variables.sh 
scripts=/user/scripts
source ${scripts}/C_variables_hg19.sh

#gtfFile='/home/ndw31/WORKING_DATA/data/ucsc/hg19/Annotation/Genes/Gencode25/gencode.v25lift37.annotation.gtf'

#gffFile='/home/ndw31/WORKING_DATA/data/ucsc/hg19/Annotation/Genes/Gencode25/gencode.v25lift37.annotation.gff3'

#gffIndex='/home/ndw31/WORKING_DATA/data/ucsc/hg19/Annotation/Genes/Gencode25/gencode.v25lift37.annotation.gff3.index'

#index_gff --index $gffFile $gffIndex
#exon_utils --get-const-exons $gffFile --min-exon-size 1000 --output-dir /home/ndw31/WORKING_DATA/data/ucsc/hg19/Annotation/Genes/Gencode25/exons/

bamFile=$1
fileName=$( basename $bamFile)
dirName=$( dirname $bamFile)
header=$( head -1 $dirName/insert-dist/$fileName.insert_len )
pe_utils --compute-insert-len $bamFile /home/ndw31/WORKING_DATA/data/ucsc/hg19/Annotation/Genes/Gencode25/exons/gencode.v25lift37.annotation.min_1000.const_exons.gff --output-dir $dirName/insert-dist/

IFS=', ' read -r -a array <<< "$header"
echo ${array[0]}
mean=$(echo ${array[0]} | sed 's/#mean=//g' )
sd=$(echo ${array[1]} | sed 's/sdev=//g' )
miso --run $gffIndex $bamFile --output-dir $dirName/miso --read-len 100 --paired-end $mean $sd --settings /home/ndw31/WORKING_DATA/MISO/miso_settings.txt
summarize_miso --summarize-samples $dirName/miso $dirName/miso/summary_output/

#compare_miso --compare-samples control/ knockdown/ comparisons/
