#!/bin/bash 
#$ -N stargenome 
#$ -q bigmem.q,all.q 
#$ -l h_vmem=30G # memory default is 6 
#$ -l h_rt=40:00:00
#$ -o /user/logs/stargen.$JOB_ID.$TASK_ID.out
#$ -e /user/logs/stargen.$JOB_ID.$TASK_ID.err

module add apps/STAR/2.5.2b

source /user/scripts/C_variables_hg19.sh

#starGenomeDir=/sharedlustre/NICR/human_genomes/hg19/hg19_Gencode14.overhang75
#starGenomeDir=/home/nlp71/WORKING_DATA/NICR_RNAseq/rna_seq/test_data/fastqs/new_scripts
#genomeFasta=/opt/databases/GATK_bundle/2.8/hg19/ucsc.hg19.fasta


STAR --runMode genomeGenerate --genomeDir $starGenomeDir --genomeFastaFiles $genomeFasta --runThreadN 40 --genomeLoad LoadAndRemove
#STAR --runMode genomeGenerate --genomeDir $GenomeDirstar --genomeFastaFiles $genomeFasta --runThreadN 40 --genomeLoad LoadAndRemove

