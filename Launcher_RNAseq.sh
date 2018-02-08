#!/bin/bash
#$ -cwd
#$ -V
#$ -l h_vmem=100G
#$ -pe smp 1-5 # CPUs default is 1
#$ -l h_rt=80:00:00
#$ -N RNAseq_launch
#$ -q bigmem.q,all.q
#$ -o /user/logs/RNAseq_launch.$JOB_ID.$TASK_ID.out
#$ -e /user/logs/RNAseq_launch.$JOB_ID.$TASK_ID.err


###module add apps/samtools/1.2
#module add apps/samtools/1.3.1
###module add apps/STAR/2.4.0j
module add apps/STAR/2.5.2b 
module add apps/picard/1.130 
#module add apps/gatk/3.2-protected 
module add apps/fastqc/0.11.5 
module add apps/python27/2.7.8 
module add libs/python/numpy/1.9.1-python27-2.7.8 
###module add libs/python/numpy/1.8.0-python27
module add libs/python/pysam/0.10.0 
module add apps/miso/0.5.3
module add apps/bedtools/2.25
module add libs/python/matplotlib/1.3.1-python27
module add apps/HTSeq/0.6.1

##/home/nlp71/WORKING_DATA/NICR_RNAseq/rna_seq/test_data/fastqs/new_scripts

scripts=/user/scripts

source ${scripts}/C_variables_hg19.sh

readarray -t FILENAME < ${scripts}/samples.txt

#FASTQ1=${WORK}\${FILENAME[$SGE_TASK_ID-1]}_1.fastq.gz
#FASTQ2=${WORK}\${FILENAME[$SGE_TASK_ID-1]}_2.fastq.gz
#echo ${FILENAME[$SGE_TASK_ID-1]}
#echo $readFile1
#echo $readFile2
#echo $FASTQ1
#echo $FASTQ2
################## SETUP SCRIPTS TO GENERATE AND SET DIRECTORIES AND VARIABLES FOR FILES AND SCRIPTS ###################### 
#### change the filepath in the sedsetup script before you run it
#qsub A_copy_scripts.sh

qsub B_sedsetup

#### Do not run the variables scripts they will be sourced in subsequent scripts
####qsub -t 1:${SAMPLE_NUMBER} C_variables_hg38.sh -hold_jid setup
####qsub -t 1:${SAMPLE_NUMBER} C_variables_hg19.sh -hold_jid setup
 
qsub -t 1:${SAMPLE_NUMBER} D_prepare_folders.sh -hold_jid variables

######################## PREPAROTORY SCRIPTS FOR TRIMMING, GENOME CREATION AND ALIGNMENT ################################### 
#### only needed if trimming fails, check phred scores are in the right format
####qsub prep_determine_quality_scores.sh

qsub -t 1:${SAMPLE_NUMBER} prep0_fastqc.sh -hold_jid prep_folders

qsub -t 1:${SAMPLE_NUMBER} prep1_trim.sh -hold_jid fastqc

#### only run once to generate the star genome, not required as star genome is on sharedlustre and set as a variable
####qsub prep2_star.sh -hold_jid prep_folders

#### run alignments with the star genome, generate the 2pass directory and save sam/bam files to the 2pass directory
qsub -t 1:${SAMPLE_NUMBER} gatk_00_run_star.sh -hold_jid starprep

########################################### GATK VARIANT CALLING SCRIPTS ###################################################
qsub -t 1:${SAMPLE_NUMBER} gatk_rename_sam.sh 

qsub -t 1:${SAMPLE_NUMBER} prep_rename_sam.sh

qsub -t 1:${SAMPLE_NUMBER} sam_name.sh

qsub -t 1:${SAMPLE_NUMBER} gatk_01_mark_duplicates.sh 

qsub -t 1:${SAMPLE_NUMBER} gatk_02_split_n_trim.sh -hold_jid markdups

qsub -t 1:${SAMPLE_NUMBER} gatk_03_realigner_target_creator.sh

qsub -t 1:${SAMPLE_NUMBER} gatk_04_indel_realigner.sh

qsub -t 1:${SAMPLE_NUMBER} gatk_05_base_recalibrator.sh

qsub -t 1:${SAMPLE_NUMBER} gatk_06_print_reads.sh   ##### check error logs

qsub -t 1:${SAMPLE_NUMBER} gatk_07_haplotype_caller.sh ##### check error logs

qsub -t 1:${SAMPLE_NUMBER} gatk_08_variant_filtering.sh

qsub -t 1:${SAMPLE_NUMBER} gatk_09_select_variants.sh

qsub -t 1:${SAMPLE_NUMBER} gatk_10_combine_variants.sh

########################### DIFFERENTIAL EXPRESSION AND EXON USAGE SCRIPTS USING DEXSEQ AND HTSEQ #########################

#### first script not needed as done above by #qsub -t 1:${SAMPLE_NUMBER} gatk_00_run_star.sh -hold_jid starprep
#### only included her if just DE is to be run 
####qsub -t 1:${SAMPLE_NUMBER} de00_star.sh

qsub -t 1:${SAMPLE_NUMBER} de01_mvbam.sh

qsub -t 1:${SAMPLE_NUMBER} de02_sortbamname.sh

qsub -t 1:${SAMPLE_NUMBER} de03_sortbamcoord.sh

qsub -t 1:${SAMPLE_NUMBER} de04_indexbam.sh

qsub -t 1:${SAMPLE_NUMBER} de05_htseqcoord.sh

qsub -t 1:${SAMPLE_NUMBER} de06_htseqname.sh

########################################################
#### Currently unused as depreciated, modular runs by array above
####qsub -t 1:${SAMPLE_NUMBER} sortsamname_func.sh

####qsub -t 1:${SAMPLE_NUMBER} process_function.sh

################################## FUSION GENE IDENTIFICATION USING FUSIONCATCHER ##########################################

qsub -t 1:${SAMPLE_NUMBER} fu00_fusioncatcher.sh

################################## ISOFORM IDENTIFICATION USING EXPRESS, MISO, HOMER, rMATS ################################

qsub -t 1:${SAMPLE_NUMBER} iso00_express_gsub.sh 

qsub -t 1:${SAMPLE_NUMBER} iso01_homer.sh

qsub -t 1:${SAMPLE_NUMBER} iso02_rMATs.sh

qsub -t 1:${SAMPLE_NUMBER} iso03_miso.sh

qsub -t 1:${SAMPLE_NUMBER} miso01_exonutils.sh

qsub -t 1:${SAMPLE_NUMBER} miso02_peutils.sh

qsub -t 1:${SAMPLE_NUMBER} miso03_indexgff.sh

qsub -t 1:${SAMPLE_NUMBER} miso04_miso.sh

#### Not required as this is done in the prep above 
#####qsub Generate_genome.sh -hold_jid D_prep_folders


