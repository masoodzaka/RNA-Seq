#!/bin/bash
#$ -N variables
###### This script is set to run using the hg19 genome build
##### leave home/nlp71/WORKING_DATA/NICR_RNAseq/rna_seq/test_data/fastqs/fastqs2 as home/nlp71/WORKING_DATA/NICR_RNAseq/rna_seq/test_data/fastqs/fastqs2 sedsetup script will ensure your home/nlp71/WORKING_DATA/NICR_RNAseq/rna_seq/test_data/fastqs/fastqs2 details are replaced with the ones set on every run
#### this avoids unnecessary editing of scripts for each home/nlp71/WORKING_DATA/NICR_RNAseq/rna_seq/test_data/fastqs/fastqs2, the sedsetup sets the home/nlp71/WORKING_DATA/NICR_RNAseq/rna_seq/test_data/fastqs/fastqs2 details and these are carried forward as a variable
### the modified directory path for the working directory becomes

WORK=/user/

scripts=${WORK}/scripts
#scripts=/home/nlp71/WORKING_DATA/NICR_RNAseq/rna_seq/test_data/fastqs/new_scripts
#scripts=$WORK/'scripts'

#ls $WORK/T* | cut -d "_" -f 1 | uniq -d | tr " " "\n" > ${scripts}/samples.txt
ls $WORK/T* | egrep '*fastq.gz' | cut -d "_" -f 1 | uniq -d | tr " " "\n" > ${scripts}/samples.txt
readarray -t FILENAME < ${scripts}/samples.txt

SAMPLE_NUMBER=$(wc -l ${scripts}/samples.txt | cut -d " " -f 1)
cores=6

gtfFile=/sharedlustre/NICR/human_genomes/hg19/gencode.v19.annotation.gtf
#gtfFile=/sharedlustre/NICR/human_genomes/hg19/gencode.v27lift37.annotation.gtf
gtfFlatFile=/sharedlustre/NICR/human_genomes/hg19/gencode.v19.annotation.gff
GFF=/sharedlustre/NICR/human_genomes/hg19/gencode.v19.annotation.gff
#gtfFlatFile=/sharedlustre/NICR/human_genomes/hg19/gencode.v27lift37.annotation.gff3
houseBed=/sharedlustre/NICR/human_genomes/hg19/hg19.HouseKeepingGenes.bed
refBed=/sharedlustre/NICR/human_genomes/hg19/hg19.refseq.bed.12
##genomeDir=/sharedlustre/NICR/human_genomes/hg19/hg19_Gencode14.overhang75 #### need to get this genome version on sharedlustre
#stargenomeDir=/sharedlustre/NICR/human_genomes/hg38/GRCh38_76bp  #### remember to change this back once tests complete
genomeDirstar=/sharedlustre/NICR/human_genomes/hg19/STAR_hg19_genome
#genomeDirstar=/home/nlp71/WORKING_DATA/NICR_RNAseq/rna_seq/test_data/fastqs
dexSeq=/sharedlustre/NICR/RNAseq_pipelines/DW_RNAseq/dexseq_count.py
#star=/opt/STAR_2.3.0e.Linux_x86_64/STAR
genomeFasta=/opt/databases/GATK_bundle/2.8/hg19/ucsc.hg19.fasta
#picard=/opt/picard-tools-1.118/picard-tools-1.118
picard=/opt/software/bsu/bin/picard-tools-1.118/picard-tools-1.118
#gatk=/opt/GenomeAnalysisTK.jar
gatk=/opt/software/bsu/bin/GenomeAnalysisTK.jar
BUNDLE_DIR="/opt/databases/GATK_bundle/2.8/hg19"
goldStan="$BUNDLE_DIR/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
dbSNP="$BUNDLE_DIR/dbsnp_138.hg19.vcf"
#express=/home/ndw31/bin/express-1.5.1-linux_x86_64/express  #### this is installed on the nicrcompsvr need an installation on fmscluster
#express=/sharedlustre/users/ndw31
#transcriptFA=/home/ndw31/WORKING_DATA/data/ucsc/hg19/Annotation/Genes/Gencode25/gencode.v25lift37.transcripts.fa  #### still need this file
transcriptsFA=/sharedlustre/NICR/human_genomes/hg19/gencode.v25lift37.transcripts.fa 
CHROMSIZES="/sharedlustre/NICR/human_genomes/hg19/hg19.chrom.size"
bowtieInd="/sharedlustre/NICR/human_genomes/hg19"
#EXPECT_READ_LENGTH_AFTER_QC=84
#DEXSEQ_GTF="/sharedlustre/NICR/human_genomes/hg19/gencode.v19.annotation.4DEXSeq.gtf"
Tophat_transcriptome_index="/sharedlustre/NICR/human_genomes/hg19/gencode.v19.annotation"
output_int=$gatk_out
GATK_BUNDLE=/opt/databases/GATK_bundle/2.8/hg19
HG19=${GATK_BUNDLE}/ucsc.hg19.fasta
NSLOTS=5
fusioncatcher=/sharedlustre/users/ndw31/fusioncatcher/fusioncatcher/bin/fusioncatcher


#alignments=$WORK/${FILENAME[$SGE_TASK_ID-1]}/'alignments'
processedData==$WORK/${FILENAME[$SGE_TASK_ID-1]}/'processedData'
fastqcOutput=$WORK/${FILENAME[$SGE_TASK_ID-1]}/processedData/'fastqc'
trimmedFastqs=$WORK/${FILENAME[$SGE_TASK_ID-1]}/processedData/'trimmedfastqs'
htseqOutput=$WORK/${FILENAME[$SGE_TASK_ID-1]}/processedData/'gene_counts'
bamOutput=$WORK/${FILENAME[$SGE_TASK_ID-1]}/processedData/'bam_files'
dexOutput=$WORK/${FILENAME[$SGE_TASK_ID-1]}/processedData/'exon_counts'
cuffQuantOut=$WORK/${FILENAME[$SGE_TASK_ID-1]}/processedData/'cuffquant'
rseqcOutput=$WORK/${FILENAME[$SGE_TASK_ID-1]}/processedData/'rseqc'
vcf=$WORK/${FILENAME[$SGE_TASK_ID-1]}/processedData/'vcf'
twopass=$WORK/${FILENAME[$SGE_TASK_ID-1]}/processedData/bam_files/'2pass'
secPassDir=$WORK/${FILENAME[$SGE_TASK_ID-1]}/processedData/'_2pass'  ### Are 2pass and _2pass different things?  Used intermittently 
#sampleGenomeDir=$WORK/${FILENAME[$SGE_TASK_ID-1]}/processedData/'2pass'
output_gatk=$WORK/${FILENAME[$SGE_TASK_ID-1]}/processedData/'gatk'
fusioncatcher_out=$WORK/${FILENAME[$SGE_TASK_ID-1]}/processedData/'fussion_out'
miso_out=$WORK/${FILENAME[$SGE_TASK_ID-1]}/processedData/'miso_isoforms'
express_out=$WORK/${FILENAME[$SGE_TASK_ID-1]}/processedData/'express_isoforms'
outCuffdiff=$WORK/${FILENAME[$SGE_TASK_ID-1]}/processedData/'Cuffdiff_output'
outCufflinks=$WORK/${FILENAME[$SGE_TASK_ID-1]}/processedData/'Cufflinks_output'
#fussionCatcherInput=$WORK/${FILENAME[$SGE_TASK_ID-1]}/procesedData/'fussioncatcher_input'  ### this si equal to the fastqs so variable is $WORK, folder not needed
#dexSeq=$WORK/${FILENAME[$SGE_TASK_ID-1]}/procesedData/'dexSeq_output'   ### this is the Dexseq python script not a folder, is output folder needed?
gatk_out=$WORK/${FILENAME[$SGE_TASK_ID-1]}/processedData/'gatk_ouput'
BowtieOut=$WORK/${FILENAME[$SGE_TASK_ID-1]}/processedData/'Bowtie_aligned'
HomerOut=$WORK/${FILENAME[$SGE_TASK_ID-1]}/processedData/'Homer_Out'
#starGenomeDir=$WORK/'STAR_genome'
rMATs_out=$WORK/${FILENAME[$SGE_TASK_ID-1]}/processedData/'rMATs_isoforms'
