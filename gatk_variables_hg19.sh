#!/bin/bash

###### This script is set to run using the hg19 genome build
##### leave home/nlp71/WORKING_DATA/NICR_RNAseq/rna_seq/test_data/fastqs/fastqs2 as home/nlp71/WORKING_DATA/NICR_RNAseq/rna_seq/test_data/fastqs/fastqs2 sedsetup script will ensure your home/nlp71/WORKING_DATA/NICR_RNAseq/rna_seq/test_data/fastqs/fastqs2 details are replaced with the ones set on every run
#### this avoids unnecessary editing of scripts for each home/nlp71/WORKING_DATA/NICR_RNAseq/rna_seq/test_data/fastqs/fastqs2, the sedsetup sets the home/nlp71/WORKING_DATA/NICR_RNAseq/rna_seq/test_data/fastqs/fastqs2 details and these are carried forward as a variable
### the modified directory path for the working directory becomes

#WORK=/home/nlp71/WORKING_DATA/NICR_RNAseq/rna_seq/test_data/fastqs/fastqs2/
GATK_BUNDLE=/opt/databases/GATK_bundle/2.8/hg19
HG19=${GATK_BUNDLE}/ucsc.hg19.fasta
#scripts=$WORK'/scripts'
NSLOTS=5
#ls $WORK | cut -d "_" -f 1 | uniq -d | tr " " "\n" > $scripts/samples.txt

#SAMPLE_NUMBER=$(wc -l $scripts/samples.txt | cut -d " " -f 1)

gtfFile=/sharedlustre/NICR/human_genomes/hg19/gencode.v19.annotation.gtf
#gtfFile=/sharedlustre/NICR/human_genomes/hg19/gencode.v27lift37.annotation.gtf
gtfFlatFile=/sharedlustre/NICR/human_genomes/hg19/gencode.v19.annotation.gff
#gtfFlatFile=/sharedlustre/NICR/human_genomes/hg19/gencode.v27lift37.annotation.gff3
genomeFasta=/opt/databases/GATK_bundle/2.8/hg19/ucsc.hg19.fasta
houseBed=/sharedlustre/NICR/human_genomes/hg19/hg19.HouseKeepingGenes.bed
refBed=/sharedlustre/NICR/human_genomes/hg19/hg19.refseq.bed.12
genomeDirstar=/sharedlustre/NICR/human_genomes/hg19/STAR_hg19_genome
genomeDir=/sharedlustre/NICR/human_genomes/hg19/hg19_Gencode14.overhang75
dexSeq=/sharedlustre/NICR/RNAseq_pipelines/DW_RNAseq/dexseq_count.py
star=/opt/STAR_2.3.0e.Linux_x86_64/STAR
genomeFasta=/opt/databases/GATK_bundle/2.8/hg19/ucsc.hg19.fasta
picard=/opt/picard-tools-1.118/picard-tools-1.118
gatk=/opt/GenomeAnalysisTK.jar
BUNDLE_DIR="/opt/databases/GATK_bundle/2.8/hg19"
goldStan="$BUNDLE_DIR/Mills_and_1000G_gold_standard.indels.hg19.vcf"
dbSNP="$BUNDLE_DIR/dbsnp_138.hg19.vcf"


alignments=$WORK'/alignments'
processedData=$alignments'/${FILENAME[$SGE_TASK_ID-1]}'/processedData'
fastqcOutput=$processedData'/fastqc'
htseqOutput=$processedData'/gene_counts'
bamOutput=$processedData'/bam_files'
dexOutput=$processedData'/exon_counts'
cuffQuantOut=$processedData'/cuffquant'
rseqcOutput=$processedData'/rseqc'
vcf=$processedData'/vcf'
secPassDir=$processedData/${FILENAME[$SGE_TASK_ID-1]}_2pass
sampleGenomeDir=$processedData'/${FILENAME[$SGE_TASK_ID-1]}2pass
output_gatk=$processedData'/gatk'
fusioncatcher_out=$processedData'/${FILENAME[$SGE_TASK_ID-1]}_fussion_out
miso_out=$processedData'/${FILENAME[$SGE_TASK_ID-1]}_miso_isoforms
express_out=$processedData'/${FILENAME[$SGE_TASK_ID-1]}_express_isoforms
outCuff=$processedData'/Cuff_output'
#fussionCatcherInput=$procesedData'/fussioncatcher_input'
