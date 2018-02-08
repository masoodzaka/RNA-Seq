#!/bin/bash
#$ -cwd
#$ -V
#$ -l h_vmem=10G # memory default is 6
#$ -pe smp 1 # CPUs default is 1
#$ -l h_rt=40:00:00
#$ -N RNAseq_func
#$ -q bigmem.q,all.q
#$ -o /home/nlp71/WORKING_DATA/NICR_RNAseq/rna_seq/test_data/fastqs/fastqs2/logs/RNAseq_func.$JOB_ID.$TASK_ID.out
#$ -e /home/nlp71/WORKING_DATA/NICR_RNAseq/rna_seq/test_data/fastqs/fastqs2/logs/RNAseq_func.$JOB_ID.$TASK_ID.err


#  gatk_rna_pipeline.sh
#
#
#  Created by Daniel Williamson on 30/07/2014.
#

NO_ARGS=0
E_OPTERROR=85

if [ $# -eq "$NO_ARGS" ]    # Script invoked with no command-line args?
then
  echo "Usage: `basename $0` options (-mnopqrs)"
  exit $E_OPTERROR          # Exit and explain usage.
                            # Usage: scriptname -options
                            # Note: dash (-) necessary
fi

##############################################
source ${scripts}/C_variables_hg38.sh

readarray -t FILENAME < ${scripts}/samples.txt

WORK=/home/nlp71/WORKING_DATA/NICR_RNAseq/rna_seq/test_data/fastqs/fastqs2

FASTQ1=${WORK}/${FILENAME[$SGE_TASK_ID-1]}_1.fastq.gz
FASTQ2=${WORK}/${FILENAME[$SGE_TASK_ID-1]}_2.fastq.gz


function pipeline_align {
#local FASTQ1=`(ls ${WORK}/${FILENAME[$SGE_TASK_ID-1]}*_*1*.fastq.gz)`
#FASTQ2=`(ls ${WORK}/${FILENAME[$SGE_TASK_ID-1]}*_2*.fastq.gz)`
#local readFile1=$1
#star=$3
#alignments=$4
#genomeDir=$5
#cores=$6
#local FASTQ2=`(ls ${WORK}/${FILENAME[$SGE_TASK_ID-1]}*_*2*.fastq.gz)`
#local readFile2=$(echo $readFile1 | sed 's/R1_001.fastq.gz/R2_001.fastq.gz/g')
echo $FASTQ1
echo $FASTQ2
#local outputName=$(basename $(echo $readFile1 | sed 's/_R1_001.fastq.gz/_/g') )
#local outputDir=$(basename $(dirname $(echo $readFile1 | sed 's/_R1_001.fastq.gz/_/g') ) )
#local fastqcOutput=$2
fastqc -o $fastqcOutput -f fastq $FASTQ1 $FASTQ2 -t $cores

echo "STAR --genomeDir ${genomeDirstar}  --readFilesIn ${FASTQ1} ${FASTQ2}"

STAR \
--genomeDir $genomeDirstar \ 
--readFilesIn $FASTQ1 $FASTQ2 \
--outFilterMismatchNmax 2 \
#####--genomeLoad LoadAndKeep \
--genomeLoad NoSharedMemory \
--outFileNamePrefix $alignments/${FILENAME[$SGE_TASK_ID-1]} \
--outSAMstrandField intronMotif \
--readFilesCommand zcat \
--runThreadN $cores

samtools view -bS $alignments/${FILENAME[$SGE_TASK_ID-1]}'Aligned.out.sam' | samtools sort -n - $alignments/${FILENAME[$SGE_TASK_ID-1]}'Aligned.out.sorted.bam'
#rm $alignments/${FILENAME[$SGE_TASK_ID-1]}'Aligned.out.sam'
}

#function pipeline_merge {
#mergeFiles=$(find ${FILENAME[$SGE_TASK_ID-1]} -name "*Aligned.out.sorted.bam" -print)
#    #local outName=$(basename $1)
#    #echo $outName
#    samtools merge -fn $alignments/${FILENAME[$SGE_TASK_ID-1]}".merged.bam" $mergeFiles
#}

#### To be edited with DAN 
#function pipeline_process {
#local outName=$(echo $(basename $1) | sed 's/.merged.bam//g' )
#samtools view ${FILENAME[$SGE_TASK_ID-1]}.bam | htseq-count -s 'reverse' -m intersection-nonempty - $alignments > $gtfFile/${FILENAME[$SGE_TASK_ID-1]}'_gene_counts.txt'
#python $6 -p 'yes' -f bam -s 'reverse' $genomeDir ${FILENAME[$SGE_TASK_ID-1]} $htseqOutput/${FILENAME[$SGE_TASK_ID-1]}'_exon_counts.txt'  #### what is variable $6
#}

#function rseqc_process {
#local bamDir=$( dirname $1 )
#local bamName=$(basename $1)
#local bamSortName=$( echo $bamName | sed 's/.merged.bam/.merged.sort/g' )
#local bamFile=$bamDir/$bamSortName.bam
#samtools sort $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.bam.sorted'
#samtools index $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.bam.sorted'
#local outDir=$2 ???  ### $bamOutput
#local outCuff=$3  ??? ## $outCuffdiff or $cuffQuantOut ### set to cuffquantout
#local gtfFile=$4
#local houseBed=$5
#local refBed=$6
#outDir=$rseqcOutput ????

#cuffquant -o $cuffQuantOut -u --library-type fr-firststrand $gtfFile $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.bam.sorted'
#geneBody_coverage.py -i $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.bam.sorted' -o $rseqcOutput/${FILENAME[$SGE_TASK_ID-1]}'.bam' -r $houseBed
#mismatch_profile.py -l 101 -i $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.bam.sorted' -o $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'_mismatches.txt'  #### check this with DAN .bam changed to should it be mismatches.txt???
#inner_distance.py -i $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.bam.sorted' -o $outDir/${FILENAME[$SGE_TASK_ID-1]}'_inner_dist.txt' -r $refBed #### .bam changed to _inner_dist.txt
#junction_annotation.py -i $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.bam.sorted' -o $outDir/${FILENAME[$SGE_TASK_ID-1]}'_junc_annot' -r $refBed  ### .bam changed to _junc_annot
#junction_saturation.py -i $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.bam.sorted' -r $refBed -o $outDir/${FILENAME[$SGE_TASK_ID-1]}'_junction_sat.txt'
#deletion_profile.py  -i $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.bam.sorted' -l 101  -o $outDir/${FILENAME[$SGE_TASK_ID-1]}'_deletions.txt' ### deletions added check this 
#clipping_profile.py -i $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.bam.sorted'-s "PE" -o $outDir/${FILENAME[$SGE_TASK_ID-1]}'_clipping.txt'  ### clipping added check this 
#bam_stat.py -i $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.bam.sorted' > $outDir/${FILENAME[$SGE_TASK_ID-1]}'_stats.txt'
#infer_experiment.py -r $refBed -i ${FILENAME[$SGE_TASK_ID-1]}'.bam' > $outDir/${FILENAME[$SGE_TASK_ID-1]}'_infer.txt'
#insertion_profile.py -s "PE" -i ${FILENAME[$SGE_TASK_ID-1]}'.bam'-o $outDir/${FILENAME[$SGE_TASK_ID-1]}'_insertions.txt'  #### insertions added to file extension, should sorted bam be used? 
#read_distribution.py  -i ${FILENAME[$SGE_TASK_ID-1]}'.bam' -r $refBed > $outDir/${FILENAME[$SGE_TASK_ID-1]}'_dist.txt'
#read_duplication.py -i ${FILENAME[$SGE_TASK_ID-1]}'.bam' -o $outDir/${FILENAME[$SGE_TASK_ID-1]}
#read_GC.py -i ${FILENAME[$SGE_TASK_ID-1]}'.bam' -o $outDir/${FILENAME[$SGE_TASK_ID-1]}'_GC.txt'
#read_NVC.py -i ${FILENAME[$SGE_TASK_ID-1]}'.bam' -o $outDir/${FILENAME[$SGE_TASK_ID-1]}'_NVC.txt'
#read_quality.py -i ${FILENAME[$SGE_TASK_ID-1]}'.bam' -o $outDir/${FILENAME[$SGE_TASK_ID-1]}'_quality'
#RNA_fragment_size.py -r $refBed -i ${FILENAME[$SGE_TASK_ID-1]}'.bam'  > $outDir/${FILENAME[$SGE_TASK_ID-1]}'_fragsize.txt'
#RPKM_saturation.py -r $refBed -d '1+-,1-+,2++,2--' -i $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.bam' -o $outDir/${FILENAME[$SGE_TASK_ID-1]}'_RPKM_sats.txt'
#cd $outDir
#tin.py -i ${FILENAME[$SGE_TASK_ID-1]}'.bam' -r $refBed
#cd $bamDir  ### $bamOutput
#mv $bamDir/abundances.cxb $outCuff/$bamName.cxb
#mv $bamOutput/abundances.cxb $cuffQuantOut/${FILENAME[$SGE_TASK_ID-1]}'.bam.cxb' ### cuffQuantOut or cuffdiffOut
#rm $bamFile
#rm $( echo $bamFile | sed 's/.bam/.bai/g' )
#}


function two_pass_align {
  #local readFile1=$1
  #local star=$2
  #local genomeFasta=$3
  #local alignments=$4
  #local readFile2=$(echo $readFile1 | sed 's/R1_001.fastq.gz/R2_001.fastq.gz/g')
  #local cores=$5
  #echo $readFile1
  #echo $readFile2
  #local outputName=$(basename $(echo $readFile1 | sed 's/_R1_001.fastq.gz/_/g') )
  #local baseDir=$( dirname $readFile1)
  #local outputDir=$(basename $(dirname $(echo $readFile1 | sed 's/_R1_001.fastq.gz/_/g') ) )
  #local secPassDir=$alignments/$outputDir/$outputName'2pass'
  #mkdir $secPassDir

 #local FASTQ1=`(ls ${WORK}/${FILENAME[$SGE_TASK_ID-1]}*_1*.fastq.gz)`
  #local FASTQ2=`(ls ${WORK}/${FILENAME[$SGE_TASK_ID-1]}*_2*.fastq.gz)`
  
  $star --runMode genomeGenerate   \
  --genomeDir $secPassDir   \ ### should this be the star genome directory???
  --genomeFastaFiles $genomeFasta   \
  --sjdbFileChrStartEnd $alignments/${FILENAME[$SGE_TASK_ID-1]}'SJ.out.tab'   \
  --sjdbOverhang 75   \
  --runThreadN $cores   
  --genomeLoad LoadAndRemove

  cd $secPassDir
  $star --genomeDir $secPassDir   \
  --readFilesIn $FASTQ1 $FASTQ2   \
  --runThreadN $cores   \
  --readFilesCommand zcat   \
  --genomeLoad LoadAndRemove

  samtools view -bS 'Aligned.out.sam' | samtools sort - 'Aligned.out.sorted'  #### check this should it be named according the sample name i.e. ${FILENAME[$SGE_TASK_ID-1]}'Aligned.out.sam' and ${FILENAME[$SGE_TASK_ID-1]}'Aligned.out.sorted'
  rm chrLength.txt
  rm chrNameLength.txt
  rm chrName.txt
  rm chrStart.txt
  rm Genome
  rm genomeParameters.txt
  rm Log.final.out
  rm Log.out
  rm Log.progress.out
  rm SA
  rm SAindex
  rm sjdbInfo.txt
  rm SJ.out.tab
  rm Aligned.out.sam
}

export -f two_pass_align
### check if this is reminent of the old multiple lane files, safe to remove, or to be run as shown
function pipeline_merge_2pass {
    #local mergeFiles=$(find $1 -name "Aligned.out.sorted.bam" -print)
    #local outName=$(basename $1)
    #echo $outName
	samtools merge -f $alignments/${FILENAME[$SGE_TASK_ID-1]}/${FILENAME[$SGE_TASK_ID-1]}'Aligned.out.sorted.bam' $alignments/${FILENAME[$SGE_TASK_ID-1]}'merged.2.pass.bam'
    #samtools merge -f $1"/"$outName".merged.2.pass.bam" $mergeFiles
    #find $1 -name "Aligned.out.sorted.bam" -exec rm {} \;
}



function readNmark {
  #local bamDir=$( dirname $1 )
  #local bamName=$(basename $1)
  #local bamSortName=$( echo $bamName | sed 's/.merged.2.pass.bam//g' )
  #local bamFile=$bamDir/$bamSortName.rg_added.bam
  #local dedupFile=$bamDir/$bamSortName.dedupped.bam
  #local splitFile=$bamDir/$bamSortName.split.bam
  #local realignFile=$bamDir/$bamSortName.split.realigned.bam
  #local recalTable=$bamDir/$bamSortName.recal_data.table
  #local postRecalTable=$bamDir/$bamSortName.post_recal_data.table
  #local recalFile=$bamDir/$bamSortName.recal_reads.bam
  #local recalPlot=$bamDir/$bamSortName.recalibration_plots.pdf
  #local vcfFile=$bamDir/$bamSortName.output.vcf  #### need to set this $vcf is the folder variable 
  #local vcfFiltFile=$bamDir/$bamSortName.output.filtered.vcf
  #local picard=$2
  #local genomeFasta=$3
  #local goldStan=$4
  #local dbSNP=$5
  #local vcfOutput=$6
  #local gatk=$7
  #
  #cd $bamOutput
  #2. Add read groups, sort, mark duplicates, and create index
  #The above step produces a SAM file, which we then put through the usual Picard processing steps: adding read group #information, sorting, marking duplicates and indexing.
  java -jar $picard/AddOrReplaceReadGroups.jar I=$bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.bam' O=$bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.rg_added.bam' SO=coordinate RGID=id RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=sample
  java -jar $picard/MarkDuplicates.jar I=$bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.bam' O=$bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.dedupped.bam' CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics
  rm ${FILENAME[$SGE_TASK_ID-1]}'.bam'  #### rm which bam file, possible to keep until the end, incase any jobs fail, that way can re-run.  Run a cleanup after run_audit.

  java -jar $gatk -T SplitNCigarReads -R $genomeFasta -I $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.dedupped.bam' -o $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.split.bam' -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
  rm ${FILENAME[$SGE_TASK_ID-1]}'.dedupped.bam' #### possible to keep until the end, incase any jobs fail, that way can re-run.  Run a cleanup after run_audit.

  #The RealignerTargetCreator step would need to be done just once for a single set of indels; so as long as the set of known indels doesn't change, the output.intervals file from below would never need to be recalculated.

  #java -Xmx1g -jar $gatk \
  #  -T RealignerTargetCreator \
  #  -R $genomeFasta \
  #  -o /home/dan/gatk/output.intervals \
  #  -known $goldStan -nt 24
mkdir $bamOutput/temp
  java -Xmx4g -Djava.io.tmpdir=$bamOutput/temp \
    -jar $gatk \
    -I $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.split.bam'\
    -R $genomeFasta \
    -T IndelRealigner \
    -targetIntervals /home/dan/gatk/output.intervals \ #### check with Dan where this is derived from ### from STAR genome generation if so is it now in $secPassDir???
    -o $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.split.realigned.bam' \
    -known $goldStan \
    --consensusDeterminationModel KNOWNS_ONLY  -LOD 0.4
	
	#The IndelRealigner step needs to be run on every bam file.  ### extracted from the gatk_gsub.sh script
#mkdir $1/temp
#java -Xmx4g -Djava.io.tmpdir=$1/temp \
  -jar $gatk \
  -I split.bam \
  -R $genomeFasta \
  -T IndelRealigner \
  -targetIntervals $output_int \
  -o split.realigned.bam \
  -known $goldStan 
  #--consensusDeterminationModel KNOWNS_ONLY  -LOD 0.4

#rm split.bam

  rm $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.split.bam'

  java -jar $gatk -T BaseRecalibrator -R $genomeFasta -I $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.split.realigned.bam' -knownSites $dbSNP -knownSites $goldStan -o $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.recal_data.table' #-nct 8
  java -jar $gatk -T BaseRecalibrator -R $genomeFasta -I $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.split.realigned.bam' -knownSites $dbSNP -knownSites $goldStan -BQSR $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.recal_data.table' -o $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.post_recal_data.table' #-nct 8
  #java -jar $gatk -T BaseRecalibrator -R $genomeFasta -I $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.split.realigned.bam' -knownSites $dbSNP -knownSites $goldStan -BQSR $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.recal_data.table' -o $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.post_recal_data.table' -nct $4 #### what is $4 this is from gatk below
  #java -jar $gatk -T BaseRecalibrator -R $genomeFasta -I split.realigned.bam -knownSites $dbSNP -knownSites $goldStan -o recal_data.table -nct $4
#java -jar $gatk -T AnalyzeCovariates -R $genomeFasta -before $recalTable -after $postRecalTable -plots $recalPlot
  java -jar $gatk -T PrintReads -R $genomeFasta -I $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.split.realigned.bam' -BQSR $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.recal_data.table' -o $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.recal_reads.bam' #-nct 8

  rm $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.split.realigned.bam'

  java -jar $gatk -T HaplotypeCaller -R $genomeFasta -I $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.recal_reads.bam' -recoverDanglingHeads -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o $vcf/${FILENAME[$SGE_TASK_ID-1]}'.output.vcf' #-nct 8
  java -jar $gatk -T VariantFiltration -R $genomeFasta -V $vcf/${FILENAME[$SGE_TASK_ID-1]}'.output.vcf' -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o $vcf/${FILENAME[$SGE_TASK_ID-1]}'.output.filtered.vcf' 
#### done to here
 #mv ${FILENAME[$SGE_TASK_ID-1]}'.output.vcf' $vcf/

  #mv $vcfFile $vcfOutput/$(basename $vcfFile)
  #mv $vcfFiltFile $vcfOutput/$(basename $vcfFiltFile)
  #mv $recalFile $vcfOutput/$(basename $recalFile)
  #mv $( echo $recalFile| sed 's/.bam/.bai/g' ) $vcfOutput/$(basename $recalFile | sed 's/.bam/.bai/g')
  #mv $vcfFile.idx $vcfOutput/$(basename $vcfFile).idx
  #mv $vcfFiltFile.idx $vcfOutput/$(basename $vcfFiltFile).idx
  rm $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.split.bam'  #### leave remove commands until a cleanup run, after audit incase jobs fail this will prevent re-running things??
  rm $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.dedupped.bam'
  rm $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.split.realigned.bam' 
  #rm $recalTable
  #rm $postRecalTable
  find . -name "*_2pass" -exec rm -r {} \;  #### check this remove statement as new twopass directory is _2pass, should this be removed??
}

function run_gatk {   #### Check this with dan as running GATK called in previous function readNmark which runs on sorted BAM files 

#wd=$1
#outputDir=$2

#sampleName=$( basename $wd )
#finalOut=$outputDir/$sampleName
#mkdir $finalOut

#mkdir $TMPDIR/tempFastq
#mkdir $TMPDIR/output
#temp_dir=$TMPDIR/output

#cp -Lr $wd $TMPDIR/tempFastq
#cd $TMPDIR/tempFastq/$sampleName

local FASTQ1=`(ls ${WORK}/${FILENAME[$SGE_TASK_ID-1]}*_1*.fastq.gz)`  ### why sorting the fastq files to run GATK here but ReadNmark function running GTAK on sorted BAM, is this step needed?
local FASTQ2=`(ls ${WORK}/${FILENAME[$SGE_TASK_ID-1]}*_2*.fastq.gz)`  #### Also should the output from this be sent to $gatk_out, what about in previous readNmark as currently sent to $bamOutput
#read1=( $(find $TMPDIR/tempFastq/$sampleName -name "*R1.fastq.gz") )
#read2=( $(find $TMPDIR/tempFastq/$sampleName -name "*R2.fastq.gz") )

sortRead1=( $(
    for el in "${FASTQ1[@]}"
    do
        echo "$el"
    done | sort) )

sortRead2=( $(
    for el in "${FASTQ2[@]}"
    do
        echo "$el"
    done | sort) )

read1CMD=$(IFS=,; echo "${sortRead1[*]}")
read2CMD=$(IFS=,; echo "${sortRead2[*]}")

run_gatk_rna_seq $WORK/$read1CMD $WORK/$read2CMD 5 $gatk_out/${FILENAME[$SGE_TASK_ID-1]}  ## check with DAN function has been removed as repeating previous steps
}
export -f run_gatk


##function run_gatk_rna_seq {   ### check with DAN as this seems to repeat the previous steps of genome alignments in pipeline align and 2passalign and gatk is also run in the readNmark function 

#one=$temp_dir

#htseqOutput=$1'/gene_counts'
  #htseqOutput=$one'/gene_counts'
  #mkdir $htseqOutput
 # dexOutput=$1'/exon_counts'
  #dexOutput=$one'/exon_counts'
  #mkdir $dexOutput


# 1) STAR uses genome index files that must be saved in unique directories. The human genome index was built from the FASTA file hg19.fa as follows:

#one='/home/dan/gatk/NMB108'
#two='/home/dan/gatk/NMB108/NMB108_GAGTGG_L006_R1_001.fastq.gz'
#three='/home/dan/gatk/NMB108/NMB108_GAGTGG_L006_R2_001.fastq.gz'
#four=20

#genomeDir=/home/dan/gatk_genome
#mkdir $genomeDir
#star --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles /home/data/genomes/Homo_sapiens/gatk/hg19/ucsc.hg19.fasta\  --runThreadN 40 --genomeLoad LoadAndRemove


#2) Alignment jobs were executed as follows:
#runDir=$1
#runDir=$one
#mkdir $runDir
#cd $runDir
#$star --genomeDir $genomeDir --genomeLoad Remove
#$star --genomeDir $genomeDir --readFilesIn $2 $3 --runThreadN $4 --readFilesCommand zcat --genomeLoad NoSharedMemory
#star --genomeDir $genomeDir --readFilesIn $read1CMD $read2CMD --runThreadN 2 --readFilesCommand zcat --genomeLoad LoadAndRemove


#3) For the 2-pass STAR, a new index is then created using splice junction information contained in the file SJ.out.tab from the first pass:

#genomeDir=$1/2pass
#genomeDir=$one/2pass
#mkdir $genomeDir
#$star --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $genomeFasta --sjdbFileChrStartEnd $1/SJ.out.tab --sjdbOverhang 75 --runThreadN $4 --genomeLoad NoSharedMemory

#star --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles /home/data/genomes/Homo_sapiens/gatk/hg19/ucsc.hg19.fasta --sjdbFileChrStartEnd $one/SJ.out.tab --sjdbOverhang 75 --runThreadN $four --genomeLoad LoadAndRemove

#4) The resulting index is then used to produce the final alignments as follows:

#runDir=$1/2pass
##runDir=$one/2pass
#cd $runDir
##$star --genomeDir $genomeDir --readFilesIn $2 $3 --runThreadN $4 --readFilesCommand zcat --genomeLoad NoSharedMemory --outSAMstrandField intronMotif
#star --genomeDir $genomeDir --readFilesIn $two $three --runThreadN $four --readFilesCommand zcat --genomeLoad LoadAndRemove

#samtools view -bS 'Aligned.out.sam' | samtools sort -n - 'Aligned.out.sorted.temp'
#pipeline_process Aligned.out.sorted.temp.bam $htseqOutput $dexOutput $gtfFile $gtfFlatFile $dexSeq
#rm Aligned.out.sorted.temp.bam
#samtools view -bS 'Aligned.out.sam' | samtools sort - 'Aligned.out.sorted'

#2. Add read groups, sort, mark duplicates, and create index
#The above step produces a SAM file, which we then put through the usual Picard processing steps: adding read group #information, sorting, marking duplicates and indexing.

#java -jar $picard AddOrReplaceReadGroups I=Aligned.out.sam O=rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=sample

#rm Aligned.out.sam

#java -jar $picard MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics

#rm rg_added_sorted.bam

#java -jar /home/dan/Software/picard-tools-1.118/CreateSequenceDictionary.jar R= /home/#data/genomes/Homo_sapiens/gatk/hg19/ucsc.hg19.fasta O=/home/data/genomes/Homo_sapiens/#UCSC/hg19/Sequence/genome/hg19.dict


#In the future we plan to integrate this into the GATK engine so that it will be done automatically where appropriate, but for now it needs to be run as a separate step.

#At this step we also add one important tweak: we need to reassign mapping qualities, because STAR assigns good alignments a MAPQ of 255 (which technically means “unknown” and is therefore meaningless to GATK). So we use the GATK’s ReassignOneMappingQuality read filter to reassign all good alignments to the default value of 60. This is not ideal, and we hope that in the future RNAseq mappers will emit meaningful quality scores, but in the meantime this is the best we can do. In practice we do this by adding the ReassignOneMappingQuality read filter to the splitter command.

#Please note that we recently (6/11/14) edited this to fix a documentation error regarding the filter to use. See this announcement for details.

#Finally, be sure to specify that reads with N cigars should be allowed. This is currently still classified as an #"unsafe" option, but this classification will change to reflect the fact that this is now a supported option for RNAseq #processing.
#java -jar $gatk -T SplitNCigarReads -R $genomeFasta -I $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.dedupped.bam' -o $bamOutput/${FILENAME[$SGE_TASK_ID-1]}'.split.bam' -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
#java -jar $gatk -T SplitNCigarReads -R $genomeFasta -I dedupped.bam -o split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

#rm dedupped.bam

#indel Realignment (optional) as per dnaseq

#Running the Indel Realigner only at known sites

#While we advocate for using the Indel Realigner over an aggregated bam using the full Smith-Waterman alignment algorithm, it will work for just a single lane of sequencing data when run in -knownsOnly mode. Novel sites obviously won't be cleaned up, but the majority of a single individual's short indels will already have been seen in dbSNP and/or 1000 Genomes. One would employ the known-only/lane-level realignment strategy in a large-scale project (e.g. 1000 Genomes) where computation time is severely constrained and limited. We modify the example arguments from above to reflect the command-lines necessary for known-only/lane-level cleaning.

#The RealignerTargetCreator step would need to be done just once for a single set of indels; so as long as the set of known indels doesn't change, the output.intervals file from below would never need to be recalculated.

# java -Xmx1g -jar $gatk \
#  -T RealignerTargetCreator \
#  -R $genomeFasta \
#  -o /home/ndw31/WORKING_DATA/gatk/output.intervals \
#  -known goldStan -nt 24

#The IndelRealigner step needs to be run on every bam file.
#mkdir $1/temp
#java -Xmx4g -Djava.io.tmpdir=$1/temp \
  -jar $gatk \
  -I split.bam \
  -R $genomeFasta \
  -T IndelRealigner \
  -targetIntervals $output_int \
  -o split.realigned.bam \
  -known $goldStan --consensusDeterminationModel KNOWNS_ONLY  -LOD 0.4

#rm split.bam

#Base Recalibration as per dnaseq

#perl -pi -e 's/^/chr/' dbSNP_135.no1000GProduction.vcf
#perl vcfsorter.pl /home/data/genomes/Homo_sapiens/UCSC/hg19/Sequence/genome/hg19.dict dbSNP_135.no1000GProduction.vcf > dbSNP_135.no1000GProduction.sort.vcf

#Creating a recalibrated BAM
############################# ALREADY IN A FUNCTION IN THIS FILE HASHED OUT #########################################
#java -jar $gatk -T BaseRecalibrator -R $genomeFasta -I split.realigned.bam -knownSites $dbSNP -knownSites $goldStan -o recal_data.table -nct $4


#2. Do a second pass to analyze covariation remaining after recalibration
#Action
#Run the following GATK command:

#java -jar $gatk -T BaseRecalibrator -R $genomeFasta -I split.realigned.bam -knownSites $dbSNP -knownSites $goldStan -BQSR recal_data.table -o post_recal_data.table -nct $4


#java -jar /home/dan/Software/GenomeAnalysisTK.jar -T AnalyzeCovariates -R /home/data/genomes/Homo_sapiens/gatk/hg19/ucsc.hg19.fasta -before recal_data.table -after post_recal_data.table -plots recalibration_plots.pdf


#To create a recalibrated BAM you can use GATK's PrintReads with the engine on-the-fly recalibration capability. Here is a typical command line to do so:

#java -jar $gatk -T PrintReads -R $genomeFasta -I split.realigned.bam -BQSR recal_data.table -o recal_reads.bam -nct $4

#rm split.realigned.bam

#6. Variant calling

#Finally, we have arrived at the variant calling step! Here, we recommend using HaplotypeCaller because it is performing much better in our hands than UnifiedGenotyper (our tests show that UG was able to call less than 50% of the true positive indels that HC calls). We have added some functionality to the variant calling code which will intelligently take into account the information about intron-exon split regions that is embedded in the BAM file by SplitNCigarReads. In brief, the new code will perform “dangling head merging” operations and avoid using soft-clipped bases (this is a temporary solution) as necessary to minimize false positive and false negative calls. To invoke this new functionality, just add -recoverDanglingHeads -dontUseSoftClippedBases to your regular HC command line. Also, we found that we get better results if we lower the minimum phred-scaled confidence threshold for calling variants on RNAseq data, so we use a default of 20 (instead of 30 in DNA-seq data).

#java -jar $gatk -T HaplotypeCaller -R $genomeFasta -I recal_reads.bam -recoverDanglingHeads -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o output.vcf -nct $4

#rm recal_reads.bam


#7. Variant filtering

#To filter the resulting callset, you will need to apply hard filters, as we do not yet have the RNAseq training/truth resources that would be needed to run variant recalibration (VQSR).

#We recommend that you filter clusters of at least 3 SNPs that are within a window of 35 bases between them by adding -window 35 -cluster 3 to your command. This filter recommendation is specific for RNA-seq data.

#As in DNA-seq, we recommend filtering based on Fisher Strand values (FS > 30.0) and Qual By Depth values (QD < 2.0).


#java -jar $gatk -T VariantFiltration -R $genomeFasta -V output.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o output.filtered.vcf

#cp output.vcf $1/output.vcf
#cp output.filtered.vcf $1/output.filtered.vcf
#cp recal_reads.bam $1/recal_reads.bam
#cp output.vcf.idx $1/output.vcf.idx
#cp output.filtered.vcf.idx $1/output.filtered.vcf.idx

#file=Aligned.out.sorted.bam
#minimumsize=92
#actualsize=$(wc -c <"$file")
#if [ $actualsize -eq $minimumsize ]; then
#    echo size equals $minimumsize bytes
#    exit 1
#else
#    echo size is greater than $minimumsize bytes
#fi

#mv Aligned.out.sorted.bam $1/$6.bam
#samtools index $1/$6.bam


#cp Aligned.out.sorted.bam $1/Aligned.out.sorted.bam
##samtools index $1/Aligned.out.sorted.bam
#cd $1
#rm -r $1/2pass
#rm -r $1/temp
#rm Log.final.out
#rm Log.progress.out
#rm Aligned.out.sam
#rm Log.out
#rm SJ.out.tab
#cp -r $1 $5
}

#export -f run_gatk_rna_seq

#$wd=$WORK

T="$(date +%s)"
echo "Aligning reads now..."
export -f pipeline_align
$WORK pipeline_align {} $FASTQ1 $FASTQ2 $star $alignments $genomeDirstar 3
#find $wd -name "*R1_001.fastq.gz" | parallel --gnu -j 10 pipeline_align {} $fastqcOutput $star $alignments $genomeDir 3  #### should fastqcOutput be the fastq file folder???  Should both fastq files be used in the alignments???
T="$(($(date +%s)-T))"
echo "Time in seconds: ${T}"

#T="$(date +%s)"
#echo "Merging reads now..."
#export -f pipeline_merge
#find $wd/alignments -mindepth 1 -maxdepth 1 -type d -not \( -path "*/processedData" -prune \) -not \( -path "*/alignments" -prune \) | parallel --gnu -j 5 pipeline_merge {} $alignments
#T="$(($(date +%s)-T))"
#echo "Time in seconds: ${T}"

T="$(date +%s)"
echo "processing pipeline..."
export -f pipeline_process
$WORK pipeline_process {} $alignments $htseqOutput $dexOutput $gtfFile $gtfFlatFile $dexSeq
#find $wd/alignments -name "*.merged.bam" | parallel --gnu -j 10 pipeline_process {} $htseqOutput $dexOutput $gtfFile $gtfFlatFile $dexSeq
T="$(($(date +%s)-T))"
echo "Time in seconds: ${T}"

T="$(date +%s)"
echo "processing rseqc pipeline..."
export -f rseqc_process
$alignments rseqc_process {} $rseqcOutput $outCuffdiff $gtfFile $houseBed $refBed
#find $wd/alignments -name "*.merged.bam" | parallel --gnu -j 10 rseqc_process {} $rseqcOutput $outCuffdiff $gtfFile $houseBed $refBed ####$cuffOutput removed as made more specific to outCuffdiff
T="$(($(date +%s)-T))"
echo "Time in seconds: ${T}"

T="$(date +%s)"
echo "Merging reads now..."
export -f two_pass_align
$WORK two_pass_align {} $FASTQ1 $FASTQ2 $star $genomeFasta $alignments 3
#find $wd -name "*R1_001.fastq.gz" | parallel --gnu -j 10 two_pass_align {} $star $genomeFasta $alignments 3
two_pass_align $FASTQ1 $star $genomeFasta  #### $readFile1 changed to FASTQ1, should FASTQ2 be in here also?   
T="$(($(date +%s)-T))"
echo "Time in seconds: ${T}"

T="$(date +%s)"
echo "Merging reads now..."
export -f pipeline_merge_2pass
$alignments pipeline_merge_2pass {}
find $wd/alignments -mindepth 1 -maxdepth 1 -type d -not \( -path "*/processedData" -prune \) | parallel --gnu -j 5 pipeline_merge_2pass {}
T="$(($(date +%s)-T))"
echo "Time in seconds: ${T}"

  #java -jar /home/dan/Software/picard-tools-1.118/CreateSequenceDictionary.jar R= /home/#data/genomes/Homo_sapiens/gatk/hg19/ucsc.hg19.fasta O=/home/data/genomes/Homo_sapiens/#UCSC/hg19/Sequence/genome/hg19.dict

T="$(date +%s)"
echo "processing vcf now..."
export -f readNmark
#mkdir -p $wd'/vcf'
#vcfOutput=$wd'/vcf'
readNmark {} $alignments/${FILENAME[$SGE_TASK_ID-1]}'merged.2.pass.bam' $picard $genomeFasta $goldStan $dbSNP $vcfOutput $gatk
#find $wd/alignments -name "*.merged.2.pass.bam" | parallel --gnu -j 10 readNmark {} $picard $genomeFasta $goldStan $dbSNP $vcf $gatk
T="$(($(date +%s)-T))"
echo "Time in seconds: ${T}"



