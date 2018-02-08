


module add apps/samtools/1.2 
module add apps/STAR/2.5.2b 
module add apps/picard/1.130 
module add apps/gatk/3.2-protected 
module add apps/python27/2.7.8 
module add libs/python/numpy/1.9.1-python27-2.7.8 
module add libs/python/pysam/0.10.0 
module add libs/python/matplotlib/1.3.1-python27 
module add apps/HTSeq/0.6.1 

#gtfFile='/home/ndw31/WORKING_DATA/data/ucsc/hg19/Annotation/Genes/Gencode25/gencode.v25lift37.annotation.gtf' 
#gtfFlatFile='/home/ndw31/WORKING_DATA/data/ucsc/hg19/Annotation/Genes/Gencode25/gencode.v25lift37.annotation.gff'
#python /home/ndw31/script/dexseq_prepare_annotation.py $gtfFile $gtfFlatFile 
#houseBed=/data/rseqc/hg19.HouseKeepingGenes.bed refBed=/data/rseqc/hg19.refseq.bed.12
#genomeDir='/home/ndw31/WORKING_DATA/data/star/genome_files/star2.5.2b_hg19' 
#dexSeq='/home/ndw31/script/dexseq_count.py' star='STAR' 
#genomeFasta='/home/ndw31/WORKING_DATA/data/genomes/Homo_sapiens/gatk/hg19/ucsc.hg19.fasta'
#picard='runpicard.sh'
picard='/opt/software/clustersoft/install/picard-tools/picard-tools-1.130/picard.jar' 
gatk='/opt/software/clustersoft/install/gatk/GenomeAnalysisTK-3.2-protected/GenomeAnalysisTK.jar' 
goldStan='/home/ndw31/WORKING_DATA/data/genomes/Homo_sapiens/gatk/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf' 
dbSNP='/home/ndw31/WORKING_DATA/data/genomes/Homo_sapiens/gatk/hg19/dbsnp_138.hg19.vcf' 
output_path='/home/ndw31/WORKING_DATA/output_gatk' 
search_path='/home/nlp71/WORKING_DATA/NICR_RNAseq/rna_seq/test_data/fastqs/' 
#output_int='/home/ndw31/WORKING_DATA/gatk/output.intervals'
#mkdir -p /home/ndw31/WORKING_DATA/data/ucsc/hg19/Annotation/Genes/Gencode17 rsync -avzhe ssh dan@instinct:/data/ucsc/hg19/Annotation/Genes/Gencode17/gencode.v17.annotation.gtf 
#/home/ndw31/WORKING_DATA/data/ucsc/hg19/Annotation/Genes/Gencode17/gencode.v17.annotation.gtf rsync -avzhe ssh dan@instinct:/data/ucsc/hg19/Annotation/Genes/Gencode17/gencode.v17.annotation.flat.gtf 
#/home/ndw31/WORKING_DATA//data/ucsc/hg19/Annotation/Genes/Gencode17/gencode.v17.annotation.flat.gtf mkdir -p /home/ndw31/WORKING_DATA/data/genomes/Homo_sapiens/gatk/hg19 rsync -avzhe ssh 
#dan@instinct:/data/genomes/Homo_sapiens/gatk/hg19 /home/ndw31/WORKING_DATA/data/genomes/Homo_sapiens/gatk/hg19 mkdir -p /home/ndw31/WORKING_DATA/data/star/genome_files/hg19_Gencode14.overhang75 rsync -avzhe ssh 
#dan@instinct:/data/star/genome_files/hg19_Gencode14.overhang75 /home/ndw31/WORKING_DATA/data/star/genome_files/hg19_Gencode14.overhang75 mkdir -p /home/ndw31/WORKING_DATA/fastq/A1905/A1905_B5/Sample_NMB419 rsync -avzhe 
#ssh dan@nicrcompsvr:/media/rna_seq/rawData/A1905/A1905_B5/Sample_NMB419 /home/ndw31/WORKING_DATA/fastq/A1905/A1905_B5/Sample_NMB419
function run_gatk { wd=$1 outputDir=$2 sampleName=$( basename $wd ) finalOut=$outputDir/$sampleName mkdir $finalOut mkdir $TMPDIR/tempFastq mkdir $TMPDIR/output temp_dir=$TMPDIR/output cp -Lr $wd $TMPDIR/tempFastq cd 
$TMPDIR/tempFastq/$sampleName read1=( $(find $TMPDIR/tempFastq/$sampleName -name "*R1.fastq.gz") ) read2=( $(find $TMPDIR/tempFastq/$sampleName -name "*R2.fastq.gz") ) sortRead1=( $(
    for el in "${read1[@]}"
    do
        echo "$el"
    done | sort) ) sortRead2=( $(
    for el in "${read2[@]}"
    do
        echo "$el"
    done | sort) ) read1CMD=$(IFS=,; echo "${sortRead1[*]}") read2CMD=$(IFS=,; echo "${sortRead2[*]}") run_gatk_rna_seq $temp_dir $read1CMD $read2CMD 5 $finalOut $sampleName
}
export -f run_gatk
##
function pipeline_process { local outName=$(echo $(basename $1) | sed 's/.temp.bam//g' ) samtools view $1 | htseq-count -s 'reverse' -m intersection-nonempty - $4 > $2/$outName'_gene_counts.txt' python $6 -p 'yes' -f bam 
-s 'reverse' $5 $1 $3/$outName'_exon_counts.txt'
}
export -f pipeline_process
#pipeline_process $bam $htseqOutput $dexOutput $gtfFile $gtfFlatFile $dexSeq
function run_gatk_rna_seq {
#one=$temp_dir
  htseqOutput=$1'/gene_counts'
  #htseqOutput=$one'/gene_counts'
  mkdir $htseqOutput
  dexOutput=$1'/exon_counts'
  #dexOutput=$one'/exon_counts'
  mkdir $dexOutput
# 1) STAR uses genome index files that must be saved in unique directories. The human genome index was built from the FASTA file hg19.fa as follows: one='/home/dan/gatk/NMB108' 
#two='/home/dan/gatk/NMB108/NMB108_GAGTGG_L006_R1_001.fastq.gz' three='/home/dan/gatk/NMB108/NMB108_GAGTGG_L006_R2_001.fastq.gz' four=20 genomeDir=/home/dan/gatk_genome mkdir $genomeDir star --runMode genomeGenerate 
#--genomeDir $genomeDir --genomeFastaFiles /home/data/genomes/Homo_sapiens/gatk/hg19/ucsc.hg19.fasta\ --runThreadN 40 --genomeLoad LoadAndRemove 2) Alignment jobs were executed as follows:
runDir=$1
#runDir=$one mkdir $runDir
cd $runDir $star --genomeDir $genomeDir --genomeLoad Remove $star --genomeDir $genomeDir --readFilesIn $2 $3 --runThreadN $4 --readFilesCommand zcat --genomeLoad NoSharedMemory
#star --genomeDir $genomeDir --readFilesIn $read1CMD $read2CMD --runThreadN 2 --readFilesCommand zcat --genomeLoad LoadAndRemove 3) For the 2-pass STAR, a new index is then created using splice junction information 
#contained in the file SJ.out.tab from the first pass:
genomeDir=$1/2pass
#genomeDir=$one/2pass
mkdir $genomeDir $star --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $genomeFasta --sjdbFileChrStartEnd $1/SJ.out.tab --sjdbOverhang 75 --runThreadN $4 --genomeLoad NoSharedMemory
#star --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles /home/data/genomes/Homo_sapiens/gatk/hg19/ucsc.hg19.fasta --sjdbFileChrStartEnd $one/SJ.out.tab --sjdbOverhang 75 --runThreadN $four --genomeLoad 
#LoadAndRemove 4) The resulting index is then used to produce the final alignments as follows:
runDir=$1/2pass
#runDir=$one/2pass
cd $runDir $star --genomeDir $genomeDir --readFilesIn $2 $3 --runThreadN $4 --readFilesCommand zcat --genomeLoad NoSharedMemory --outSAMstrandField intronMotif
#star --genomeDir $genomeDir --readFilesIn $two $three --runThreadN $four --readFilesCommand zcat --genomeLoad LoadAndRemove
samtools view -bS 'Aligned.out.sam' | samtools sort -n - 'Aligned.out.sorted.temp' pipeline_process Aligned.out.sorted.temp.bam $htseqOutput $dexOutput $gtfFile $gtfFlatFile $dexSeq rm Aligned.out.sorted.temp.bam samtools 
view -bS 'Aligned.out.sam' | samtools sort - 'Aligned.out.sorted'
#2. Add read groups, sort, mark duplicates, and create index The above step produces a SAM file, which we then put through the usual Picard processing steps: adding read group #information, sorting, marking duplicates and 
#indexing.
java -jar $picard AddOrReplaceReadGroups I=Aligned.out.sam O=rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=sample rm Aligned.out.sam java -jar $picard MarkDuplicates 
I=rg_added_sorted.bam O=dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics rm rg_added_sorted.bam
#java -jar /home/dan/Software/picard-tools-1.118/CreateSequenceDictionary.jar R= /home/#data/genomes/Homo_sapiens/gatk/hg19/ucsc.hg19.fasta O=/home/data/genomes/Homo_sapiens/#UCSC/hg19/Sequence/genome/hg19.dict In the 
#future we plan to integrate this into the GATK engine so that it will be done automatically where appropriate, but for now it needs to be run as a separate step. At this step we also add one important tweak: we need to 
#reassign mapping qualities, because STAR assigns good alignments a MAPQ of 255 (which technically means “unknown” and is therefore meaningless to GATK). So we use the GATK’s ReassignOneMappingQuality read 
#filter to reassign all good alignments to the default value of 60. This is not ideal, and we hope that in the future RNAseq mappers will emit meaningful quality scores, but in the meantime this is the best we can do. In 
#practice we do this by adding the ReassignOneMappingQuality read filter to the splitter command. Please note that we recently (6/11/14) edited this to fix a documentation error regarding the filter to use. See this 
#announcement for details. Finally, be sure to specify that reads with N cigars should be allowed. This is currently still classified as an #"unsafe" option, but this classification will change to reflect the fact that 
#this is now a supported option for RNAseq #processing.
java -jar $gatk -T SplitNCigarReads -R $genomeFasta -I dedupped.bam -o split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS rm dedupped.bam
#indel Realignment (optional) as per dnaseq Running the Indel Realigner only at known sites While we advocate for using the Indel Realigner over an aggregated bam using the full Smith-Waterman alignment algorithm, it will 
#work for just a single lane of sequencing data when run in -knownsOnly mode. Novel sites obviously won't be cleaned up, but the majority of a single individual's short indels will already have been seen in dbSNP and/or 
#1000 Genomes. One would employ the known-only/lane-level realignment strategy in a large-scale project (e.g. 1000 Genomes) where computation time is severely constrained and limited. We modify the example arguments from 
#above to reflect the command-lines necessary for known-only/lane-level cleaning. The RealignerTargetCreator step would need to be done just once for a single set of indels; so as long as the set of known indels doesn't 
#change, the output.intervals file from below would never need to be recalculated.
# java -Xmx1g -jar $gatk \
#  -T RealignerTargetCreator \ -R $genomeFasta \ -o /home/ndw31/WORKING_DATA/gatk/output.intervals \ -known goldStan -nt 24 The IndelRealigner step needs to be run on every bam file.
mkdir $1/temp java -Xmx4g -Djava.io.tmpdir=$1/temp \
  -jar $gatk \
  -I split.bam \
  -R $genomeFasta \
  -T IndelRealigner \
  -targetIntervals $output_int \
  -o split.realigned.bam \
  -known $goldStan --consensusDeterminationModel KNOWNS_ONLY -LOD 0.4 rm split.bam
#Base Recalibration as per dnaseq perl -pi -e 's/^/chr/' dbSNP_135.no1000GProduction.vcf perl vcfsorter.pl /home/data/genomes/Homo_sapiens/UCSC/hg19/Sequence/genome/hg19.dict dbSNP_135.no1000GProduction.vcf > 
#dbSNP_135.no1000GProduction.sort.vcf Creating a recalibrated BAM
java -jar $gatk -T BaseRecalibrator -R $genomeFasta -I split.realigned.bam -knownSites $dbSNP -knownSites $goldStan -o recal_data.table -nct $4
#2. Do a second pass to analyze covariation remaining after recalibration Action Run the following GATK command:
java -jar $gatk -T BaseRecalibrator -R $genomeFasta -I split.realigned.bam -knownSites $dbSNP -knownSites $goldStan -BQSR recal_data.table -o post_recal_data.table -nct $4
#java -jar /home/dan/Software/GenomeAnalysisTK.jar -T AnalyzeCovariates -R /home/data/genomes/Homo_sapiens/gatk/hg19/ucsc.hg19.fasta -before recal_data.table -after post_recal_data.table -plots recalibration_plots.pdf To 
#create a recalibrated BAM you can use GATK's PrintReads with the engine on-the-fly recalibration capability. Here is a typical command line to do so:
java -jar $gatk -T PrintReads -R $genomeFasta -I split.realigned.bam -BQSR recal_data.table -o recal_reads.bam -nct $4 rm split.realigned.bam
#6. Variant calling Finally, we have arrived at the variant calling step! Here, we recommend using HaplotypeCaller because it is performing much better in our hands than UnifiedGenotyper (our tests show that UG was able 
#to call less than 50% of the true positive indels that HC calls). We have added some functionality to the variant calling code which will intelligently take into account the information about intron-exon split regions 
#that is embedded in the BAM file by SplitNCigarReads. In brief, the new code will perform “dangling head merging” operations and avoid using soft-clipped bases (this is a temporary solution) as necessary to 
#minimize false positive and false negative calls. To invoke this new functionality, just add -recoverDanglingHeads -dontUseSoftClippedBases to your regular HC command line. Also, we found that we get better results if we 
#lower the minimum phred-scaled confidence threshold for calling variants on RNAseq data, so we use a default of 20 (instead of 30 in DNA-seq data).
java -jar $gatk -T HaplotypeCaller -R $genomeFasta -I recal_reads.bam -recoverDanglingHeads -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o output.vcf -nct $4 rm recal_reads.bam
#7. Variant filtering To filter the resulting callset, you will need to apply hard filters, as we do not yet have the RNAseq training/truth resources that would be needed to run variant recalibration (VQSR). We recommend 
#that you filter clusters of at least 3 SNPs that are within a window of 35 bases between them by adding -window 35 -cluster 3 to your command. This filter recommendation is specific for RNA-seq data. As in DNA-seq, we 
#recommend filtering based on Fisher Strand values (FS > 30.0) and Qual By Depth values (QD < 2.0).
java -jar $gatk -T VariantFiltration -R $genomeFasta -V output.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o output.filtered.vcf cp output.vcf $1/output.vcf cp 
output.filtered.vcf $1/output.filtered.vcf
#cp recal_reads.bam $1/recal_reads.bam cp output.vcf.idx $1/output.vcf.idx cp output.filtered.vcf.idx $1/output.filtered.vcf.idx
file=Aligned.out.sorted.bam minimumsize=92 actualsize=$(wc -c <"$file") if [ $actualsize -eq $minimumsize ]; then
    echo size equals $minimumsize bytes
    exit 1 else
    echo size is greater than $minimumsize bytes fi mv Aligned.out.sorted.bam $1/$6.bam samtools index $1/$6.bam
#cp Aligned.out.sorted.bam $1/Aligned.out.sorted.bam samtools index $1/Aligned.out.sorted.bam
cd $1 rm -r $1/2pass rm -r $1/temp rm Log.final.out rm Log.progress.out rm Aligned.out.sam rm Log.out
#rm SJ.out.tab
cp -r $1 $5
}
export -f run_gatk_rna_seq date hostname run_gatk $1 $2
# Store sample id (if in file name) as $BASE BASE=`basename $1 .vcf` Change as needed OUTDIR="~/WORKING_DATA/project/stage_1/output" 2nd Arg to pass NO=$2 # you can place 10E6 etc on command line as second arg or hardwire 
#it here
# Sometimes useful to log these and to submit each job Submit_jobs.pl (could easily be shell if you prefer)
