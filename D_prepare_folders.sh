
#!/bin/bash
#$ -N prep_folders
#### use this script to create the samples.txt file and folders that the launcher script and function files will read
#### all scripts in this pipeline are lettered alphabetically so that the order in which they should be run is apparent
#### The working directory is set in the variables script you can source that and it will be carried to all subsequent scripts

#### Source the variables script relevant to the genomes you want to use for alignments, unhash the one you want to use
scripts=/user/scripts

source ${scripts}/C_variables_hg19.sh
#source ${scripts}/C_variables_hg38.sh

#### list all files in the directory, in the directory, take the first field of the sample names, pipe to unique to get unique IDs and inseert into the samples.txt file one per line

#ls $WORK/*fastq* | cut -d "_" -f 1 | uniq -d | tr " " "\n" > ${scripts}/samples.txt

### This part of the script will assess if multiple lanes have been uased and move all files with a name pattern into sample specific directories
#ls $WORK/*fastq* > ${scripts}/all_fqs.txt
#all_fqs=$(wc -l ${scripts}/all_fqs.txt | cut -f 1 -d " ")
#samps=$(wc -l ${scripts}/samples.txt | cut -f 1 -d " ")
#if [ "$all_fqs" > "$samps"*2]; then
###1# INALL is an array of fastq file names that are put into per sample directories      
	INALL=( "${WORK}"/* )
		#fsALL=$(find -name "*.fastq.gz")
		for pattern in "${INALL[@]}"; do
		#dir="${pattern%_*L00*_R*_001.fastq.gz}"
		#### unhash this line if files are not gziped	
		#dir="${pattern%_*L00*_R*_001.fastq}"
		#### unhash the line below if format of filename is different e.g. Cancer1_L001_1_001.fastq.gz
		#dir="${pattern%_*L00*_*_001.fastq.gz}"
   		dir="${pattern%_*.fastq.gz}"

			if [ ! -d  ${dir} ]

                then

                                mkdir ${dir}

                else

                                echo "${dir} already exists"

                fi

                mv ${pattern} ${dir}
	#fi
done
#done


#### read the samples in the samples.txt file into an array that will be used to rename output files
readarray -t FILENAME < ${scripts}/samples.txt
#### create the logs folder that error and output files (.e and .o) will be sent to 
logs=$WORK/'logs'
if [ ! -d ${logs} ]; then
        mkdir -p ${logs}
else
        echo "${logs} exists"
fi

#### create the alignments folder for the processed data
#alignments=$WORK/'alignments' #### set in the variables scripts
#if [ ! -d ${alignments} ]; then
#        mkdir -p ${alignments}
#else
#        echo "${alignments} exists"
#fi

#### create a processed data folder specific to each of the sample IDs so that results are correctly ordered
processedData=$WORK/${FILENAME[$SGE_TASK_ID-1]}/'processedData'
if [ ! -d ${processedData} ]; then
        mkdir -p ${processedData}
else
        echo "${processedData} exists"
fi

#fastqcOutput=$processedData'/fastqc'
#htseqOutput=$processedData'/gene_counts'
#bamOutput=$processedData'/bam_files'
#dexOutput=$processedData'/exon_counts'
#cuffQuantOut=$processedData'/cuffquant'
#rseqcOutput=$processedData'/rseqc'
##alignments=$WORK'/alignments'
#vcf=$processedData'/vcf'
#secPassDir=$processedData/${FILENAME[$SGE_TASK_ID-1]}_2pass
#sampleGenomeDir=$processedData'/${FILENAME[$SGE_TASK_ID-1]}2pass
#output_gatk=$processedData'/gatk'
#fusioncatcher_out=$processedData'/${FILENAME[$SGE_TASK_ID-1]}_fussion_out
#miso_out=$processedData'/${FILENAME[$SGE_TASK_ID-1]}_miso_isoforms
#express_out=£processedData'/${FILENAME[$SGE_TASK_ID-1]}_express_isoforms
#outCuffdiff=$processedData'/Cuffdiff_output'
#outCufflinks=$processedData'/Cufflinks_output'
#fussionCatcherInput=$procesedData'/fussioncatcher_input'


#### create the subdirectories of the processedData folder for each of the samples
for j in $fastqcOutput $trimmedFastqs $htseqOutput $bamOutput $dexOutput $cuffQuantOut $rseqcOutput $vcf $secPassDir $output_gatk $fusioncatcher_out $miso_out $express_out $outCuffdiff $outCufflinks $gatk_out $BowtieOut $HomerOut $starGenomeDir;
do
if [ ! -d ${j} ]; then
        mkdir -p ${j}
else
        echo "${j} exists"
fi
done



#for j in $fastqcOutput $htseqOutput $bamOutput $dexOutput $cuffOutput $rseqcOutput $alignments $vcf $secPassDir;
#do
#if [ ! -d ${j}/${FILENAME[$SGE_TASK_ID-1]} ]; then
#        mkdir -p ${j}/${FILENAME[$SGE_TASK_ID-1]}
#else
#        echo "${j}/${FILENAME[$SGE_TASK_ID-1]} exists"
#fi
#done



### the sample_names.txt document is read by the variables.sh script
## use the command wc -l sample_names.txt to determine the number of samples, replace the number of samples (n) in the qsub -t 1:n
# when you submit the scripts use qsub -t 1:73 (there were 73 samples in this run)

