# prep1_trim.sh version 1
# Jannetta's bioinformatics scripts repurposed for NICR by Louise Pease October 2017
# The purpose of this script is to trim the reads in the fastq files.
# To determine the appropriate lengths you need to look at the fastqc
# report prepared by the fastqc program (script prep0_fastqc)
#
#!/bin/bash -e
#$ -cwd -V
#$ -pe smp 1-5
#$ -l h_vmem=40G
# $ -M louise.pease2@newcastle.ac.uk
#$ -m a
#$ -N trim
#$ -o /user/logs/trim.$JOB_ID.$TASK_ID.out
#$ -e /user/logs/trim.$JOB_ID.$TASK_ID.err

module add apps/Trimmomatic/0.3.6

start_time=$(date)
start_secs=$(date +%s)
echo "### Logging Info ###"
echo "### Job: ${JOB_ID}: ${JOB_NAME} ###"
echo "### Array ID: ${SGE_TASK_ID} ###"
echo "### Job Occupying ${NSLOTS} slots ###"
echo "### Job running on ${HOSTNAME} ###"
echo "### Started at: ${start_time} ###"
echo

source /user/scripts/C_variables_hg19.sh

readarray -t FILENAME < ${WORK}/new_scripts/samples.txt

FASTQ1=${WORK}/${FILENAME[$SGE_TASK_ID-1]}/${FILENAME[$SGE_TASK_ID-1]}_1.fastq.gz
FASTQ2=${WORK}/${FILENAME[$SGE_TASK_ID-1]}/${FILENAME[$SGE_TASK_ID-1]}_2.fastq.gz
FASTQtrimmed1=${FILENAME[$SGE_TASK_ID-1]}_1.fastq.gz
FASTQtrimmed2=${FILENAME[$SGE_TASK_ID-1]}_2.fastq.gz 
HEADCROP=11
CROP=88
MINLEN=15


runDir=${TMPDIR}/${FILENAME[$SGE_TASK_ID-1]}
mkdir $runDir
cd $runDir

echo java -Xmx4g -Xms4g -jar /opt/software/clustersoft/install/Trimmomatic/Trimmomatic-0.35/trimmomatic-0.35.jar PE ${FASTQ1} ${FASTQ2} ${trimmedFastqs}/${FASTQtrimmed1}.trimmed.fastq.gz ${trimmedFastqs}/${FASTQtrimmed1}.unpaired ${trimmedFastqs}/${FASTQtrimmed2}.trimmed.fastq.gz ${trimmedFastqs}/${FASTQtrimmed2}.unpaired ILLUMINACLIP: TruSeq3-PE.fa:2:30:10 HEADCROP:$HEADCROP MINLEN:$MINLEN
java -Xmx4g -Xms4g -jar /opt/software/clustersoft/install/Trimmomatic/Trimmomatic-0.35/trimmomatic-0.35.jar PE \
	${FASTQ1} ${FASTQ2} \
	${trimmedFastqs}/${FASTQtrimmed1}.trimmed.fastq.gz \
	${trimmedFastqs}/${FASTQtrimmed1}.unpaired \
	${trimmedFastqs}/${FASTQtrimmed2}.trimmed.fastq.gz \
	${trimmedFastqs}/${FASTQtrimmed2}.unpaired \
	ILLUMINACLIP:/opt/software/clustersoft/install/Trimmomatic/Trimmomatic-0.35/adapters/TruSeq3-PE.fa:2:30:10:8:TRUE \
	HEADCROP:$HEADCROP \
	CROP:$CROP \
	MINLEN:$MINLEN

mv $runDir ${trimmedFastqs}

end_time=$(date)
end_secs=$(date +%s)
time_elapsed=$(echo "${end_secs} - ${start_secs}" | bc)
echo
echo "### Ended at: ${end_time} ###"
echo "### Time Elapsed: ${time_elapsed} ###"
