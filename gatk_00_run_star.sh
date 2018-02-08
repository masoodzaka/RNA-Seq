#!/bin/bash
#$ -cwd -V
#$ -pe smp 20
#$ -l h_vmem=40G
# $ -M louise.pease2@ncl.ac.uk
#$ -m a
#$ -N star
#$ -hold_jid starprep
#$ -o /user/logs/star.$JOB_ID.$TASK_ID.out
#$ -e /user/logs/star.$JOB_ID.$TASK_ID.err

module add apps/STAR/2.5.2b

start_time=$(date)
start_secs=$(date +%s)
echo "### Logging Info ###"
echo "### Job: ${JOB_ID}: ${JOB_NAME} ###"
echo "### Array ID: ${SGE_TASK_ID} ###"
echo "### Job Occupying ${NSLOTS} slots ###"
echo "### Job running on ${HOSTNAME} ###"
echo "### Started at: ${start_time} ###"
echo

source /user/scripts/C_variables_hg19.sh

readarray -t FILENAME < ${scripts}/samples.txt

FASTQ1=${WORK}/${FILENAME[$SGE_TASK_ID-1]}/${FILENAME[$SGE_TASK_ID-1]}_1.fastq.gz
FASTQ2=${WORK}/${FILENAME[$SGE_TASK_ID-1]}/${FILENAME[$SGE_TASK_ID-1]}_2.fastq.gz

# 1) STAR uses genome index files that must be saved in unique directories. The human genome index was built from the FASTA file hg19.fa as follows:
# DO ONCE FOR ALL RUNS

#GENOME=${WORK}/STAR_genome

#mkdir ${GENOME}
#STAR --runMode genomeGenerate --genomeDir ${GENOME} \
#    --genomeFastaFiles ${HG19} \
#    --runThreadN ${NSLOTS}

# 2) Alignment jobs were executed as follows:

runDir=${TMPDIR}/1pass
mkdir $runDir
cd $runDir
STAR --genomeDir ${genomeDirstar} --readFilesIn ${FASTQ1} ${FASTQ2} --readFilesCommand zcat --runThreadN ${cores}

# 3) For the 2-pass STAR, a new index is then created using splice junction information contained in the file SJ.out.tab from the first pass:

genomeDir=${TMPDIR}/hg19_2pass
mkdir ${genomeDir}
STAR --runMode genomeGenerate --genomeDir ${genomeDir} --genomeFastaFiles ${genomeFasta} \
    --sjdbFileChrStartEnd ${TMPDIR}/1pass/SJ.out.tab --sjdbOverhang 75 --runThreadN ${cores}
    
# 4) The resulting index is then used to produce the final alignments as follows:

runDir=${TMPDIR}/2pass
mkdir $runDir
cd $runDir
STAR --genomeDir ${genomeDir} --readFilesIn ${FASTQ1} ${FASTQ2} --readFilesCommand zcat --runThreadN ${cores}

#mv ${runDir} ${WORK}/${FILENAME[$SGE_TASK_ID-1]}/processedData/bam_files
mv ${runDir} ${bamOutput}

end_time=$(date)
end_secs=$(date +%s)
time_elapsed=$(echo "${end_secs} - ${start_secs}" | bc)
echo
echo "### Ended at: ${end_time} ###"
echo "### Time Elapsed: ${time_elapsed} ###"
