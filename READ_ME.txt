#### These scripts will run an RNA-seq analysis with gatk base recalibration and realignment, STAR genome alignment
#### The results can be anlysed using a variety of methods including:
#### DexSeq = Differential Exon usage in RNASeq data, see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3460195/
#### HTseq = High Throughput sequencing, see https://academic.oup.com/bioinformatics/article/31/2/166/2366196/HTSeq-a-Python-framework-to-work-with-high
#### Miso Mixture of Isoforms probablistic model for RNASeq, isoform analysis, see http://genes.mit.edu/burgelab/miso/
#### Express Isoform analysis
#### FussionCatcher for gene fussion detection in RNAseq data, see https://omictools.com/fusioncatcher-tool
#### Variant calling on sequences (GATK)


#### The pipeline is designed to be modular so that unnecessary or unwanted steps can be removed by hashing out the relevant line of the Launcher script
#### Two variables files are provided so that different genome builds can be used, however most tools currently support the hg19 build at present

#### The preparation scripts are lettered in alphabetical order so the order in which to run the scripts is apparent
#### STEP 1 put all you fastq files in a directory and type pwd the right click to copy the filepath
#### STEP2 Clone the github directory into a scripts folder in the directory containing your fastqs, or paste the fastq file path into and run the script A_copy_scripts from the sharedlustre/NICR/RNAseq_pipelines folder, then type cd and paste the filepath, now cd into the scripts directory
#### STEP 3 In the scripts directory open the sedsetup script and change the filepath accordingly, remember to keep the forward slashes as these act as escape characters 
#### STEP 4 check the format of your FASTQs filenames (Cancer1_R1_L001_001.fastq.gz, Cancer1_R1_L002_001.fast.gz means multiple lanes used), open the D_prepare folders script and select the correct filename format, unhash the correct line as below.
for pattern in "${INALL[@]}"; do
                dir="${pattern%_R*_L00*_001.fastq.gz}"
                #### unhash this line if files are not gziped
                #dir="${pattern%_*L00*_R*_001.fastq}"
                #### unhash the line below if format of filename is different e.g. Cancer1_L001_1_001.fastq.gz
                #dir="${pattern%_*L00*_*_001.fastq.gz}"
                #dir="${pattern%_*.fastq.gz}"

#### STEP 4 open the Launcher_RNAseq.sh script with nano and hash out all analysis except ### SETUP SCRIPTS TO GENERATE AND SET DIRECTORIES AND VARIABLES FOR FILES AND SCRIPTS ###
#### STEP 5 type ./Launcher_RNAseq.sh at the command line
#0# User will be updated with the path to fastqs in all scripts
#1# Fastq files will be moved into folders named by sample
#2# subdirectories will be created within each sample folder to hold analysis outputs and a logs folder will be created to hold error and output messages
#3# a samples.txt file will be created in the scripts directory which will be read to assign filenames to samples
#4# check that the user filepath is correct in scripts, check that the samples.txt filenames are unique and correct, check that directories have been created using ls */*/*/* from the fastq containing directory
#5# Once all is correct type nano Launcher_RNAseq.sh and hash out the setup scripts job submission commands, unhash the analysis sections you want to complete (prep (fastqc is compulsory, default is to use untrimmed fastqs, but this can be changed based on fastqc results), GATK, de (differential expression)) 
#4# type ./Launcher_RNAseq.sh and array jobs will be submit to the SUN GRID ENGINE for all samples in the samples.txt file. 
