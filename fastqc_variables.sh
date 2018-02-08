#!/bin/bash 

WORK=/user/
scripts=/user/scripts
readarray -t FILENAME < ${scripts}/samples.txt 
fastqcOutput=$WORK/${FILENAME[$SGE_TASK_ID-1]}/processedData/'fastqc'

