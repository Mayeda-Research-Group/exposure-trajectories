#!/bin/bash
#$ -cwd #uses current working directory
# error = Merged with joblog
#$ -o joblog.$JOB_ID #creates a file called joblog.jobidnumber to write to. 
#$ -j y 
#$ -l h_rt=0:30:00,h_data=2G #requests 30 minutes, 2GB of data (per core)
#$ -pe shared 1 #requests 2 cores
# Email address to notify
#$ -M $USER@mail #don't change this line, finds your email in the system 
# Notify when
#$ -m bea #sends you an email (b) when the job begins (e) when job ends (a) when job is aborted (error)
## 
# load the job environment:
. /u/local/Modules/default/init/modules.sh
module load R/4.0.2 #loads R/4.0.2 for use 
## 
# run R code
echo 'Running mask_impute_pool.R' #prints this quote to joblog.jobidnumberR 
R CMD BATCH --no-save --no-restore '--args mechanism = "MNAR" method = "JMVN" mask_percent = "10%'' '  mask_impute_pool.R output.$JOB_ID