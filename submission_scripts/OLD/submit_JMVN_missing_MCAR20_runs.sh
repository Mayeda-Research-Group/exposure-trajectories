#!/bin/bash
#$ -cwd #uses current working directory
# error = Merged with joblog
#$ -o joblogs/joblog.$JOB_ID.$TASK_ID #creates a file called joblog.jobidnumber to write to. 
#$ -j y 
#$ -l h_rt=4:00:00,h_data=2G #requests 4 hours, 2GB of data (per core)
#$ -pe shared 1 #requests 1 cores
# Email address to notify
#$ -M $USER@mail #don't change this line, finds your email in the system 
# Notify when
#$ -m bea #sends you an email (b) when the job begins (e) when job ends (a) when job is aborted (error)
# submit array job:
# FULL RUN:
#$ -t 1-1:1
## 

# load the job environment:
. /u/local/Modules/default/init/modules.sh
module load R/4.0.2 #loads R/4.0.2 for use 
export OMP_NUM_THREADS=1 #uses max 1 threads (needs to match -pe shared)
## 
# run R code

#The JMVN job will just be
#
#method = "JMVN"
#mechanism = "MCAR"
#mask_percent = "20%"

echo "======"
echo SGE_TASK_ID=$SGE_TASK_ID      
R CMD BATCH --no-save --no-restore '--args mechanism="MCAR" method="JMVN" mask_percent="20%" save="no" sens="no" '  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID
echo R CMD BATCH --no-save --no-restore \'--args mechanism=\"MCAR\" method=\"JMVN\" mask_percent=\"20%\" save=\"no\" sens=\"no\" \'  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID

