#!/bin/bash
#$ -cwd #uses current working directory
# error = Merged with joblog
#$ -o joblogs/joblog.$JOB_ID.$TASK_ID #creates a file called joblog.jobidnumber to write to. 
#$ -j y 
#$ -l h_rt=24:00:00,h_data=8G #requests 24 hours, 8GB of data (per core)
#$ -pe shared 8 #requests 8 cores
# Email address to notify
#$ -M $USER@mail #don't change this line, finds your email in the system 
# Notify when
##$ -m n
#$ -m bea #sends you an email (b) when the job begins (e) when job ends (a) when job is aborted (error)
# submit array job:
# TO TEST MATTER RUN ONLY TWICE:
#$ -t 1-2:1
# FOR THE FULL RUN USE INSTEAD:
##$ -t 1-100:1
## 

# load the job environment:
. /u/local/Modules/default/init/modules.sh
module load R/4.0.2 #loads R/4.0.2 for use 
export OMP_NUM_THREADS=8 #uses max 8 threads (needs to match -pe shared)
## 
# run R code

#The LMM job will just be
#
#method = "LMM"
#mechanism = "MAR"
#mask_percent = "30%"

echo "======"
echo SGE_TASK_ID=$SGE_TASK_ID      
R CMD BATCH --no-save --no-restore '--args mechanism="MAR" method="LMM" mask_percent="30%" '  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID
echo R CMD BATCH --no-save --no-restore \'--args mechanism=\"MAR\" method=\"LMM\" mask_percent=\"30%\" \'  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID
