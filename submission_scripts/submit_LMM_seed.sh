#!/bin/bash
#$ -cwd #uses current working directory
# error = Merged with joblog
#$ -o joblogs/joblog.$JOB_ID.$TASK_ID #creates a file called joblog.jobidnumber to write to. 
#$ -j y 
#$ -l highp,h_rt=30:00:00,h_data=4G #requests 30 hours, 4GB of data (per core) on biostat high performance nodes
#$ -pe shared 1 #requests 1 core
# Email address to notify
#$ -M $USER@mail #don't change this line, finds your email in the system 
# Notify when
#$ -m bea #sends you an email (b) when the job begins (e) when job ends (a) when job is aborted (error)
# submit array job:
# TO TEST MATTER RUN ONLY TWICE:
##$ -t 1-2:1
# FOR THE FULL RUN USE INSTEAD:
#$ -t 1-998:1
## 

# load the job environment:
. /u/local/Modules/default/init/modules.sh
module load R/4.0.2 #loads R/4.0.2 for use 
export OMP_NUM_THREADS=1 #uses max 1 thread (needs to match -pe shared)
## 
# run R code

#The LMM job will just be
#
#method = "LMM"
#mechanism = "MAR"
#mask_percent = "30%"

echo "======"
echo SGE_TASK_ID=$SGE_TASK_ID      
R CMD BATCH --no-save --no-restore '--args mechanism="MAR" method="LMM" mask_percent="30%" save="no" sens="no" '  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID
echo R CMD BATCH --no-save --no-restore \'--args mechanism=\"MAR\" method=\"LMM\" mask_percent=\"30%\" seed=123 save=\"no\" sens=\"no\" \'  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID

