#!/bin/bash
#$ -cwd #uses current working directory
# error = Merged with joblog
#$ -o joblogs/joblog.$JOB_ID.$TASK_ID #creates a file called joblog.jobidnumber to write to. 
#$ -j y 
#$ -l highp,h_rt=60:00:00,h_data=4G #requests 60 hours, 4GB of data (per core) on biostat high performance nodes
#$ -pe shared 1 #requests 1 core
# Email address to notify
#$ -M $USER@mail #don't change this line, finds your email in the system 
# Notify when
#$ -m bea #sends you an email (b) when the job begins (e) when job ends (a) when job is aborted (error)
# submit array job:
# TO TEST MATTER RUN ONLY TWICE:
##$ -t 1-2:1
# FOR THE FULL RUN USE INSTEAD (length of the array of missing jobs):
#$ -t 1-28:1
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
#mechanism = "MNAR"
#mask_percent = "30%"

echo "======"
echo SGE_TASK_ID=$SGE_TASK_ID 

# SET THE ARRAY WITH THE $SGE_TASK_IDs that failed 
### (note: UPDATE line "#$ -t 1-1000:1" ACCORDINGLY!!!) 

array=(949 950 951 952 975 978 979 980 981 982 983 984 985 986 987 988 989 990 991 992 993 994 995 996 997 998 999 1000)

# SET THE SEED
###(note: array index starts at 0 so we need to scale back the current $SGE_TASK_ID)
myindex=${array[$SGE_TASK_ID-1]}

R CMD BATCH --no-save --no-restore "--args mechanism=\"MNAR\" method=\"LMM\" mask_percent=\"30%\" seed=$myindex save=\"no\" sens=\"no\" "  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID
echo R CMD BATCH --no-save --no-restore \'--args mechanism=\"MNAR\" method=\"LMM\" mask_percent=\"30%\" seed=$myindex save=\"no\" sens=\"no\" \'  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID

