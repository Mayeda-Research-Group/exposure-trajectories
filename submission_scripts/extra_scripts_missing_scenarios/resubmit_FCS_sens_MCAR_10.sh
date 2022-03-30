###########################################################
# USE THIS SCRIPT TO RESUBMIT TASKS THAT DID NOT COMPLETE #
# DATE/AUTHOR: 2022MAR29RD RD (hpc@ucla.edu)              #
# LAST UPDATED: 2022MAR29RD (hpc@ucla.edu)                #
###########################################################
#!/bin/bash
#$ -cwd #uses current working directory
# error = Merged with joblog
#$ -o joblogs/joblog.$JOB_ID.$TASK_ID #creates a file called joblog.jobidnumber to write to. 
#$ -j y 
#$ -l h_rt=18:00:00,h_data=3G #requests 18 hours, 3GB of data (per core)
#$ -pe shared 1 #requests 1 cores
# Email address to notify
#$ -M $USER@mail #don"t change this line, finds your email in the system 
# Notify when
#$ -m bea #sends you an email (b) when the job begins (e) when job ends (a) when job is aborted (error)
# Run a job array on the cases that did not complete (need to save the task id that failed in an array)
#$ -t 1-30:1
# testing on 3 cases:
##$ -t 1-3:1

# load the job environment:
. /u/local/Modules/default/init/modules.sh
module load R/4.0.2 #loads R/4.0.2 for use 
export OMP_NUM_THREADS=1 #uses max 1 threads (needs to match -pe shared)
echo ""
module li
echo OMP_NUM_THREADS=$OMP_NUM_THREADS
echo ""

echo "TASK $SGE_TASK_ID of JOB $JOB_ID run on " `hostname`
echo "TASK $SGE_TASK_ID of JOB $JOB_ID ended at " `date`
echo ""

# SET THE ARRAY WITH THE $SGE_TASK_IDs that failed 
### (note: UPDATE line "#$ -t 1-3:1" ACCORDINGLY!!!) 

array=(1351 2260 2323 2710 2917 3142 3178 3682 3781 4123 4537 4897 4978 5077 5095 5248 5365 5869 5923 6166 6598 6670 6751 6760 6832 6967 7147 7183 7516 8056)

# SET THE SEED
###(note: array index starts at 0 so we need to scale back the current $SGE_TASK_ID)

seed=${array[$SGE_TASK_ID-1]}


# run R code
echo "R CMD BATCH --no-save --no-restore "--args mechanism=\"MCAR\" method=\"FCS\" mask_percent=\"10%\" seed=$seed save=\"no\" sens=\"yes\" "  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID"
R CMD BATCH --no-save --no-restore "--args mechanism=\"MCAR\" method=\"FCS\" mask_percent=\"10%\" seed=$seed save=\"no\" sens=\"yes\" "  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID

echo ""
echo "TASK $SGE_TASK_ID of JOB $JOB_ID run on " `hostname`
echo "TASK $SGE_TASK_ID of JOB $JOB_ID ended at " `date`
echo ""
