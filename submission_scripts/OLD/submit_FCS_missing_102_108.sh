################################################################################
### this submission script was created on JAN2022
### with ref: https://support.idre.ucla.edu/helpdesk/Ticket/40244340
### author RD for inquiries: hpc@ucla.edu
### it submits an array job of 9000 tasks for:
### FCS jobs are all combos of these variables:
### method = "FCS"
### mechanism = c("MCAR", "MAR", "MNAR")
### mask_percent = c("10%", "20%", "30%")
### 9 sets of parameters that each need to run 1000 times (total 9000 jobs)
################################################################################
#!/bin/bash
#$ -cwd #uses current working directory
# error = Merged with joblog
#$ -o joblogs/joblog.$JOB_ID.$TASK_ID #creates a file called joblog.jobidnumber to write to. 
#$ -j y 
#$ -l h_rt=16:00:00,h_data=3G #requests 16 hours, 3GB of data (per core)
#$ -pe shared 1 #requests 1 cores
# Email address to notify
#$ -M $USER@mail #don"t change this line, finds your email in the system 
# Notify when
#$ -m bea #sends you an email (b) when the job begins (e) when job ends (a) when job is aborted (error)
# submit array job:
# TEST OUT BY RUNNING ONLY 9 CASES:
##$ -t 1-9:1
# FOR THE FULL RUN USE INSTEAD (9*num runs):
#$ -t 102-108:3
## 

# load the job environment:
. /u/local/Modules/default/init/modules.sh
module load R/4.0.2 #loads R/4.0.2 for use 
export OMP_NUM_THREADS=1 #uses max 1 threads (needs to match -pe shared)
## 
# run R code

#The sets of FCS jobs are all combos of these variables:
#
#method = "FCS"
#mechanism = c("MCAR", "MAR", "MNAR")
#mask_percent = c("10%", "20%", "30%")
#
#So this makes 9 sets of parameters that each need to run 1000 times (total 9000 jobs)

echo "======"
echo "TASK $SGE_TASK_ID of JOB $JOB_ID started on " `hostname`
echo "TASK $SGE_TASK_ID of JOB $JOB_ID started at " `date`

myindex=$((SGE_TASK_ID%9))
mylastindex=$((myindex%3))
if [ $myindex == 0 ]; then
  myindex=9
fi

if [ $myindex -gt 0 ] && [ $myindex -lt 4 ]; then

    m=1
    ##mechanism="MCAR"

    if [ $mylastindex == 0 ]; then
        n=3
        ##mask_percent="30%"
        echo SGE_TASK_ID=$SGE_TASK_ID, MYINDEX=$myindex, MYLASTINDEX=$mylastindex, l=$myindex, m=$m, n=$n      
        echo R CMD BATCH --no-save --no-restore \"--args mechanism=\"MCAR\" method=\"FCS\" mask_percent=\"30%\" seed=$SGE_TASK_ID save=\"no\" sens=\"no\" \"  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID
        R CMD BATCH --no-save --no-restore "--args mechanism=\"MCAR\" method=\"FCS\" mask_percent=\"30%\" seed=$SGE_TASK_ID save=\"no\" sens=\"no\" "  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID
    elif [ $mylastindex == 1 ]; then
        n=$mylastindex
        ##mask_percent="10%"
        echo SGE_TASK_ID=$SGE_TASK_ID, MYINDEX=$myindex, MYLASTINDEX=$mylastindex, l=$myindex, m=$m, n=$n      
        echo R CMD BATCH --no-save --no-restore \"--args mechanism=\"MCAR\" method=\"FCS\" mask_percent=\"10%\" seed=$SGE_TASK_ID save=\"no\" sens=\"no\" \"  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID
        R CMD BATCH --no-save --no-restore "--args mechanism=\"MCAR\" method=\"FCS\" mask_percent=\"10%\" seed=$SGE_TASK_ID save=\"no\" sens=\"no\" "  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID
    elif [ $mylastindex == 2 ]; then
        n=$mylastindex
        ##mask_percent="20%"
        echo SGE_TASK_ID=$SGE_TASK_ID, MYINDEX=$myindex, MYLASTINDEX=$mylastindex, l=$myindex, m=$m, n=$n      
        echo R CMD BATCH --no-save --no-restore \"--args mechanism=\"MCAR\" method=\"FCS\" mask_percent=\"20%\" seed=SGE_TASK_ID save=\"no\" sens=\"no\" \"  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID
        R CMD BATCH --no-save --no-restore "--args mechanism=\"MCAR\" method=\"FCS\" mask_percent=\"20%\" seed=$SGE_TASK_ID save=\"no\" sens=\"no\" "  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID
    fi   


elif [ $myindex -ge 4 ] && [ $myindex -lt 7 ] ; then


    m=2
    ##mechanism="MAR"

    if [ $mylastindex == 0 ]; then
        n=3
        ##mask_percent="30%"
        echo SGE_TASK_ID=$SGE_TASK_ID, MYINDEX=$myindex, MYLASTINDEX=$mylastindex, l=$myindex, m=$m, n=$n      
        echo R CMD BATCH --no-save --no-restore \"--args mechanism=\"MAR\" method=\"FCS\" mask_percent=\"30%\" seed=$SGE_TASK_ID save=\"no\" sens=\"no\" \"  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID
        R CMD BATCH --no-save --no-restore "--args mechanism=\"MAR\" method=\"FCS\" mask_percent=\"30%\" seed=$SGE_TASK_ID save=\"no\" sens=\"no\" "  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID
    elif [ $mylastindex == 1 ]; then
        n=$mylastindex
        ##mask_percent="10%"
        echo SGE_TASK_ID=$SGE_TASK_ID, MYINDEX=$myindex, MYLASTINDEX=$mylastindex, l=$myindex, m=$m, n=$n      
        echo R CMD BATCH --no-save --no-restore \"--args mechanism=\"MAR\" method=\"FCS\" mask_percent=\"10%\" seed=$SGE_TASK_ID save=\"no\" sens=\"no\" \"  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID
        R CMD BATCH --no-save --no-restore "--args mechanism=\"MAR\" method=\"FCS\" mask_percent=\"10%\" seed=$SGE_TASK_ID save=\"no\" sens=\"no\" "  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID
    elif [ $mylastindex == 2 ]; then
        n=$mylastindex
        ##mask_percent="20%"
        echo SGE_TASK_ID=$SGE_TASK_ID, MYINDEX=$myindex, MYLASTINDEX=$mylastindex, l=$myindex, m=$m, n=$n      
        echo R CMD BATCH --no-save --no-restore \"--args mechanism=\"MAR\" method=\"FCS\" mask_percent=\"20%\" seed=$SGE_TASK_ID save=\"no\" sens=\"no\" \"  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID
        R CMD BATCH --no-save --no-restore "--args mechanism=\"MAR\" method=\"FCS\" mask_percent=\"20%\" seed=$SGE_TASK_ID save=\"no\" sens=\"no\" "  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID
    fi   


elif [ $myindex -ge 7 ] &&[ $myindex -le 9 ]; then


    m=3
    ##mechanism="MNAR"
    
    if [ $mylastindex == 0 ]; then
        n=3
        ##mask_percent="30%"
        echo SGE_TASK_ID=$SGE_TASK_ID, MYINDEX=$myindex, MYLASTINDEX=$mylastindex, l=$myindex, m=$m, n=$n      
        echo R CMD BATCH --no-save --no-restore \"--args mechanism=\"MNAR\" method=\"FCS\" mask_percent=\"30%\" seed=$SGE_TASK_ID save=\"no\" sens=\"no\" \"  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID
        R CMD BATCH --no-save --no-restore "--args mechanism=\"MNAR\" method=\"FCS\" mask_percent=\"30%\" seed=$SGE_TASK_ID save=\"no\" sens=\"no\" "  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID
    elif [ $mylastindex == 1 ]; then
        n=$mylastindex
        ##mask_percent="10%"
        echo SGE_TASK_ID=$SGE_TASK_ID, MYINDEX=$myindex, MYLASTINDEX=$mylastindex, l=$myindex, m=$m, n=$n      
        echo R CMD BATCH --no-save --no-restore \"--args mechanism=\"MNAR\" method=\"FCS\" mask_percent=\"10%\" seed=$SGE_TASK_ID save=\"no\" sens=\"no\" \"  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID
        R CMD BATCH --no-save --no-restore "--args mechanism=\"MNAR\" method=\"FCS\" mask_percent=\"10%\" seed=$SGE_TASK_ID save=\"no\" sens=\"no\" "  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID
    elif [ $mylastindex == 2 ]; then
        n=$mylastindex
        ##mask_percent="20%"
        echo SGE_TASK_ID=$SGE_TASK_ID, MYINDEX=$myindex, MYLASTINDEX=$mylastindex, l=$myindex, m=$m, n=$n      
        echo R CMD BATCH --no-save --no-restore \"--args mechanism=\"MNAR\" method=\"FCS\" mask_percent=\"20%\" seed=$SGE_TASK_ID save=\"no\" sens=\"no\" \"  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID     
        R CMD BATCH --no-save --no-restore "--args mechanism=\"MNAR\" method=\"FCS\" mask_percent=\"20%\" seed=$SGE_TASK_ID save=\"no\" sens=\"no\" "  mask_impute_pool.R output/output.$JOB_ID.$SGE_TASK_ID
    fi          

fi

echo "TASK $SGE_TASK_ID of JOB $JOB_ID ended at " `date`
echo ""

