################################################################################
### this submission script was created on NOV2021
### with ref: https://support.idre.ucla.edu/helpdesk/Ticket/40244340
### author RD for inquiries: hpc@ucla.edu
### it submits an array job of 900 tasks for:
### CC jobs are all combos of these variables:
### method = "CC"
### mechanism = c("MCAR", "MAR", "MNAR")
### mask_percent = c("10%", "20%", "30%")
### 9 sets of parameters that each need to run 100 times (total 900 jobs)
################################################################################
#!/bin/bash
#$ -cwd #uses current working directory
# error = Merged with joblog
#$ -o joblog.$JOB_ID.$TASK_ID #creates a file called joblog.jobidnumber to write to. 
#$ -j y 
#$ -l h_rt=2:00:00,h_data=2G #requests 2 hours, 2GB of data (per core)
#$ -pe shared 1 #requests 2 cores
# Email address to notify
#$ -M $USER@mail #don't change this line, finds your email in the system 
# Notify when
##$ -m n
#$ -m bea #sends you an email (b) when the job begins (e) when job ends (a) when job is aborted (error)
# submit array job:
# TEST OUT BY RUNNING ONLY 9 CASES:
#$ -t 1-9:1
# FOR THE FULL RUN USE INSTEAD:
##$ -t 1-900:1
## 

# load the job environment:
. /u/local/Modules/default/init/modules.sh
module load R/4.0.2 #loads R/4.0.2 for use 
## 
# run R code

#The sets of CC jobs are all combos of these variables:
#
#method = "CC"
#mechanism = c("MCAR", "MAR", "MNAR")
#mask_percent = c("10%", "20%", "30%")
#
#So this makes 9 sets of parameters that each need to run 100 times (total 900 jobs) [is my math correct lol]
#
#The LMM job will just be
#
#method = "LMM"
#mechanism = "MAR"
#mask_percent = "30%"


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
        echo R CMD BATCH --no-save --no-restore \'--args mechanism=\"MCAR\" method=\"CC\" mask_percent=\"30%\" \'  mask_impute_pool.R output.$JOB_ID.$SGE_TASK_ID
        R CMD BATCH --no-save --no-restore '--args mechanism="MCAR" method="CC" mask_percent="30%" '  mask_impute_pool.R output.$JOB_ID.$SGE_TASK_ID
    elif [ $mylastindex == 1 ]; then
        n=$mylastindex
        ##mask_percent="10%"
        echo SGE_TASK_ID=$SGE_TASK_ID, MYINDEX=$myindex, MYLASTINDEX=$mylastindex, l=$myindex, m=$m, n=$n      
        echo R CMD BATCH --no-save --no-restore \'--args mechanism=\"MCAR\" method=\"CC\" mask_percent=\"10%\" \'  mask_impute_pool.R output.$JOB_ID.$SGE_TASK_ID
        R CMD BATCH --no-save --no-restore '--args mechanism="MCAR" method="CC" mask_percent="10%" '  mask_impute_pool.R output.$JOB_ID.$SGE_TASK_ID
    elif [ $mylastindex == 2 ]; then
        n=$mylastindex
        ##mask_percent="20%"
        echo SGE_TASK_ID=$SGE_TASK_ID, MYINDEX=$myindex, MYLASTINDEX=$mylastindex, l=$myindex, m=$m, n=$n      
        echo R CMD BATCH --no-save --no-restore \'--args mechanism=\"MCAR\" method=\"CC\" mask_percent=\"20%\" \'  mask_impute_pool.R output.$JOB_ID.$SGE_TASK_ID
        R CMD BATCH --no-save --no-restore '--args mechanism="MCAR" method="CC" mask_percent="20%" '  mask_impute_pool.R output.$JOB_ID.$SGE_TASK_ID
    fi   


elif [ $myindex -ge 4 ] && [ $myindex -lt 7 ] ; then


    m=2
    ##mechanism="MAR"

    if [ $mylastindex == 0 ]; then
        n=3
        ##mask_percent="30%"
        echo SGE_TASK_ID=$SGE_TASK_ID, MYINDEX=$myindex, MYLASTINDEX=$mylastindex, l=$myindex, m=$m, n=$n      
        echo R CMD BATCH --no-save --no-restore \'--args mechanism=\"MAR\" method=\"CC\" mask_percent=\"30%\" \'  mask_impute_pool.R output.$JOB_ID.$SGE_TASK_ID
        R CMD BATCH --no-save --no-restore '--args mechanism="MAR" method="CC" mask_percent="30%" '  mask_impute_pool.R output.$JOB_ID.$SGE_TASK_ID
    elif [ $mylastindex == 1 ]; then
        n=$mylastindex
        ##mask_percent="10%"
        echo SGE_TASK_ID=$SGE_TASK_ID, MYINDEX=$myindex, MYLASTINDEX=$mylastindex, l=$myindex, m=$m, n=$n      
        echo R CMD BATCH --no-save --no-restore \'--args mechanism=\"MAR\" method=\"CC\" mask_percent=\"10%\" \'  mask_impute_pool.R output.$JOB_ID.$SGE_TASK_ID
        R CMD BATCH --no-save --no-restore '--args mechanism="MAR" method="CC" mask_percent="10%" '  mask_impute_pool.R output.$JOB_ID.$SGE_TASK_ID
    elif [ $mylastindex == 2 ]; then
        n=$mylastindex
        ##mask_percent="20%"
        echo SGE_TASK_ID=$SGE_TASK_ID, MYINDEX=$myindex, MYLASTINDEX=$mylastindex, l=$myindex, m=$m, n=$n      
        echo R CMD BATCH --no-save --no-restore \'--args mechanism=\"MAR\" method=\"CC\" mask_percent=\"20%\" \'  mask_impute_pool.R output.$JOB_ID.$SGE_TASK_ID
        R CMD BATCH --no-save --no-restore '--args mechanism="MAR" method="CC" mask_percent="20%" '  mask_impute_pool.R output.$JOB_ID.$SGE_TASK_ID
    fi   


elif [ $myindex -ge 7 ] &&[ $myindex -le 9 ]; then


    m=3
    ##mechanism="MNAR"
    
    if [ $mylastindex == 0 ]; then
        n=3
        ##mask_percent="30%"
        echo SGE_TASK_ID=$SGE_TASK_ID, MYINDEX=$myindex, MYLASTINDEX=$mylastindex, l=$myindex, m=$m, n=$n      
        echo R CMD BATCH --no-save --no-restore \'--args mechanism=\"MNAR\" method=\"CC\" mask_percent=\"30%\" \'  mask_impute_pool.R output.$JOB_ID.$SGE_TASK_ID
        R CMD BATCH --no-save --no-restore '--args mechanism="MNAR" method="CC" mask_percent="30%" '  mask_impute_pool.R output.$JOB_ID.$SGE_TASK_ID
    elif [ $mylastindex == 1 ]; then
        n=$mylastindex
        ##mask_percent="10%"
        echo SGE_TASK_ID=$SGE_TASK_ID, MYINDEX=$myindex, MYLASTINDEX=$mylastindex, l=$myindex, m=$m, n=$n      
        echo R CMD BATCH --no-save --no-restore \'--args mechanism=\"MNAR\" method=\"CC\" mask_percent=\"10%\" \'  mask_impute_pool.R output.$JOB_ID.$SGE_TASK_ID
        R CMD BATCH --no-save --no-restore '--args mechanism="MNAR" method="CC" mask_percent="10%" '  mask_impute_pool.R output.$JOB_ID.$SGE_TASK_ID
    elif [ $mylastindex == 2 ]; then
        n=$mylastindex
        ##mask_percent="20%"
        echo SGE_TASK_ID=$SGE_TASK_ID, MYINDEX=$myindex, MYLASTINDEX=$mylastindex, l=$myindex, m=$m, n=$n      
        echo R CMD BATCH --no-save --no-restore \'--args mechanism=\"MNAR\" method=\"CC\" mask_percent=\"20%\" \'  mask_impute_pool.R output.$JOB_ID.$SGE_TASK_ID     
        R CMD BATCH --no-save --no-restore '--args mechanism="MNAR" method="CC" mask_percent="20%" '  mask_impute_pool.R output.$JOB_ID.$SGE_TASK_ID
    fi          

fi

echo "TASK $SGE_TASK_ID of JOB $JOB_ID ended at " `date`
echo ""

