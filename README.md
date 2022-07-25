# exposure-trajectories

This repo contains all the code necessary to replicate optimization routines and analyses in the paper

Comparison of imputation strategies for incomplete longitudinal data in lifecourse epidemiology

Dataset: HRS (1998-2008; 2018) <br>
Exposure: Center for Epidemiological Studies-Depression (CES-D) measured longitudinally <br>
Outcome: Mortality by 2018 <br>

### Directories
RScripts: all the RScripts required to clean data, run optimization routines, mask and impute data, and summarize and visualize results. RScripts are named in the order they should be run. <br> 
&nbsp;&nbsp;&nbsp;&nbsp;A_find_subset.R identifies the optimal complete subset of 1998-2008 HRS for analysis <br>
&nbsp;&nbsp;&nbsp;&nbsp;The remaining scripts are run in numerical order <br>

RScripts/functions: custom functions that are called into analysis scripts <br>

submission_scripts: Analyses for this project were run on UCLA's Hoffman computing cluster (https://qcb.ucla.edu/collaboratory/hoffman2-cluster-user-guide/). These are submission scripts for those jobs.

All RScripts have descriptions of what they do at the top.
