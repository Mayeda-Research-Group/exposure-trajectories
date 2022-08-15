# exposure-trajectories

This repo contains all the code necessary to replicate optimization routines and analyses in the paper

Comparison of imputation strategies for incomplete longitudinal data in lifecourse epidemiology (under review)

Dataset: Health and Retirement Study (HRS waves 1998-2008; 2018) <br>
Exposure: Center for Epidemiological Studies-Depression (CES-D), measured longitudinally <br>
Outcome: Mortality by 2018 <br>

Methods: Multiple imputation (MI) using the MICE R package specifying three different MI methods (Joint Multivariate Normal Modeling, Predictive Mean Matching, Fully Conditional Specification) and a case study on MI using Linear Mixed Models

### Directory Descriptions
**Prelim Analyses**: preliminary analyses that helped us identify a useful exposure-outcome example in HRS for this project <br>

**RScripts**: all the RScripts required to clean data, run optimization routines, mask and impute data, and summarize and visualize results. RScripts are named in the order they should be run <br> 
&nbsp;&nbsp;&nbsp;&nbsp;A_find_subset.R identifies the optimal complete subset of 1998-2008 HRS for analysis <br>
&nbsp;&nbsp;&nbsp;&nbsp;The remaining scripts are run in numerical order <br>

&nbsp;&nbsp;&nbsp;&nbsp;**RScripts/OLD**: saving old code just in case it's needed-- will delete once paper is published <br>

&nbsp;&nbsp;&nbsp;&nbsp;**RScripts/Preliminary analysis**: additional analyses that helped us determine how to operationalize exposures, etc. <br>

&nbsp;&nbsp;&nbsp;&nbsp;**RScripts/Troubleshooting**: digging deeper when we had issues identifying problems with modeling missingness, etc. <br>
&nbsp;&nbsp;&nbsp;&nbsp;These issues have been resolved and this folder serves as helpful reminders for future projects <br>

&nbsp;&nbsp;&nbsp;&nbsp;**RScripts/functions**: custom functions that are called into analysis scripts <br>

**submission_scripts**: Analyses for this project were run on UCLA's Hoffman computing cluster (https://qcb.ucla.edu/collaboratory/hoffman2-cluster-user-guide/). These are submission scripts for those jobs.

### RScript Descriptions
All RScripts have descriptions of what they do at the top.
