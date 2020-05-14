#Start with a simple model CysC ~ age + sex/gender + race/ethnicity
#Figure out what the right SD is to add to the predicted values from the lmm

#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "lme4", "tidyverse", "magrittr")

#---- Read in data ----
analytic_df <- read_csv(here::here("Data", "analytic_df.csv"))

#---- wide --> long ----
long_df <- analytic_df %>% 
  dplyr::select("HHIDPN", "female", "hispanic", "black", "other", 
                contains("AGE", ignore.case = FALSE), contains("CYSC_ADJ")) %>% 
  pivot_longer(cols = c(contains("AGE"), contains("CYSC_ADJ")), 
               names_to = c("wave", ".value"), 
               names_pattern = "(.)(.)") %>% 
  set_colnames(c("HHIDPN", "female", "hispanic", "black", "other", "wave", 
                 "Age", "CYSC"))

#---- look at the outcome variable ----
hist(long_df$CYSC)
boxplot(long_df$CYSC)

#CYSC is really skewed-- take the log of the outcome
long_df[, "log_CSYC"] <- log(long_df$CYSC)

#Look at transformed outcome-- super symmetric!
hist(long_df$log_CSYC)
boxplot(long_df$log_CSYC)

#---- Model ----

