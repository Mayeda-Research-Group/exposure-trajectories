#---- Package loading + seed setting ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr")

set.seed(20200517)

#---- Load data ----
#This dataset has the HRS data joined with the biomarker data
analytic_df <- read_csv(here::here("Data", "analytic_df.csv"))

#---- Wide --> long ----
long_df <- analytic_df %>% 
  dplyr::select("HHIDPN", "female", "hispanic", "black", "other", 
                contains("AGE", ignore.case = FALSE), contains("CYSC_ADJ")) %>% 
  pivot_longer(cols = c(contains("AGE"), contains("CYSC_ADJ")), 
               names_to = c("wave", ".value"), 
               names_pattern = "(.)(.)") %>% 
  set_colnames(c("HHIDPN", "female", "hispanic", "black", "other", "wave", 
                 "Age", "CYSC"))

long_df[, "log_CYSC"] <- log(long_df$CYSC)

#---- Complete data only ----
complete_data <- long_df %>% na.omit()

#---- Induce missingness ----
degree_of_missingness <- c(0.25)

#MCAR
for(prop in degree_of_missingness){
  test_ind <- sample(seq(1, nrow(complete_data)), 
                     floor(prop*nrow(complete_data)))
  assign(paste0("MCAR_", prop*100, "_train"), complete_data[-test_ind, ])
  assign(paste0("MCAR_", prop*100, "_test"), complete_data[test_ind, ])
  assign(paste0("MCAR_", prop*100, "_test"), 
         `[[<-`(get(paste0("MCAR_", prop*100, "_test")), 'Missingness', 
                value = paste0(prop*100, "% Missingness")))}