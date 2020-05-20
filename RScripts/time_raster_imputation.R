#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here")

set.seed(20200520)

#---- Read in and format data ----
analytic_df <- read_csv(here::here("Data", "analytic_df.csv"))

#wide --> long
long_df <- analytic_df %>% 
  dplyr::select("HHIDPN", "female", "hispanic", "black", "other", 
                contains("AGE", ignore.case = FALSE), contains("CYSC_ADJ")) %>% 
  pivot_longer(cols = c(contains("AGE"), contains("CYSC_ADJ")), 
               names_to = c("wave", ".value"), 
               names_pattern = "(.)(.)") %>% 
  set_colnames(c("HHIDPN", "female", "hispanic", "black", "other", "wave", 
                 "Age", "CYSC"))

long_df[, "log_CYSC"] <- log(long_df$CYSC)

#Get rid of all of the missing observations 
#We have 15,653 people with at least one measure
complete_data <- long_df %>% na.omit()
#length(unique(complete_data$HHIDPN))

#

