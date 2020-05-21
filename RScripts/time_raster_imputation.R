#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr")

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

#---- Induce MCAR missingness ----
degree_of_missingness <- c(0.25)

for(prop in degree_of_missingness){
  test_ind <- sample(seq(1, nrow(complete_data)), 
                     floor(prop*nrow(complete_data)))
  assign(paste0("MCAR_", prop*100, "_train"), complete_data[-test_ind, ])
  assign(paste0("MCAR_", prop*100, "_test"), complete_data[test_ind, ])
  assign(paste0("MCAR_", prop*100, "_test"), 
         `[[<-`(get(paste0("MCAR_", prop*100, "_test")), 'Missingness', 
                value = paste0(prop*100, "% Missingness")))
}

#Make sure that individuals have at least one observation in training set
slot = 1
df_list <- list()
for(prop in degree_of_missingness){
  df_list[[slot]] <- 
    eval(parse(text = paste0("MCAR_", prop*100, "_test")))[which(
      eval(parse(text = paste0("MCAR_", prop*100, "_test")))$HHIDPN %in% 
        eval(parse(text = paste0("MCAR_", prop*100, "_train")))$HHIDPN), ]
  slot = slot + 1
}


