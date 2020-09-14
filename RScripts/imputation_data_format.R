#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr")

#No scientific notation
options(scipen = 999)

#---- Read in analytical sample ----
analytical_sample <- read_csv(paste0("/Users/CrystalShaw/Dropbox/Projects/", 
                                     "exposure_trajectories/data/", 
                                     "hrs_samp_alive_70_cysc_60_70.csv"))

#---- Wave indicators ----
letter_waves <- LETTERS[seq(from = 11, to = 15)] #biomarker sample

#---- Select variables ----
keep <- c("HHIDPN", paste0(head(letter_waves), "age_y_int"), "female", 
          "hispanic", "black", "other", "raedyrs", "cses_index", "death", 
          "age_death_y", "smoker", "A1C_ADJ", "TC_ADJ", "HDL_ADJ", "height", 
          "weight", "BMI", "bpsys", "bpdia", "CYSC_ADJ")

impute <- analytical_sample %>% dplyr::select(contains(keep))

#---- Age range ----
age_range <- c(min(impute %>% dplyr::select(contains("age_y"))), 
               max(impute %>% dplyr::select(contains("age_y"))))

# #Checking if any of the in between ages are missing (none)
# #There are no observations at age 60
# sum(!which(seq(age_range[1], age_range[2]) %in% 
#              (impute %>% dplyr::select(contains("age_y")))))

#---- long data ----
impute_long <- impute %>% 
  pivot_longer(cols = contains(c("age_y", "A1C_ADJ", "TC_ADJ", "HDL_ADJ", 
                                 "weight", "BMI", "bpsys", "bpdia", "CYSC_ADJ"), 
                               ignore.case = FALSE), 
               names_to = c("Wave", ".value"),
               names_pattern = "(.)(.*)") %>% 
  arrange(HHIDPN) 

#---- Add rows for age bracket of interest ----
#Exposure ages: [60, 69]
invariant <- c("female", "hispanic", "black", "other", "raedyrs", "cses_index", 
               "death", "age_death_y", "smoker", "height_measured", "height")
measured_indicators <- c("weight_measured", "BMI_measured")

exp_ages <- seq(60, 69, by = 1)

HHIDPNs <- unique(impute_long$HHIDPN)
for(ID in HHIDPNs){
  subset <- impute_long %>% dplyr::filter(HHIDPN == ID)
  obs_ages <- subset$age_y_int
  missing_ages <- exp_ages[which(!exp_ages %in% obs_ages)]
  
  #Create new rows
  new_rows <- as.data.frame(matrix(nrow = length(missing_ages), 
                                   ncol = ncol(impute_long))) %>% 
    set_colnames(colnames(impute_long))
  
  new_rows$HHIDPN <- ID
  new_rows$age_y_int <- missing_ages
  new_rows[, invariant] <- subset[1, invariant]
  new_rows[, measured_indicators] <- 0
  
  #Add new rows to long data
  impute_long %<>% rbind(new_rows)
}

# #Sanity Check-- num measures at each age (only 60 should be 0)
# measures_by_age <- impute_long %>% dplyr::group_by(age_y_int) %>%
#   summarize_at("CYSC_ADJ", ~sum(!is.na(.)))
# 
# #If measures are not observed, they should be set to 0, not NA
# sum(is.na(impute_long$weight_measured))
# sum(is.na(impute_long$BMI_measured))

#---- log CysC ----
impute_long %<>% mutate("log_CysC" = log(CYSC_ADJ)) 

# #Look at the distributions of CysC and log(CysC)
# ggplot(impute_long %>%
#          dplyr::select("CYSC_ADJ", "age_y_int"), aes(x = CYSC_ADJ)) +
#   geom_histogram() + facet_wrap(~ age_y_int)
# 
# ggplot(impute_long %>%
#          dplyr::select("log_CysC", "age_y_int"), aes(x = log_CysC)) +
#   geom_histogram() + facet_wrap(~ age_y_int)

#---- save datasets ----
write_csv(impute_long, paste0("/Users/CrystalShaw/Dropbox/Projects/",
                              "exposure_trajectories/data/",
                              "imputation_data_long.csv"))




