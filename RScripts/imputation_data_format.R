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

#---- CysC by age ----
impute_long <- impute %>% 
  pivot_longer(cols = contains(c("age_y", "CYSC_ADJ"), ignore.case = FALSE), 
               names_to = c("Wave", ".value"),
               names_pattern = "(.)(.*)") %>% 
  arrange(age_y_int) %>% filter(!is.na(CYSC_ADJ))
  
#Youngest age is 61, so will need to add a column for age 60
age_range <- c(min(impute_long$age_y_int), max(impute_long$age_y_int))
#Checking if any of the in between ages are missing (none)
sum(!which(seq(age_range[1], age_range[2]) %in% impute_long$age_y_int))

impute <- impute_long %>%
  mutate("label" = "CYSC_ADJ") %>% 
  unite("names", c("label", "age_y_int"), sep = "_") %>% 
  #get rid of Wave variable
  dplyr::select(-Wave) %>%
  #get columns of CysC by age
  pivot_wider(names_from = "names", values_from = "CYSC_ADJ") %>% 
  mutate("CYSC_ADJ_60" = NA) %>% mutate_at("CYSC_ADJ_60", as.numeric) %>%
  relocate(CYSC_ADJ_60, .before = CYSC_ADJ_61) 

# #Sanity Check-- num measures at each age (only 60 should be 0)
# colSums(1 - is.na(impute %>% dplyr::select(contains("CYSC_ADJ"))))
  
#---- log CysC ----
log_CysC <- log(impute %>% dplyr::select(contains("CYSC"))) 
colnames(log_CysC) <- paste0("log", colnames(log_CysC))
impute %<>% cbind(log_CysC) 

# #Look at the distributions of CysC and log(CysC)
# ggplot(impute %>% 
#          dplyr::select(paste0("CYSC_ADJ_", seq(60, age_range[2]))) %>% 
#          pivot_longer(cols = everything(), 
#                       names_to = "Age", values_to = "CysC"), aes(x = CysC)) + 
#   geom_histogram() + facet_wrap(~ Age)
# 
# ggplot(impute %>% 
#          dplyr::select(paste0("logCYSC_ADJ_", seq(60, age_range[2]))) %>% 
#          pivot_longer(cols = everything(), 
#                       names_to = "Age", values_to = "CysC"), aes(x = CysC)) + 
#   geom_histogram() + facet_wrap(~ Age)

#---- remove extra variables ----
impute %<>% dplyr::select(-paste0("CYSC_ADJ_", seq(60, age_range[2])))

#---- create long data ----
impute_long <- impute %>% 
  pivot_longer(cols = contains(c("logCYSC_ADJ"), ignore.case = FALSE), 
               names_to = c(".value", "Age"),
               names_pattern = "(.*).{1}(..$)") 
  
#---- save datasets ----
write_csv(impute, paste0("/Users/CrystalShaw/Dropbox/Projects/",
                         "exposure_trajectories/data/",
                         "imputation_data_wide.csv"))

write_csv(impute_long, paste0("/Users/CrystalShaw/Dropbox/Projects/",
                              "exposure_trajectories/data/",
                              "imputation_data_long.csv"))




