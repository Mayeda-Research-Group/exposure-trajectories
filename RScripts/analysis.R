#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse")

#No scientific notation
options(scipen = 999)

#---- Read in analytical sample ----
analytical_sample <- read_csv(paste0("/Users/CrystalShaw/Dropbox/Projects/", 
                                     "exposure_trajectories/data/", 
                                     "hrs_samp_alive_70_cysc_60_70.csv"))

#---- Wave indicators ----
letter_waves <- LETTERS[seq(from = 11, to = 16)] #biomarker sample + 2016 HRS

#---- Imputation dataset ----
keep <- c("HHIDPN", "raedyrs", paste0("r", head(letter_waves, -1), "smoken"), 
          "cses_index", "CYSC_ADJ", "A1C_ADJ", "TC_ADJ", "HDl_ADJ", "death", 
          "age_death_y", paste0(head(letter_waves, -1), "age_y"), "female", 
          "hispanic", "black", "other", "BMI", "sbp_avg", "dbp_avg")

impute <- analytical_sample %>% dplyr::select(contains(keep))

impute_long <- 

#---- Observed Cystatin C ----
analytical_sample_long <- analytical_sample %>%
  pivot_longer(cols = contains("CYSC_ADJ"), 
               names_to = "CysC_wave", values_to = "CysC")

cc_long <- analytical_sample_long %>% filter(!is.na(CysC))

cc_wide <- cc_long %>% pivot_wider(names_from = "CysC_wave", 
                                   values_from = "CysC")

#---- Missing data in predictors ----
#Predictors of Cystatin C: 
# baseline: Sex/gender, race/ethnicity, cSES
# time-varying: Age, smoking status, BMI, total cholesterol, HDL, HbA1c, sbp, 
#               dbp


