#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "mice")

#No scientific notation
options(scipen = 999)

#---- Read in analytical sample ----
imputation_data <- read_csv(paste0("/Users/CrystalShaw/Dropbox/Projects/", 
                                   "exposure_trajectories/data/", 
                                   "imputation_data.csv"))

#---- Wave indicators ----
letter_waves <- LETTERS[seq(from = 11, to = 16)] #biomarker sample + 2016 HRS

#---- Imputation dataset ----
keep <- c("HHIDPN", "raedyrs", paste0(head(letter_waves, -1), "smoken"), 
          "cses_index", "CYSC_ADJ", "A1C_ADJ", "TC_ADJ", "HDL_ADJ", "death", 
          "age_death_y", paste0(head(letter_waves, -1), "age_y"), "female", 
          "hispanic", "black", "other", "ht", "wt", "BMI", "sbp_avg", "dbp_avg")

impute <- analytical_sample %>% dplyr::select(contains(keep)) %>% 
  #leftover because of "contains"
  dplyr::select(-c("age_death_d"))

impute_long <- impute %>% 
  pivot_longer(cols = starts_with(head(letter_waves, -1), ignore.case = FALSE), 
               names_to = c("Wave", ".value"),
               names_pattern = "(.)(.*)")

#---- Missing data in predictors ----
#Predictors of Cystatin C: 
# baseline: Sex/gender, race/ethnicity, cSES
# time-varying: Age, smoking status, BMI, total cholesterol, HDL, HbA1c, sbp, 
#               dbp

cc_impute_long <- impute_long %>% filter(!is.na(CYSC_ADJ))
colSums(is.na(cc_impute_long))

#Does anyone have height but not weight or vice-versa-- Yep!
table(is.na(analytical_sample$Kht), is.na(analytical_sample$Kwt))

#---- cc_long --> cc_wide ----
cc_impute_wide <- cc_impute_long %>% 
  pivot_wider(names_from = Wave, 
              values_from = colnames(cc_impute_long)[11:ncol(cc_impute_long)],
              names_glue = "{Wave}{.value}")

#Sanity check
View(cc_impute_wide %>% dplyr::select(starts_with("K")) %>% 
       filter(is.na(KCYSC_ADJ)))

#---- Create missingness ----
#Put -999 in for those whose values of Cystatin C are missing to begin with 
#because we don't want them in the imputation model at alland should not be
#imputed

impute_long_temp <- impute_long




