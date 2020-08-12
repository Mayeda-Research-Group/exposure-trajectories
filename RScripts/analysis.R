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
                                   "imputation_data.csv"), 
                            col_types = cols(.default = col_double(), 
                                             HHIDPN = col_character()))

#---- Missing data in predictors ----
#Predictors of Cystatin C: 
# baseline: Sex/gender, race/ethnicity, cSES, death, age at death
# time-varying: 

# #Looking at what else needs to be imputed-- no other missing data at this time
# colSums(is.na(imputation_data))






