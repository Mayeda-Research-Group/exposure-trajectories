#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "mice")

#No scientific notation
options(scipen = 999)

#---- Read in analytical sample ----
imputation_data <- 
  read_csv(paste0("/Users/CrystalShaw/Dropbox/Projects/", 
                  "exposure_trajectories/data/", 
                  "imputation_data.csv"), 
           col_types = cols(.default = col_double(), HHIDPN = col_character(), 
                            death = col_factor(), female = col_factor(), 
                            hispanic = col_factor(), black = col_factor(), 
                            other = col_factor()))

#---- Missing data in predictors ----
#Predictors of Cystatin C: 
# baseline: Sex/gender, race/ethnicity, cSES, death, age at death
# time-varying: 

#Define where we want data imputed
impute_here <- is.na(imputation_data) %>% 
  set_colnames(colnames(imputation_data))*1

#Check where there is missing data
colSums(impute_here)

impute_here[, c("age_death_y", paste0("logCYSC_ADJ_", seq(70, 78)))] <- 0

#---- predictor matrix ----
predictors <- matrix(1, ncol = ncol(imputation_data), 
                     nrow = ncol(imputation_data)) %>% 
  set_colnames(colnames(imputation_data)) %>% 
  set_rownames(colnames(imputation_data))
diag(predictors) <- 0

predictors[, "HHIDPN"] <- 0

#---- MICE ----
# #Look at missing data pattern
# md.pattern(imputation_data %>% dplyr::select(contains("CYSC")))

#Imputation
#Want 25 imputations 
#maxit seems to be the number of iterations for the trace plot
imputations <- mice(imputation_data, m = 2, maxit = 25, where = impute_here,
                    defaultMethod = rep("pmm", 4), seed = 20200812)

#check diagnostics
imputations$loggedEvents
plot(imputations)
#densityplot(imp_pmm_all_brfss, ~income)




