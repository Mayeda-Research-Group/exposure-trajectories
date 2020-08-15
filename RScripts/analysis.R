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
                            other = col_factor())) %>% 
  dplyr::select(-age_death_y)

#---- Missing data in predictors ----
#Predictors of Cystatin C: 
# baseline: Sex/gender, race/ethnicity, cSES, death, age at death
# time-varying: 

#Define where we want data imputed
impute_here <- is.na(imputation_data) %>% 
  set_colnames(colnames(imputation_data))*1

#Check where there is missing data
colSums(impute_here)/3204

impute_here[, c("age_death_y", paste0("logCYSC_ADJ_", c(60, seq(70, 78))))] <- 0

#---- predictor matrix ----
predictors <- matrix(1, ncol = ncol(imputation_data), 
                     nrow = ncol(imputation_data)) %>% 
  set_colnames(colnames(imputation_data)) %>% 
  set_rownames(colnames(imputation_data))
diag(predictors) <- 0

#Don't use these as predictors
predictors[, "HHIDPN"] <- 0
predictors[, "logCYSC_ADJ_60"] <- 0

#Don't predict these
predictors[colnames(imputation_data)[-which(colnames(imputation_data) %in% 
                                   paste0("logCYSC_ADJ_", seq(60, 69)))], ] <- 0

#---- MICE ----
# #Look at missing data pattern
# md.pattern(imputation_data %>% dplyr::select(contains("CYSC")))

md.pattern(complete(imputations, action = 2)[, c("age_death_y", "logCYSC_ADJ_62")])

#Imputation
#Want 25 imputations 
#maxit seems to be the number of iterations for the trace plot
imputations <- mice(imputation_data, m = 2, maxit = 25, 
                    predictorMatrix = predictors, 
                    #where = impute_here,
                    defaultMethod = rep("norm.predict", 4), seed = 20200812)

#check diagnostics
View(imputations$loggedEvents)
plot(imputations)
densityplot(imputations, ~ age_death_y)

#Checking
sample_complete <- complete(imputations, action = 1)
max(table(sample_complete$age_death_y))




