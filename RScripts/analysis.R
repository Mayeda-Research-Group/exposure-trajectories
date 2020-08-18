#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "mice")

#No scientific notation
options(scipen = 999)

#---- Read in analytical sample ----
imputation_data_wide <- 
  read_csv(paste0("/Users/CrystalShaw/Dropbox/Projects/", 
                  "exposure_trajectories/data/", 
                  "imputation_data_wide.csv"), 
           col_types = cols(.default = col_double(), HHIDPN = col_character(), 
                            death = col_factor(), female = col_factor(), 
                            hispanic = col_factor(), black = col_factor(), 
                            other = col_factor())) 
imputation_data_long <- 
  read_csv(paste0("/Users/CrystalShaw/Dropbox/Projects/", 
                  "exposure_trajectories/data/", 
                  "imputation_data_long.csv"), 
           col_types = cols(.default = col_double(), HHIDPN = col_character(), 
                            death = col_factor(), female = col_factor(), 
                            hispanic = col_factor(), black = col_factor(), 
                            other = col_factor())) 

#---- Missing data in predictors ----
#Predictors of Cystatin C: 
# baseline: Sex/gender, race/ethnicity, cSES, death, age at death
# time-varying: 

#Indicate where there is missing data in the wide data
impute_here_wide <- is.na(imputation_data_wide) %>% 
  set_colnames(colnames(imputation_data_wide))*1
colSums(impute_here_wide)/nrow(impute_here_wide)

#Get rid of lines with missing data in ages [70, max_age]
max_age <- max(imputation_data_long$Age, na.rm = TRUE)
imputation_data_long %<>% 
  mutate("remove" = ifelse(Age %in% seq(70, max_age) & 
                             is.na(logCYSC_ADJ), 1, 0)) %>% 
  filter(remove == 0) %>% 
  dplyr::select(-"remove")

#Indicate where there is missing data in the long data
impute_here_long <- is.na(imputation_data_long) %>% 
  set_colnames(colnames(imputation_data_long))*1

colSums(impute_here_long)/nrow(impute_here_long)


# #Define where we want imputations in wide data
# impute_here_wide[, c("age_death_y", 
#                      paste0("logCYSC_ADJ_", c(60, seq(70, 78))))] <- 0


#---- predictor matrix ----
predictors <- matrix(1, ncol = ncol(imputation_data_long), 
                     nrow = ncol(imputation_data_long)) %>% 
  set_colnames(colnames(imputation_data_long)) %>% 
  set_rownames(colnames(imputation_data_long))
diag(predictors) <- 0

#Don't use these as predictors
predictors[, "HHIDPN"] <- 0

#Don't predict these
predictors[colnames(predictors)[-which(colnames(predictors) == 
                                         "logCYSC_ADJ")], ] <- 0

# #Don't predict these in wide dataset
# predictors[colnames(imputation_data)[-which(colnames(imputation_data) %in% 
#                                    paste0("logCYSC_ADJ_", seq(60, 69)))], ] <- 0

#---- MICE ----
# #Look at missing data pattern
# md.pattern(imputation_data %>% dplyr::select(contains("CYSC")))

md.pattern(complete(imputations, 
                    action = 2)[, c("age_death_y", "logCYSC_ADJ_62")])

#Imputation
#Want 25 imputations 
#maxit seems to be the number of iterations for the trace plot
imputations <- mice(imputation_data_long, m = 5, maxit = 100, 
                    predictorMatrix = predictors, 
                    where = impute_here_long,
                    defaultMethod = rep("norm", 4), seed = 20200812)

#check diagnostics
View(imputations$loggedEvents)
plot(imputations)
densityplot(imputations, ~ age_death_y)

#Checking
sample_complete <- complete(imputations, action = 1)
max(table(sample_complete$age_death_y))

#---- Analytic model ----
for(i in 0:2){
  View(complete(imputations, action = i))
}




