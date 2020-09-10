#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "mice", "ghibli")

#No scientific notation
options(scipen = 999)

set.seed(20200819)

#---- Read in analytical sample ----
# imputation_data_wide <- 
#   read_csv(paste0("/Users/CrystalShaw/Dropbox/Projects/", 
#                   "exposure_trajectories/data/", 
#                   "imputation_data_wide.csv"), 
#            col_types = cols(.default = col_double(), HHIDPN = col_character(), 
#                             death = col_factor(), female = col_factor(), 
#                             hispanic = col_factor(), black = col_factor(), 
#                             other = col_factor())) 
imputation_data_long <- 
  read_csv(paste0("/Users/CrystalShaw/Dropbox/Projects/", 
                  "exposure_trajectories/data/", 
                  "imputation_data_long.csv"), 
           col_types = cols(.default = col_double(), HHIDPN = col_factor(), 
                            Wave = col_factor(), death = col_factor(), 
                            female = col_factor(), hispanic = col_factor(), 
                            black = col_factor(), other = col_factor(), 
                            smoker = col_integer())) 

#---- Indicate observed Cystatin C ----
#only interested in observed if they are in age range [60, 69]
imputation_data_long %<>% 
  mutate("observed" = ifelse(!is.na(log_CysC) & age_y_int < 70, 1, 0))

# #Sanity check-- only ages [61, 69] should have a 1
# obs_by_age <- imputation_data_long %>% group_by(age_y_int) %>% 
#   summarize_at("observed", ~sum(. == 1))

#---- Induce missingness ----
#which values [60, 69] are observed
obs_in_range <- which(imputation_data_long$observed == 1)

#create missing indicator by scenario
mcar10 <- sample(obs_in_range, size = floor(0.10*length(obs_in_range)))

imputation_data_long[, "mcar10"] <- 0
imputation_data_long[mcar10, "mcar10"] <- 1

# #Sanity check
# imputation_data_long %>% filter(age_y_int >= 70) %>% summarise_at("mcar10", sum)
# sum(imputation_data_long$mcar10)

#mask values based on missing value indicator
imputation_data_long %<>% 
  mutate("log_CysC_masked" = ifelse(mcar10 == 1, NA, log_CysC)) 

# #Sanity check
# View(imputation_data_long[, c("age_y_int", "log_CysC", "log_CysC_masked",
#                               "mcar10")])

#---- Remove people with no Cystatin C measures ----
#On this run, I've removed 78 people
no_cysc <- imputation_data_long %>% group_by(HHIDPN) %>% 
  summarise_at("log_CysC_masked", function(x) sum(!is.na(x))) %>% 
  filter(log_CysC_masked == 0)

imputation_data_long %<>% filter(!HHIDPN %in% no_cysc$HHIDPN)

#---- Missing data in predictors ----
#Predictors of Cystatin C: 
# baseline: Sex/gender, race/ethnicity, cSES, death, age at death, 
#           smoking status
# time-varying: 

#Indicate where there is missing data in the long data
impute_here_long <- is.na(imputation_data_long) %>% 
  set_colnames(colnames(imputation_data_long))*1

missingness <- t(colSums(impute_here_long)/nrow(impute_here_long)) %>% 
  as.data.frame() %>% round(., 2) %>%
  dplyr::select(-c("HHIDPN", "height_measured", "Wave", "weight_measured", 
                   "BMI_measured", "observed", "log_CysC", "mcar10", 
                   "log_CysC_masked"))

write_csv(missingness, 
          paste0("/Users/CrystalShaw/Dropbox/Projects/", 
                 "exposure_trajectories/manuscript/", 
                 "tables/missingness.csv"))

#---- check col types of dataframe ----
sapply(imputation_data_long, class)

imputation_data_long %<>% 
  mutate_at(c("observed", "mcar10"), as.factor)

#---- predictor matrix ----
predictors <- matrix(1, ncol = ncol(imputation_data_long), 
                     nrow = ncol(imputation_data_long)) %>% 
  set_colnames(colnames(imputation_data_long)) %>% 
  set_rownames(colnames(imputation_data_long))
diag(predictors) <- 0

#Don't use these as predictors
predictors[, c("HHIDPN", "log_CysC", "observed", "mcar10", "height_measured", 
               "Wave", "weight_measured", "BMI_measured", "CYSC_ADJ")] <- 0

#Don't predict these
predictors[colnames(predictors)[-which(colnames(predictors) == 
                                         "log_CysC_masked")], ] <- 0


#---- MICE ----
# #Look at missing data pattern
# md.pattern(imputation_data %>% dplyr::select(contains("CYSC")))

#FCS Imputation
#Want 25 imputations 
#maxit seems to be the number of iterations for the trace plot
num_impute = 3
imputations <- mice(imputation_data_long, m = num_impute, maxit = 5, 
                    predictorMatrix = predictors, 
                    where = impute_here_long,
                    defaultMethod = rep("norm", 4), seed = 20200812)

#check diagnostics
View(imputations$loggedEvents)
plot(imputations)
densityplot(imputations, ~ age_death_y)

#Checking
sample_original <- complete(imputations, action = 0)
sample_complete <- complete(imputations, action = 3)

colSums(is.na(sample_original))
colSums(is.na(sample_complete))

#LMM Imputation




#---- Observed vs Predicted ----
plot_data <- data.frame(matrix(nrow = nrow(imputation_data_long), 
                        ncol = num_impute)) %>% 
  set_colnames(paste0("impute", seq(1:num_impute))) %>% 
  mutate("Observed" = imputation_data_long$logCYSC_ADJ, 
         "mcar10" = imputation_data_long$mcar10)

for(i in 1:num_impute){
  plot_data[, paste0("impute", i)] <- 
    complete(imputations, action = i)[, "logCYSC_ADJ_masked"]
}

#Subset to those masked in the sample
plot_data %<>% 
  mutate("Imputed" = plot_data %>% 
           dplyr::select(contains("impute")) %>% rowMeans()) %>% 
  filter(mcar10 == 1) 
  
#plot
ggplot(data = plot_data, aes(x = Observed, y = Imputed)) + 
  geom_point(color = "#B4DAE5FF") + 
  geom_smooth(method = lm, se = FALSE, color = "#F0D77BFF") +
  geom_abline(slope = 1, intercept = 0, color = "#5C5992FF", lty = "dashed", 
              size = 1) + 
  ggtitle(paste0("Missingness Pattern: MCAR 10% \n", 
                 "Imputation Strategy: Fully Conditional Specification")) + 
  theme_minimal()  

ggsave(paste0("/Users/CrystalShaw/Dropbox/Projects/exposure_trajectories/",
              "manuscript/figures/mcar10_obs_pred.jpeg"), 
       device = "jpeg", width = 7, height = 4.5, units = "in", dpi = 300)

#---- Analytic model ----


#---- Saving output ----
#trace plots
png(paste0("/Users/CrystalShaw/Dropbox/Projects/exposure_trajectories/",
           "manuscript/figures/mcar10_traceplot.png"), 
    width = 7, height = 4.5, units = "in", res = 300)
plot(imputations)
dev.off()





