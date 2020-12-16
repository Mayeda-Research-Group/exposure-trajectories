#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "mice", "broom", "ghibli", 
       "ResourceSelection", "survival", "openxlsx")

#No scientific notation
options(scipen = 999)

set.seed(20200819)

#---- note ----
# Since the difference between win and OS, put substituted directory here
# Yingyan's directory: C:/Users/yingyan_wu
#                      C:/Users/yingyan_wu/Dropbox
# Crystal's directory: /Users/CrystalShaw
#                     ~/Dropbox/Projects

#Changing directories here will change them throughout the script
path_to_box <- "/Users/CrystalShaw"
path_to_dropbox <- "~/Dropbox/Projects"

#---- read in analytical sample ----
CESD_data_wide <- 
  read_csv(paste0(path_to_dropbox, 
                  "/exposure_trajectories/data/", 
                  "CESD_data_wide.csv"), 
           col_types = cols(.default = col_double(), HHIDPN = col_character(), 
                            death2018 = col_integer(), DOD = col_character(), 
                            Bday = col_character(), ed_cat = col_factor(), 
                            drop = col_logical(), r4mstat_cat = col_factor(), 
                            r9mstat_cat = col_factor(),
                            drinking4_cat_impute = col_factor(),
                            drinking9_cat_impute = col_factor(),
                            female = col_factor(), hispanic = col_factor(), 
                            black = col_factor(), other = col_factor(), 
                            smoker = col_integer()))

#---- select variables ----
CESD_vars <- c("HHIDPN", "4age_y", "9age_y", "female", "hispanic", "black", 
               "other", "smoker", "drinking4_cat_impute", 
               "drinking9_cat_impute", "r4mstat_cat", "r9mstat_cat", 
               "r4cesd_elevated", "r9cesd_elevated", "avg_cesd_elevated", 
               "total_elevated_cesd", "survtime", "observed")

CESD_subset <- CESD_data_wide %>% dplyr::select(all_of(CESD_vars))

#---- Check missingness ----
colSums(is.na(CESD_subset))

#---- make long dataset ----
CESD_only_long <- CESD_data_wide %>% 
  dplyr::select("HHIDPN", paste0("r", seq(4, 9), "cesd")) %>% 
  pivot_longer(contains("cesd"), values_to = "cesd", names_to = "wave")

#---- induce missingness ----
mcar10 <- sample(obs_in_range, size = floor(0.10*length(obs_in_range)))

imputation_data_long[, "mcar10"] <- 0
imputation_data_long[mcar10, "mcar10"] <- 1

mcar10 <- CESD_only
for(i in 2:ncol(mcar10)){
  mask <- sample(seq(1, nrow(mcar10)), size = ceiling(0.10*nrow(mcar10)), 
                 replace = FALSE)
  mcar10[mask, i] <- NA
}

# #Sanity check
# sum(is.na(mcar10))
# floor(0.10*nrow(mcar10)*(ncol(mcar10) - 1))
# table(rowSums(is.na(mcar10)))



# #Sanity check
# imputation_data_long %>% filter(age_y_int >= 70) %>% summarise_at("mcar10", sum)
# sum(imputation_data_long$mcar10)

#mask values based on missing value indicator
imputation_data_long %<>% 
  mutate("log_CysC_masked" = ifelse(mcar10 == 1, NA, log_CysC)) %>%
  mutate("CysC_masked" = exp(log_CysC_masked))

# #Sanity check
# View(imputation_data_long[, c("age_y_int", "log_CysC", "log_CysC_masked",
#                               "mcar10")])

#---- Remove people with no Cystatin C measures ----
#On this run, I've removed 78 people
no_cysc <- imputation_data_long %>% group_by(HHIDPN) %>% 
  summarise_at("log_CysC_masked", function(x) sum(!is.na(x))) %>% 
  filter(log_CysC_masked == 0)

imputation_data_long %<>% filter(!HHIDPN %in% no_cysc$HHIDPN)

#---- check col types of dataframe ----
sapply(imputation_data_long, class)

imputation_data_long %<>% 
  mutate_at(c("observed", "mcar10"), as.factor)

#---- Specify formulas ----
meth <- make.method(imputation_data_long)
meth["BMI"] <- "~I(weight / height^2)"
meth["CysC_masked"] <- "~I(exp(log_CysC_masked))"

#---- predictor matrix ----
pred <- make.predictorMatrix(imputation_data_long)

#Don't use these as predictors
pred[, c("HHIDPN", "log_CysC", "observed", "mcar10", "height_measured", 
         "Wave", "weight_measured", "BMI_measured", "CYSC_ADJ", 
         "CysC_masked")] <- 0

#Formulas are already specified for these
pred[c("BMI", "CysC_masked"), ] <- 0

#Do not need imputations for these-- do I even need to specify this?
pred[c("HHIDPN", "female", "hispanic", "black", "other", "raedyrs", "death", 
       "height", "height_measured", "Wave", "age_y_int", "weight_measured", 
       "BMI_measured", "CYSC_ADJ", "log_CysC", "observed", "mcar10"), ] <- 0

#---- Missing data in predictors ----
#Predictors of Cystatin C: 
# baseline: Sex/gender, race/ethnicity, cSES, death, age at death, 
#           smoking status
# time-varying: 

#Indicate where there is missing data in the long data
impute_here_long <- is.na(imputation_data_long) %>% 
  set_colnames(colnames(imputation_data_long))*1

missingness <- t(colSums(impute_here_long)/nrow(impute_here_long)) %>% 
  as.data.frame() %>% round(., 2) 

#Indicate where we don't want imputations-- original data
impute_here_long[, c("CYSC_ADJ", "log_CysC")] <- 0

missingness_table <- missingness %>%
  dplyr::select(-c("HHIDPN", "height_measured", "Wave", "weight_measured", 
                   "BMI_measured", "observed", "log_CysC", "mcar10", 
                   "log_CysC_masked", "CysC_masked"))

write_csv(missingness_table, 
          paste0("/Users/CrystalShaw/Dropbox/Projects/", 
                 "exposure_trajectories/manuscript/", 
                 "tables/missingness.csv"))

#---- MICE ----
# #Look at missing data pattern
# md.pattern(imputation_data %>% dplyr::select(contains("CYSC")))

#FCS Imputation
#Want 25 imputations 
#maxit seems to be the number of iterations for the trace plot
num_impute = 3
imputations <- mice(imputation_data_long, m = num_impute, maxit = 5, 
                    predictorMatrix = pred, 
                    where = impute_here_long,
                    defaultMethod = rep("norm", 4), seed = 20200812)

# #check diagnostics
# View(imputations$loggedEvents)
# plot(imputations)
# densityplot(imputations, ~ age_death_y)

# #Checking
# sample_original <- complete(imputations, action = 0)
# sample_complete <- complete(imputations, action = 3)
# 
# colSums(is.na(sample_original))
# colSums(is.na(sample_complete))


#---- MI ----


#---- imputed exposure defs ----


