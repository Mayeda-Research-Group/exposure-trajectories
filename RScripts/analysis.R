#---- package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "mice", "broom", "ResourceSelection", 
       "survival", "openxlsx")

#No scientific notation
options(scipen = 999)

set.seed(20200819)

#---- source scripts ----
source(here("RScripts", "mask_impute_pool.R"))

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
                            death2018 = col_integer(), 
                            ed_cat = col_factor(), 
                            r4mstat_cat = col_factor(), 
                            r9mstat_cat = col_factor(),
                            drinking4_cat_impute = col_factor(),
                            drinking9_cat_impute = col_factor(),
                            female = col_factor(), hispanic = col_factor(), 
                            black = col_factor(), other = col_factor(), 
                            smoker = col_integer()))

#---- Table 2 shell: Effect Estimates ----
exposures <- c("CES-D Wave 4", "CES-D Wave 9", "Elevated CES-D Count", 
               "Elevated Average CES-D")
methods <- c("JMVN", "FCS", "JMVN Long", "FCS Long")
mechanisms <- c("MCAR", "MAR", "NMAR")
mask_props <- c(.10, .25, .50)

table_effect_ests <- 
  data.frame("Exposure" = rep(rep(exposures, 13), 3),
             "beta" = NA, "LCI" = NA, "UCI" = NA, 
             "Method" = rep(c(rep("Truth", 4), rep(methods, each = 12)), 3), 
             "Missingness" = rep(c(rep("0%", 4), 
                                   rep(rep(paste0(mask_props*100, "%"), each = 4), 
                                       4)), 3), 
             "Type" = c(c(rep("Truth", 4), rep("MCAR", 48)), 
                        c(rep("Truth", 4), rep("MAR", 48)), 
                        c(rep("Truth", 4), rep("NMAR", 48))))

#---- truth ----
#---- **CES-D Wave 4 ----
TTEmodel_CESD4 <- 
  coxph(Surv(survtime, observed) ~ r4age_y_int + female + hispanic + black + 
          other + ed_cat + r4mstat_cat + ever_mem + ever_arthritis + 
          ever_stroke + ever_heart + ever_lung + ever_cancer + ever_hibp + 
          ever_diabetes + r4BMI + drinking4_cat_impute + smoker + 
          r4cesd_elevated, data = CESD_data_wide)

#summary(TTEmodel_CESD4)

TTEmodel_CESD4_results <- tidy(TTEmodel_CESD4, 
                               exponentiate = TRUE, conf.int = TRUE)

table_effect_ests[which(table_effect_ests$Exposure == "CES-D Wave 4" & 
                          table_effect_ests$Method == "Truth"), 
                  c("beta", "LCI", "UCI")] <- 
  c(TTEmodel_CESD4_results[nrow(TTEmodel_CESD4_results), 
                           c("estimate", "conf.low", "conf.high")])

#---- **CES-D Wave 9 ----
TTEmodel_CESD9 <- 
  coxph(Surv(survtime, observed) ~ r9age_y_int + female + hispanic + black + 
          other + ed_cat + r9mstat_cat + ever_mem + ever_arthritis + 
          ever_stroke + ever_heart + ever_lung + ever_cancer + ever_hibp + 
          ever_diabetes + r9BMI + drinking9_cat_impute + smoker + 
          r9cesd_elevated, data = CESD_data_wide)

#summary(TTEmodel_CESD9)

TTEmodel_CESD9_results <- tidy(TTEmodel_CESD9, 
                               exponentiate = TRUE, conf.int = TRUE)

table_effect_ests[which(table_effect_ests$Exposure == "CES-D Wave 9" & 
                          table_effect_ests$Method == "Truth"), 
                  c("beta", "LCI", "UCI")] <- 
  c(TTEmodel_CESD9_results[nrow(TTEmodel_CESD9_results), 
                           c("estimate", "conf.low", "conf.high")])

#---- **Total Count Elevated CES-D ----
TTEmodel_total_CESD <- 
  coxph(Surv(survtime, observed) ~ r4age_y_int + female + hispanic + black + 
          other + ed_cat + r4mstat_cat + ever_mem + ever_arthritis + 
          ever_stroke + ever_heart + ever_lung + ever_cancer + ever_hibp + 
          ever_diabetes + r4BMI + drinking4_cat_impute + smoker + 
          total_elevated_cesd, data = CESD_data_wide)

TTEmodel_total_CESD_results <- tidy(TTEmodel_total_CESD, 
                                    exponentiate = TRUE, conf.int = TRUE)

table_effect_ests[which(table_effect_ests$Exposure == "Elevated CES-D Count" & 
                          table_effect_ests$Method == "Truth"), 
                  c("beta", "LCI", "UCI")] <- 
  c(TTEmodel_total_CESD_results[nrow(TTEmodel_total_CESD_results), 
                                c("estimate", "conf.low", "conf.high")])

#---- **Elevated Average CES-D ----
TTEmodel_elevated_avg_CESD <- 
  coxph(Surv(survtime, observed) ~ r4age_y_int + female + hispanic + black + 
          other + ed_cat + r4mstat_cat + ever_mem + ever_arthritis + 
          ever_stroke + ever_heart + ever_lung + ever_cancer + ever_hibp + 
          ever_diabetes + r4BMI + drinking4_cat_impute + smoker + 
          avg_cesd_elevated, data = CESD_data_wide)

TTEmodel_elevated_avg_CESD_results <- tidy(TTEmodel_elevated_avg_CESD, 
                                           exponentiate = TRUE, conf.int = TRUE)

table_effect_ests[which(table_effect_ests$Exposure == "Elevated Average CES-D" & 
                          table_effect_ests$Method == "Truth"), 
                  c("beta", "LCI", "UCI")] <- 
  c(TTEmodel_elevated_avg_CESD_results[nrow(TTEmodel_elevated_avg_CESD_results), 
                                       c("estimate", "conf.low", "conf.high")])

#---- create one set of imputations for plot ----
all_combos <- expand_grid(mechanisms, methods, mask_props) %>% 
  mutate("mask_percent" = paste0(100*mask_props, "%"))

apply(all_combos[1:3, ], 1, function(x) 
  mask_impute_pool(CESD_data_wide, mechanism = x[1], method = x[2], 
                   mask_percent = x[4], 
                   num_impute = 5, save = "yes"))

#---- get pooled effect estimates ----
test <- mask_impute_pool(CESD_data_wide, mechanism = "MCAR", method = "JMVN", 
                         mask_percent = "10%", num_impute = 5, save = "no")

#---- save tables ----
#Round numbers in dataframe
table_effect_ests %<>% mutate(across(where(is.numeric), ~ round(., 2)))

#Save results
table_list <- list("Table 2" = table_effect_ests)
write.xlsx(table_list, file = paste0(path_to_dropbox, 
                                     "/exposure_trajectories/manuscript/", 
                                     "tables/main_text_tables.xlsx"))

#---- MICE long ----
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

#LMM Imputation
pred_lmm <- make.predictorMatrix(imputation_data_long)

#Don't use these as predictors
pred_lmm[, c("log_CysC", "observed", "mcar10", "height_measured", 
             "Wave", "weight_measured", "BMI_measured", "CYSC_ADJ", 
             "CysC_masked")] <- 0

#Formulas are already specified for these
pred_lmm[c("BMI", "CysC_masked"), ] <- 0

#Specify class variable
pred_lmm[, "HHIDPN"] <- -2

#Specify fixed effects
pred_lmm[, c("female", "hispanic", "black", "other", "raedyrs", "cses_index", 
             "smoker", "height")] <- 1

#Specify random effect
pred_lmm[, c("weight", "age_y_int", "A1C_ADJ", "TC_ADJ", "HDL_ADJ", 
             "bpsys", "bpdia", "BMI")] <- 2


#Do not need imputations for these-- do I even need to specify this?
pred_lmm[c("female", "hispanic", "black", "other", "raedyrs", "death", 
           "height", "height_measured", "Wave", "age_y_int", "weight_measured", 
           "BMI_measured", "CYSC_ADJ", "log_CysC", "observed", "mcar10"), ] <- 0

lmm_imputations <- mice(imputation_data_long, m = num_impute, maxit = 5, 
                        predictorMatrix = pred_lmm, 
                        where = impute_here_long,
                        defaultMethod = rep("2l.pan", 4), seed = 20200812)

# #check diagnostics
# View(lmm_imputations$loggedEvents)
# plot(lmm_imputations)

# #Checking
# sample_original <- complete(imputations, action = 0)
# sample_complete <- complete(imputations, action = 3)
#
# colSums(is.na(sample_original))
# colSums(is.na(sample_complete))








