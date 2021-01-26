#---- package loading + options ----
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
table_effect_ests <- 
  data.frame("Exposure" = c("CES-D Wave 4", "CES-D Wave 9", 
                            "Elevated CES-D Count", "Elevated Average CES-D"),
             "beta" = NA, "LCI" = NA, "UCI" = NA) 

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

table_effect_ests[which(table_effect_ests$Exposure == "CES-D Wave 4"), 
                  c("beta", "LCI", "UCI")] <- 
  TTEmodel_CESD4_results[nrow(TTEmodel_CESD4_results), 
                         c("estimate", "conf.low", "conf.high")]

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

table_effect_ests[which(table_effect_ests$Exposure == "CES-D Wave 9"), 
                  c("beta", "LCI", "UCI")] <- 
  TTEmodel_CESD9_results[nrow(TTEmodel_CESD9_results), 
                         c("estimate", "conf.low", "conf.high")]

#---- **Total Count Elevated CES-D ----
TTEmodel_total_CESD <- 
  coxph(Surv(survtime, observed) ~ r4age_y_int + female + hispanic + black + 
          other + ed_cat + r4mstat_cat + ever_mem + ever_arthritis + 
          ever_stroke + ever_heart + ever_lung + ever_cancer + ever_hibp + 
          ever_diabetes + r4BMI + drinking4_cat_impute + smoker + 
          total_elevated_cesd, data = CESD_data_wide)

TTEmodel_total_CESD_results <- tidy(TTEmodel_total_CESD, 
                                    exponentiate = TRUE, conf.int = TRUE)

table_effect_ests[which(table_effect_ests$Exposure == "Elevated CES-D Count"), 
                  c("beta", "LCI", "UCI")] <- 
  TTEmodel_total_CESD_results[nrow(TTEmodel_total_CESD_results), 
                              c("estimate", "conf.low", "conf.high")]

#---- **Elevated Average CES-D ----
TTEmodel_elevated_avg_CESD <- 
  coxph(Surv(survtime, observed) ~ r4age_y_int + female + hispanic + black + 
          other + ed_cat + r4mstat_cat + ever_mem + ever_arthritis + 
          ever_stroke + ever_heart + ever_lung + ever_cancer + ever_hibp + 
          ever_diabetes + r4BMI + drinking4_cat_impute + smoker + 
          avg_cesd_elevated, data = CESD_data_wide)

TTEmodel_elevated_avg_CESD_results <- tidy(TTEmodel_elevated_avg_CESD, 
                                  exponentiate = TRUE, conf.int = TRUE)

table_effect_ests[which(table_effect_ests$Exposure == "Elevated Average CES-D"), 
                  c("beta", "LCI", "UCI")] <- 
  TTEmodel_elevated_avg_CESD_results[nrow(TTEmodel_elevated_avg_CESD_results), 
                            c("estimate", "conf.low", "conf.high")]

#---- create incomplete data ----
#---- **MCAR ----
#it's easier to do this with my own code than the ampute function in MICE, which
# requires specifying all possible missing patterns you'd like it to consider
cesd_data_long <- CESD_data_wide %>% 
  dplyr::select("HHIDPN", paste0("r", seq(4, 9), "cesd")) %>% 
  pivot_longer(-c("HHIDPN"), names_to = "wave", values_to = "cesd")

mcar10_mask <- sample.int(n = nrow(cesd_data_long), 
                          size = floor(0.10*nrow(cesd_data_long)))

#mask these values in long data
cesd_mcar10 <- cesd_data_long
cesd_mcar10[mcar10_mask, "cesd"] <- NA

#create a masked dataset
mcar10 <- CESD_data_wide
mcar10[, c("HHIDPN", paste0("r", seq(4, 9), "cesd"))] <- 
  cesd_mcar10 %>% na.omit() %>% 
  pivot_wider(names_from = "wave", values_from = "cesd")

#---- check missings ----
#make sure no one is missing every cesd measure
num_missings_mcar10 <- table(rowSums(is.na(mcar10)))

#---- transformations ----
for(wave in 4:9){
  mcar10[, paste0("logr", wave, "cesd")] <- 
    log(1 + mcar10[, paste0("r", wave, "cesd")])
}

#---- imputation ----
#---- **JMVN ----
#Joint multivariate normal
#---- ***predictor matrix ----
predict <- matrix(1, nrow = 6, ncol = ncol(mcar10)) %>% 
  set_rownames(paste0("logr", seq(4, 9), "cesd")) %>% 
  set_colnames(colnames(mcar10))
#Don't use these as predictors
predict[, c("HHIDPN", "conde", "age_death_y", "r4cesd_elevated", 
            paste0("r", seq(4, 9), "cesd"), "r9cesd_elevated", 
            "total_elevated_cesd", "avg_cesd", "avg_cesd_elevated", 
            "observed", paste0("r", seq(5, 9), "age_y_int"))] <- 0

# #Use values at the current wave to predict-- not sure about this yet
# predict[, paste0("r", seq(4, 9), "BMI")] <- diag(x = 1, nrow = 6, ncol = 6)
# predict[, paste0("r", seq(4, 9), "age_y_int")] <- diag(x = 1, nrow = 6, ncol = 6)
# predict[, paste0("r", seq(4, 9), "shlt")] <- diag(x = 1, nrow = 6, ncol = 6)

#Exclude values that predict themselves
predict[, paste0("logr", seq(4, 9), "cesd")] <- 
  (diag(x = 1, nrow = 6, ncol = 6) == 0)*1

#---- ***run imputation ----
jmvn <- mice(data = mcar10, m = 5, method = "norm", predictorMatrix = predict, 
             where = is.na(mcar10), 
             blocks = as.list(paste0("logr", seq(4, 9), "cesd")))

#---- ***trace plots ----
#trace plots
png(paste0("/Users/CrystalShaw/Dropbox/Projects/exposure_trajectories/",
           "manuscript/figures/mcar10_jmvn_traceplot.png"), 
    width = 7, height = 4.5, units = "in", res = 300)
plot(imputations)
dev.off()

#---- ***visualize imputations ----

#---- ***effect estimates ----

#---- **FCS ----
#Fully conditional specification
fcs <- mice(data = mcar10, m = 5, method = "polr", )

#---- **JMVN long ----
#Longitudinal joint multivariate normal model

#---- **FCS long ----
#Longitudinal fully conditional specification


#---- save tables ----
#Round numbers in dataframe
table_effect_ests %<>% mutate(across(where(is.numeric), ~ round(., 2)))

#Save results
table_list <- list("Table 2" = table_effect_ests)
write.xlsx(table_list, file = paste0(path_to_dropbox, 
                                     "/exposure_trajectories/manuscript/", 
                                     "tables/main_text_tables.xlsx"))



#---- OLD CODE ----
imputation_vars <- c(paste0(seq(4, 9, by = 1), "BMI"), "9age_y_int", "female", 
                     "hispanic", "white", "black", "other", "height", 
                     "r9mstat_cat", "smoker", "drinking9_cat", "ed_cat", 
                     "cses_index", "death2018")

model_vars <- c("9age_y_int", "female", "hispanic", "white", "black", 
                "other", "r9mstat_cat","smoker", "drinking9_cat", "ed_cat", 
                "cses_index")

BMI_data_wide %<>% dplyr::select(all_of(c(ID, imputation_vars, model_vars)))

#---- OLD CODE ----
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
  mutate("log_CysC_masked" = ifelse(mcar10 == 1, NA, log_CysC)) %>%
  mutate("CysC_masked" = exp(log_CysC_masked))

# #Sanity check
# View(imputation_data_long[, c("age_y_int", "log_CysC", "log_CysC_masked",
#                               "mcar10")])


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


#---- Observed vs Predicted ----
plot_data <- data.frame(matrix(nrow = nrow(imputation_data_long), 
                        ncol = num_impute)) %>% 
  set_colnames(paste0("impute", seq(1:num_impute))) %>% 
  mutate("Observed" = imputation_data_long$log_CysC, 
         "mcar10" = imputation_data_long$mcar10)

for(i in 1:num_impute){
  plot_data[, paste0("impute", i)] <- 
    complete(lmm_imputations, action = i)[, "log_CysC_masked"]
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
#From the original data-- start with Cystatin C at age 65
CysC_65_model <- glm(death ~ female + hispanic + black + other + raedyrs + 
                       cses_index + smoker + BMI + CYSC_ADJ, 
                     family = binomial(link = "logit"), 
                     data = imputation_data_long %>% filter(age_y_int == 65))

tidy(CysC_65_model, exp = TRUE, conf.int = TRUE)

#Based on imputations
model_list <- vector(mode = "list", length = num_impute)

for(i in 1:num_impute){
  data = complete(imputations, action = i)
  model_list[[i]] <- 
    glm(death ~ female + hispanic + black + other + raedyrs + 
          cses_index + smoker + BMI + CysC_masked, 
        family = binomial(link = "logit"), 
        data = data %>% filter(age_y_int == 65))
}

pooled_models <- summary(pool(model_list))
pt_ests <- exp(pooled_models$estimate)
CIs <- exp(cbind(pooled_models$estimate - pooled_models$std.error, 
                 pooled_models$estimate + pooled_models$std.error))

#Based on imputations
lmm_model_list <- vector(mode = "list", length = num_impute)

for(i in 1:num_impute){
  data = complete(lmm_imputations, action = i)
  lmm_model_list[[i]] <- 
    glm(death ~ female + hispanic + black + other + raedyrs + 
          cses_index + smoker + BMI + CysC_masked, 
        family = binomial(link = "logit"), 
        data = data %>% filter(age_y_int == 65))
}

pooled_lmm_models <- summary(pool(lmm_model_list))
pt_ests <- exp(pooled_lmm_models$estimate)
CIs <- exp(cbind(pooled_lmm_models$estimate - pooled_lmm_models$std.error, 
                 pooled_lmm_models$estimate + pooled_lmm_models$std.error))


#---- Saving output ----
#trace plots
png(paste0("/Users/CrystalShaw/Dropbox/Projects/exposure_trajectories/",
           "manuscript/figures/mcar10_traceplot.png"), 
    width = 7, height = 4.5, units = "in", res = 300)
plot(imputations)
dev.off()





