#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "mice", "broom", "ghibli", 
       "ResourceSelection")

#No scientific notation
options(scipen = 999)

set.seed(20200819)

#---- Read in analytical sample ----
BMI_data_wide <- 
  read_csv(paste0("/Users/CrystalShaw/Dropbox/Projects/", 
                  "exposure_trajectories/data/", 
                  "hrs_samp_6BMI_waves4-9.csv"), 
           col_types = cols(.default = col_double(), HHIDPN = col_character(), 
                            death2018 = col_integer(), DOD = col_character(), 
                            Bday = col_character(), ed_cat = col_factor(), 
                            drop = col_logical(), r9mstat_cat = col_factor(),
                            drinking9_cat = col_factor(),
                            female = col_factor(), hispanic = col_factor(), 
                            black = col_factor(), other = col_factor(), 
                            smoker = col_integer()))

#---- Cap age at baseline at 90 ----
# #Sanity check
# max(BMI_data_wide$`4age_y_int`)

BMI_data_wide %<>% filter(`4age_y_int` <= 90)

# #Sanity check
# hist(BMI_data_wide$`4age_y_int`)
# test <- BMI_data_wide %>% dplyr::select(paste0(seq(4, 9, by = 1), "BMI")) %>% 
#   is.na() %>% sum()

#---- Sample sizes ----
num_people = nrow(BMI_data_wide)
num_obs = num_people*6

#stratified by age at baseline
by_age_baseline <- data.frame("start" = seq(50, 85, by = 5)) %>%
  mutate("end" = start + 4, 
         "n" = 0)
by_age_baseline[nrow(by_age_baseline), "end"] <- 90

for(i in 1:nrow(by_age_baseline)){
  by_age_baseline[i, "n"] <- BMI_data_wide %>% 
    filter(`4age_y_int` %in% 
             seq(by_age_baseline[i, "start"], 
                 by_age_baseline[i, "end"], by = 1)) %>% nrow()
}

#stratified by overall age
overall_ages <- BMI_data_wide %>% 
  dplyr::select(paste0(seq(4, 9, by = 1), "age_y_int")) %>% 
  pivot_longer(everything())

by_age_overall <- data.frame("start" = seq(50, 95, by = 5)) %>%
  mutate("end" = start + 4, 
         "n" = 0)
by_age_overall[nrow(by_age_overall), "end"] <- 100

for(i in 1:nrow(by_age_overall)){
  by_age_overall[i, "n"] <- overall_ages %>% 
    filter(value %in% 
             seq(by_age_overall[i, "start"], 
                 by_age_overall[i, "end"], by = 1)) %>% nrow()
}

# #Sanity check
# sum(by_age_baseline$n)
# sum(by_age_overall$n)

#---- E1 Def: Average BMI within 5-year age bands ----
#Effect of E1 on mortality within a decade of the end of follow-up

#---- **Format dataset ----
for_dataset <- c("HHIDPN", "4age_y_int", "4BMI", "age_death_y")

E1_BMI_data_long <- BMI_data_wide %>% 
  dplyr::select(all_of(for_dataset)) %>% 
  set_colnames(c("HHIDPN", "age_BMI_y", "BMI", "age_death_y"))

for(wave in 5:9){
  for_dataset <- c("HHIDPN", paste0(wave, "age_y_int"), paste0(wave, "BMI"), 
                   "age_death_y")
  
  subset <- BMI_data_wide %>% dplyr::select(all_of(for_dataset)) %>% 
    set_colnames(c("HHIDPN", "age_BMI_y", "BMI", "age_death_y"))
  
  E1_BMI_data_long %<>% rbind(subset)
}

# #Sanity check
# dim(E1_BMI_data_long)
# View(E1_BMI_data_long)

E1_BMI_data_wide <- E1_BMI_data_long %>% 
  pivot_wider(names_from = age_BMI_y, values_from = BMI, 
              names_prefix = "BMI_")

# #Sanity check
# dim(E1_BMI_data_wide)
# colnames(E1_BMI_data_wide)

#---- **Average BMI within age bands ----
for(i in seq(50, 80, by = 5)){
  E1_BMI_data_wide[, paste0("BMI_", i, "-", (i + 4))] = 
    apply(E1_BMI_data_wide %>% 
            dplyr::select(paste0("BMI_", seq(i, (i + 4), by = 1))), 
          1, function(x) mean(x, na.rm = TRUE))
}

#Get rid of columns for individual ages
E1_BMI_data_wide %<>% 
  dplyr::select(-paste0("BMI_", 
                        seq(min(E1_BMI_data_long$age_BMI_y), 
                            max(E1_BMI_data_long$age_BMI_y), by = 1)))

# #Sanity check
# View(E1_BMI_data_wide)




imputation_vars <- c(paste0(seq(4, 9, by = 1), "BMI"), "9age_y_int", "female", 
                     "hispanic", "white", "black", "other", "height", 
                     "r9mstat_cat", "smoker", "drinking9_cat", "ed_cat", 
                     "cses_index", "death2018")

model_vars <- c("9age_y_int", "female", "hispanic", "white", "black", 
                "other", "r9mstat_cat","smoker", "drinking9_cat", "ed_cat", 
                "cses_index")

BMI_data_wide %<>% dplyr::select(all_of(c(ID, imputation_vars, model_vars)))

#---- Table XX shell: Effect Estimates ----
tableXX_effect_ests <- 
  data.frame("Rownames" = c(NA, NA, "MCAR", "10%", "25%", "50%"),
             "Truth" = c("pt. est.", "(95% CI)", rep(NA, 4))) 

#---- E1 Def: BMI at wave 9 ----
E1_wide <- BMI_data_wide %>% 
  dplyr::select(-one_of(paste0(seq(4, 8, by = 1), "BMI")))

# #Distribution of BMI-- this is really symmetric
# hist(E1_wide$`9BMI`)

#---- E1 Truth ----
#True effect of BMI at wave 9 on mortality by 2018
E1_wide %<>% 
  mutate("9BMI_cat" = case_when(`9BMI` < 18.5 ~ "Underweight", 
                                `9BMI` >= 18.5 & `9BMI` < 25 ~ "Normal", 
                                `9BMI` >= 25 & `9BMI` < 30 ~ "Overweight", 
                                `9BMI` >= 30 ~ "Obese"), 
         "avg_BMI" = BMI_data_wide %>% 
           dplyr::select(paste0(seq(4, 9, by = 1), "BMI")) %>% 
           rowMeans(.), 
         "avg_BMI_cat" = case_when(avg_BMI < 18.5 ~ "Underweight", 
                                   avg_BMI >= 18.5 & avg_BMI < 25 ~ "Normal", 
                                   avg_BMI >= 25 & avg_BMI < 30 ~ "Overweight", 
                                   avg_BMI >= 30 ~ "Obese"))

E1_truth_cont <- glm(death2018 ~ `9age_y_int` + female + hispanic + black + 
                       other + cses_index + ed_cat + smoker + 
                       as.factor(drinking9_cat) + as.factor(r9mstat_cat) + 
                       `9BMI`, family = poisson(link = "log"), 
                     data = E1_wide)
tidy(E1_truth_cont, exponentiate = TRUE, conf.int = TRUE)

E1_truth_cat <- glm(death2018 ~ `9age_y_int` + female + hispanic + black + 
                      other + #cses_index + ed_cat + 
                      smoker + 
                      as.factor(drinking9_cat) + as.factor(r9mstat_cat) + 
                      `9BMI_cat`, family = poisson(link = "log"), 
                    data = E1_wide)
tidy(E1_truth_cat, exponentiate = TRUE, conf.int = TRUE)

test_cat <- glm(death2018 ~ `9age_y_int` + female + hispanic + black + 
                      other + cses_index + ed_cat + 
                      smoker + 
                      as.factor(drinking9_cat) + as.factor(r9mstat_cat) + 
                      `avg_BMI_cat`, family = poisson(link = "log"), 
                    data = E1_wide)
tidy(test_cat, exponentiate = TRUE, conf.int = TRUE)

# #Sanity check-- amount missingness in each variable
# colSums(is.na(E1_wide))

#Create missing indicator by variable
can_mask <- which(E1$age_y_int == 68)
mcar10 <- sample(can_mask, size = floor(0.10*length(can_mask)))
E1[, "mcar10"] <- 0
E1[mcar10, "mcar10"] <- 1

# #Sanity check
# #Should be 77
# sum(E1$mcar10)
# #1s should only show up at age 68
# E1 %>% group_by(age_y_int) %>% summarise_at("mcar10", ~sum(.))

E1 %<>% 
  dplyr::mutate("log_CysC_masked" = ifelse(mcar10 == 1, NA, log_CysC)) %>%
  dplyr::mutate("CysC_masked" = exp(log_CysC_masked))

#Remove people with no CysC measures-- on this run, I've removed 17 people
no_cysc <- E1 %>% dplyr::group_by(HHIDPN) %>% 
  summarise_at("log_CysC_masked", function(x) sum(!is.na(x))) %>% 
  filter(log_CysC_masked == 0)

E1 %<>% filter(!HHIDPN %in% no_cysc$HHIDPN)

# #Sanity check-- final sample size of complete data ~753
# length(unique(E1$HHIDPN))

#---- E1: Analytical Model ----
covariates <- c("female", "hispanic", "black", "other", "raedyrs", "cses_index", 
                "smoker", "BMI", "CYSC_ADJ")

subset <- E1 %>% filter(age_y_int == 68)
missingness <- subset %>% dplyr::select(covariates) %>% is.na() %>% colSums()

#What are we trying to recover?
E1_truth <- glm(death ~ female + hispanic + black + other + raedyrs + 
                  cses_index + smoker + BMI + CYSC_ADJ, 
                family = binomial(link = "logit"), 
                data = subset, 
                na.action = "na.omit")

# #Look at the model
# summary(E1_truth)

# #Tidy output
# tidy(E1_truth, exp = TRUE, conf.int = TRUE)

# #Goodness of fit-- Hosmer-Lemeshow Test (g should be > num covariates in model)
# model_subset <- na.omit(subset[, c(covariates, "death")])
# hoslem.test(model_subset$death, fitted(E1_truth), g = (length(covariates) + 2))

#---- E2 Def: Number of years elevated Cystatin C ----

#---- E2: Complete Data ----





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





