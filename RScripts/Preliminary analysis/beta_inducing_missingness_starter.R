# Betas for missingness models, providing basis for the fixed effect size used
# in 2a_MAR_MNAR_beta_0_optimizer.R
# Data input: trk2018tr_r.sas7bdat, randhrs1992_2018v1.dta, cses_measures.dta
# Data output: CESD_missing_model_betas_in_sampled.csv, predict_death2018_betas.csv
# Author: CS & YW

#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}
p_load("here", "readr", "tidyverse", "magrittr", "plyr", "haven", "labelled",
       "broom", "kableExtra")

#No scientific notation
options(scipen = 999)

#---- Note ----
# Since the difference between win and OS, put substituted directory here
# Yingyan's directory: C:/Users/Yingyan Wu/Box
#                      C:/Users/Yingyan Wu/Dropbox
# Crystal's directory: /Users/crystalshaw/Library/CloudStorage/Box-Box
#                     ~/Dropbox/Projects

#Changing directories here will change them throughout the script
path_to_box <- "C:/Users/Yingyan Wu/Box"
path_to_dropbox <- "C:/Users/Yingyan Wu/Dropbox"

#---- source scripts ----
source(here::here("RScripts", "functions", "impute_ages.R"))
source(here::here("Rscripts", "functions", "impute_condition.R"))

#---- wave mapping between HRS and RAND ----
#Wave Year | HRS Core Data | RAND
# 1992 | V | 1
# 1994 | W | 2
# 1996 | E | 3
# 1998 | F | 4
# 2000 | G | 5
# 2002 | H | 6
# 2004 | J | 7
# 2006 | K | 8
# 2008 | L | 9
# 2010 | M | 10
# 2012 | N | 11
# 2014 | O | 12
# 2016 | P | 13

number_waves <- seq(1, 13, by = 1) 

#---- read in HRS tracker ----
hrs_tracker <-
  read_sas(paste0(path_to_box, "/HRS/tracker/trk2018v2a/",
                  "trk2018tr_r.sas7bdat")) %>%
  select("HHID", "PN", "PIWTYPE", "PALIVE", "QIWTYPE", "QALIVE", 
         paste0(c("F", "G", "H", "J", "K", "L"), "ALIVE")) %>%
  unite("HHIDPN", c("HHID", "PN"), sep = "", remove = TRUE) %>%
  mutate_at("HHIDPN", as.numeric)

#---- read in RAND file ----
#reading in STATA file because SAS file wouldn't load
#Variables of interest: 
## Demographics: HHIDPN, gender, race, hispanic, birth month, 
##               birth year, birth date, death month, death year, death date,
##               age data in months (biomarker waves), years of education, 
##               highest degree (masked),  
## Health: weight (kg; measured; waves 8+),
##         weight (kg; self-report),
##         height (m; measured; waves 8+),
##         height (m; self-report),
##         BMI (measured; waves 8+),
##         BMI (self-report), 
##         waist circumference,
##         BP (systolic; waves 8+), 
##         BP (diastolic; waves 8+), 
##         reports high blood pressure this wave,
##         reports diabetes this wave,
##         diabetes ever/never,
##         reports cancer this wave,
##         reports stroke this wave,
##         reports heart problems this wave，
##         report memory prob this wv(4-9), 
##         report lung disease this wv,
##         count of chronic conditions
##         CESD (depression; from wave 2-13)
##         Self-reported health
## Health Behaviors: current smoker 
##                   number of days drinking per week (waves 3+)
##                   number of drinks per day (waves 3+)
## 
# Note: Dates are formatted as SAS dates (days from January 1, 1960)

rand_variables <- c("hhidpn", "ragender", "raracem", "rahispan", "rabmonth", 
                    "rabyear", "rabdate", "radmonth", "radyear", "raddate",
                    paste0("r", number_waves, "agem_e"), "raedyrs", 
                    "raedegrm", 
                    paste0("r", number_waves, "mstat"),
                    paste0("r", seq(8, 13, by = 1), "pmwght"), 
                    paste0("r", number_waves, "weight"),
                    paste0("r", seq(8, 13, by = 1), "pmhght"), 
                    paste0("r", number_waves, "height"),
                    paste0("r", number_waves, "bmi"), 
                    paste0("r", seq(8, 13, by = 1), "pmbmi"),
                    paste0("r", seq(8, 13, by = 1), "pmwaist"),
                    paste0("r", seq(8, 13, by = 1), "bpsys"), 
                    paste0("r", seq(8, 13, by = 1), "bpdia"),
                    paste0("r", number_waves, "hibp"),
                    paste0("r", number_waves, "hibpe"),
                    paste0("r", number_waves, "diab"),
                    paste0("r", number_waves, "diabe"),
                    paste0("r", number_waves, "cancr"),
                    paste0("r", number_waves, "cancre"),
                    paste0("r", number_waves, "strok"), 
                    paste0("r", number_waves, "stroke"),
                    paste0("r", number_waves, "heart"),
                    paste0("r", number_waves, "hearte"),
                    paste0("r", number_waves, "lung"),
                    paste0("r", number_waves, "lunge"),
                    paste0("r", seq(4 ,9 , by = 1), "memry"),
                    paste0("r", seq(4 ,9 , by = 1), "memrye"),
                    paste0("r", number_waves, "conde"),
                    paste0("r", number_waves, "smoken"), 
                    paste0("r", seq(3, 13, by = 1), "drinkd"),
                    paste0("r", seq(3, 13, by = 1), "drinkn"),
                    paste0("r", seq(2, 13, by = 1), "cesd"),
                    paste0("r", number_waves, "shlt"))

RAND <- read_dta(paste0(path_to_box, "/HRS/RAND_longitudinal/STATA/", 
                        "randhrs1992_2018v1.dta"), 
                 col_select = all_of(rand_variables)) 

colnames(RAND)[1] <- "HHIDPN" #For merging

#Remove labeled data format
val_labels(RAND) <- NULL

#---- read in Anusha Vable's CSES index ----
cSES <- read_dta(paste0(path_to_dropbox, "/exposure_trajectories/data/", 
                        "cSES measures/cses_measures.dta"), 
                 col_select = c("hhid", "pn", "cses_index")) %>% 
  unite(col = "HHIDPN", c("hhid", "pn"), sep = "") %>% 
  mutate_at("HHIDPN", as.numeric)

#---- merge datasets ----
#Use this to subset RAND data
hrs_samp <- 
  join_all(c(list(RAND, cSES, hrs_tracker)), by = "HHIDPN", type = "left") 

#---- Data cleaning ----
hrs_samp %<>% mutate("death2018" = ifelse(is.na(raddate), 0, 1))

#---- **age ----
age_m <- hrs_samp %>% dplyr::select(contains("agem_e")) %>% 
  apply(., 1, impute_ages) %>% t() 

#Exact ages
hrs_samp[, paste0("r", number_waves, "age_y")] <- age_m/12
#Ages rounded down to nearest year
hrs_samp[, paste0("r", number_waves, "age_y_int")] <- floor(age_m/12)

# #Sanity check
# View(hrs_samp[, c(paste0(number_waves, "age_y"), 
#                   paste0(number_waves, "age_y_int"))])

#Drop RAND age variables
hrs_samp %<>% dplyr::select(-paste0("r", number_waves, "agem_e"))

#---- **CESD ----
# Create indicators for whether CESD is missing at each wave
hrs_samp[, paste0("r", seq(2, 13), "cesd_missing")] <- NA

cesd_mat <- hrs_samp %>% select(contains("cesd"))
for (j in 1:length(seq(2, 13))) {
  cesd_mat[, j + length(seq(2, 13))] <- ifelse(is.na(cesd_mat[, j]), 1, 0)
}

hrs_samp[, colnames(cesd_mat)] <- cesd_mat

# CESD missingness proportion
cesd_missing <- cesd_mat %>%
  select(contains("missing")) %>%
  colSums(is.na(.))
cesd_missing/nrow(hrs_samp) 

# Imputing chronic conditions
#---- ** diabetes ----
hrs_samp <- impute_chronic_condition("diabe", paste0("r", seq(1, 9), "diabe"),
                                     seq(1, 9), hrs_samp)

#---- **high bp ----
hrs_samp <- impute_chronic_condition("hibpe", paste0("r", seq(1, 9), "hibpe"),
                                     seq(1, 9), hrs_samp)

#---- **cancer ----
hrs_samp <- impute_chronic_condition("cancre", paste0("r", seq(1, 9), "cancre"),
                                     seq(1, 9), hrs_samp)

#---- **lung ----
hrs_samp <- impute_chronic_condition("lunge", paste0("r", seq(1, 9), "lunge"),
                                     seq(1, 9), hrs_samp)

#---- **heart ----
hrs_samp <- impute_chronic_condition("hearte", paste0("r", seq(1, 9), "hearte"),
                                     seq(1, 9), hrs_samp)

#---- **stroke ----
hrs_samp <- impute_chronic_condition("stroke", paste0("r", seq(1, 9), "stroke"),
                                     seq(1, 9), hrs_samp)

#---- **memory ----
#For memory problems, data starts from wave 4.
hrs_samp <- impute_chronic_condition("memrye", paste0("r", seq(4, 9), "memrye"),
                                     seq(4, 9), hrs_samp)
#---- **sum of conditions ----
# wave-specific r(wave)conde
cond_mat <- hrs_samp %>%
  dplyr::select(contains("_impute"))

waves <- seq(1, 9)
for(j in 1:length(waves)){
  wave <- waves[j] 
  cond_mat[, paste0("r", wave , "conde", "_impute")] <- 
    rowSums(cond_mat %>% dplyr::select(contains(paste0("r", wave ))), 
            na.rm = TRUE)
}

hrs_samp[, colnames(cond_mat %>% select(contains("conde_impute")))] <- 
  cond_mat %>% select(contains("conde_impute"))

#---- Preliminary models ----
# DO NOT DELETE DURING CODE CLEAN-UP
# logit(P(missing CESD this wave)) =
# \beta_0 + \beta_1*age at current wave + \beta_2*value of previous CESD +
# \beta_3* chronic condition count (at last wave)

#num missing cesd among those sampled
wave4_sampled <- hrs_samp %>% filter(!is.na(FALIVE))
mean(wave4_sampled$r4cesd_missing)
summary(r4cesdmissing_mod <- glm(r4cesd_missing ~
                                   r4age_y_int + r3cesd + r3conde_impute +
                                   death2018 - 1,
                                 family = binomial(link = "logit"),
                                 data = wave4_sampled))
r4results <- tidy(r4cesdmissing_mod, exponentiate = TRUE, conf.int = TRUE)

wave5_sampled <- hrs_samp %>% filter(!is.na(GALIVE))
mean(wave5_sampled$r5cesd_missing)
summary(r5cesdmissing_mod <- glm(r5cesd_missing ~
                                   r5age_y_int + r4cesd + r4conde_impute +
                                   death2018 - 1,
                                 family = binomial(link = "logit"),
                                 data = wave5_sampled))
r5results <- tidy(r5cesdmissing_mod, exponentiate = TRUE, conf.int = TRUE)

wave6_sampled <- hrs_samp %>% filter(!is.na(HALIVE))
mean(wave6_sampled$r6cesd_missing)
summary(r6cesdmissing_mod <- glm(r6cesd_missing ~
                                   r6age_y_int + r5cesd + r5conde_impute +
                                   death2018 - 1,
                                 family = binomial(link = "logit"),
                                 data = wave6_sampled))
r6results <- tidy(r6cesdmissing_mod, exponentiate = TRUE, conf.int = TRUE)

wave7_sampled <- hrs_samp %>% filter(!is.na(JALIVE))
mean(wave7_sampled$r7cesd_missing)
summary(r7cesdmissing_mod <- glm(r7cesd_missing ~
                                   r7age_y_int + r6cesd + r6conde_impute +
                                   death2018 - 1,
                                 family = binomial(link = "logit"),
                                 data = wave7_sampled))
r7results <- tidy(r7cesdmissing_mod, exponentiate = TRUE, conf.int = TRUE)

wave8_sampled <- hrs_samp %>% filter(!is.na(KALIVE))
mean(wave8_sampled$r8cesd_missing)
summary(r8cesdmissing_mod <- glm(r8cesd_missing ~
                                   r8age_y_int + r7cesd + r7conde_impute +
                                   death2018 - 1,
                                 family = binomial(link = "logit"),
                                 data = wave8_sampled))
r8results <- tidy(r8cesdmissing_mod, exponentiate = TRUE, conf.int = TRUE)

wave9_sampled <- hrs_samp %>% filter(!is.na(LALIVE))
mean(wave9_sampled$r9cesd_missing)
summary(r9cesdmissing_mod <- glm(r9cesd_missing ~
                                   r9age_y_int + r8cesd + r8conde_impute +
                                   death2018 - 1,
                                 family = binomial(link = "logit"),
                                 data = wave9_sampled))
r9results <- tidy(r9cesdmissing_mod, exponentiate = TRUE, conf.int = TRUE)

results_tbl <- tibble(
  variables = c("age at current wave", "previous CESD value",
                "previous chronic condition count", "death2018"),
  r4beta = round(r4results$estimate, 4),
  r5beta = round(r5results$estimate, 4),
  r6beta = round(r6results$estimate, 4),
  r7beta = round(r7results$estimate, 4),
  r8beta = round(r8results$estimate, 4),
  r9beta = round(r9results$estimate, 4)
)

results_tbl %>%
  kbl(caption =
        "Exponentiated betas of the CESD missing model (wave 4 - 9)") %>%
  kable_classic(full_width = F, html_font = "Arial")

write_csv(results_tbl, paste0(path_to_dropbox,
                              "/exposure_trajectories/data/",
                              "CESD_missing_model_betas_in_sampled.csv"))

#---- **predict outcome?? ----
# DO NOT DELETE DURING CODE CLEAN-UP
summary(r4outcome_mod <- glm(death2018 ~
                               r4age_y_int + r3cesd + r3conde_impute,
                             family = binomial(link = "logit"),
                             data = hrs_samp))
r4outcome_results <- tidy(r4outcome_mod, exponentiate = TRUE, conf.int = TRUE)

summary(r5outcome_mod <- glm(death2018 ~
                               r5age_y_int + r4cesd + r4conde_impute,
                             family = binomial(link = "logit"),
                             data = hrs_samp))
r5outcome_results <- tidy(r5outcome_mod, exponentiate = TRUE, conf.int = TRUE)

summary(r6outcome_mod <- glm(death2018 ~
                               r6age_y_int + r5cesd + r5conde_impute,
                             family = binomial(link = "logit"),
                             data = hrs_samp))
r6outcome_results <- tidy(r6outcome_mod, exponentiate = TRUE, conf.int = TRUE)

summary(r7outcome_mod <- glm(death2018 ~
                                   r7age_y_int + r6cesd + r6conde_impute,
                                 family = binomial(link = "logit"),
                                 data = hrs_samp))
r7outcome_results <- tidy(r7outcome_mod, exponentiate = TRUE, conf.int = TRUE)

summary(r8outcome_mod <- glm(death2018 ~
                                   r8age_y_int + r7cesd + r7conde_impute,
                                 family = binomial(link = "logit"),
                                 data = hrs_samp))
r8outcome_results <- tidy(r8outcome_mod, exponentiate = TRUE, conf.int = TRUE)

summary(r9outcome_mod <- glm(death2018 ~
                                   r9age_y_int + r8cesd + r8conde_impute,
                                 family = binomial(link = "logit"),
                                 data = hrs_samp))
r9outcome_results <- tidy(r9outcome_mod, exponentiate = TRUE, conf.int = TRUE)

results_tbl <- tibble(
  variables = c("Intercept", "age at current wave", "previous CESD value",
                "Previous chronic condition count"),
  r4beta = round(r4outcome_results$estimate, 4),
  r5beta = round(r5outcome_results$estimate, 4),
  r6beta = round(r6outcome_results$estimate, 4),
  r7beta = round(r7outcome_results$estimate, 4),
  r8beta = round(r8outcome_results$estimate, 4),
  r9beta = round(r9outcome_results$estimate, 4)
)

results_tbl %>%
  kbl(caption = "Exponentiated betas association with death2018 missing model (wave 4 - 9)") %>%
  kable_classic(full_width = F, html_font = "Arial")

write_csv(results_tbl, paste0(path_to_dropbox,
                              "/exposure_trajectories/data/",
                              "predict_death2018_betas.csv"))