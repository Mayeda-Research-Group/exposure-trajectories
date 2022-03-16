#This script creates the truth tables for analytic models 
# (main and sensitivity analyses)
# 
# Input:CESD_data_wide.csv
# Output:truth.csv, truth_sens.csv

#---- package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "broom", "ResourceSelection", 
       "survival")

#---- note ----
# Since the difference between win and OS, put substituted directory here
# Yingyan's directory: C:/Users/yingyan_wu
#                      C:/Users/yingyan_wu/Dropbox
# Crystal's directory: /Users/CrystalShaw
#                     ~/Dropbox/Projects
# MRG desktop directory: C:/Users/cshaw/Dropbox/Projects

#Changing directories here will change them throughout the script
path_to_dropbox <- "~/Dropbox/Projects"

#---- read in analytical sample ----
CESD_data_wide <- 
  read_csv(paste0(path_to_dropbox, 
                  "/exposure_trajectories/data/", 
                  "CESD_data_wide.csv"), 
           col_types = cols(HHIDPN = col_character())) %>% as.data.frame() 

#---- truth table ----
exposures <- c("CES-D Wave 4", "CES-D Wave 9", "Elevated Average CES-D", 
               "Elevated CES-D Prop")

table_effect_ests <- 
  data.frame(expand_grid(exposures, "Truth")) %>% 
  set_colnames(c("Exposure", "Method")) %>% 
  mutate("beta" = NA, "SE" = NA, "LCI_beta" = NA, "UCI_beta" = NA)

table_effect_ests_sens <- 
  data.frame(expand_grid(exposures, "Truth")) %>% 
  set_colnames(c("Exposure", "Method")) %>% 
  mutate("beta" = NA, "SE" = NA, "LCI_beta" = NA, "UCI_beta" = NA)

#---- **CES-D Wave 4 ----
TTEmodel_CESD4 <- 
  coxph(Surv(survtime, observed) ~ r4not_married_partnered + r4widowed + 
          ed_cat + r4drinking_cat + r4memrye_impute + r4stroke_impute + 
          r4hearte_impute + r4lunge_impute + r4cancre_impute + r4hibpe_impute + 
          r4diabe_impute + smoker + r4BMI + hispanic + black + other + female + 
          r4age_y_int + r4cesd_elevated, data = CESD_data_wide)

#summary(TTEmodel_CESD4)

TTEmodel_CESD4_results <- tidy(TTEmodel_CESD4, 
                               exponentiate = FALSE, conf.int = TRUE)

table_effect_ests[which(table_effect_ests$Exposure == "CES-D Wave 4"), 
                  c("beta", "SE", "LCI_beta", "UCI_beta")] <- 
  c(TTEmodel_CESD4_results[nrow(TTEmodel_CESD4_results), 
                           c("estimate", "std.error", "conf.low", "conf.high")])

#---- **CES-D Wave 4 sens ----
TTEmodel_CESD4_sens <- 
  coxph(Surv(survtime, observed) ~ r4not_married_partnered + r4widowed + 
          ed_cat + r4drinking_cat + r4memrye_impute + r4stroke_impute + 
          r4hearte_impute + r4lunge_impute + r4cancre_impute + r4hibpe_impute + 
          r4diabe_impute + smoker + r4BMI + hispanic + black + other + female + 
          r4age_y_int + r4cesd_elevated_sens, data = CESD_data_wide)

#summary(TTEmodel_CESD4_sens)

TTEmodel_CESD4_sens_results <- tidy(TTEmodel_CESD4_sens, 
                                    exponentiate = FALSE, conf.int = TRUE)

table_effect_ests_sens[which(table_effect_ests_sens$Exposure == "CES-D Wave 4"), 
                       c("beta", "SE", "LCI_beta", "UCI_beta")] <- 
  c(TTEmodel_CESD4_sens_results[nrow(TTEmodel_CESD4_sens_results), 
                                c("estimate", "std.error", "conf.low", "conf.high")])

#---- **CES-D Wave 9 ----
TTEmodel_CESD9 <- 
  coxph(Surv(survtime, observed) ~ r9not_married_partnered + r9widowed + 
          ed_cat + r9drinking_cat + r9memrye_impute + r9stroke_impute + 
          r9hearte_impute + r9lunge_impute + r9cancre_impute + r9hibpe_impute + 
          r9diabe_impute + smoker + r9BMI + hispanic + black + other + female + 
          r9age_y_int + r9cesd_elevated, data = CESD_data_wide)

#summary(TTEmodel_CESD9)

TTEmodel_CESD9_results <- tidy(TTEmodel_CESD9, 
                               exponentiate = FALSE, conf.int = TRUE)

table_effect_ests[which(table_effect_ests$Exposure == "CES-D Wave 9"), 
                  c("beta", "SE", "LCI_beta", "UCI_beta")] <- 
  c(TTEmodel_CESD9_results[nrow(TTEmodel_CESD9_results), 
                           c("estimate", "std.error", "conf.low", "conf.high")])

#---- **CES-D Wave 9 sens ----
TTEmodel_CESD9_sens <- 
  coxph(Surv(survtime, observed) ~ r9not_married_partnered + r9widowed + 
          ed_cat + r9drinking_cat + r9memrye_impute + r9stroke_impute + 
          r9hearte_impute + r9lunge_impute + r9cancre_impute + r9hibpe_impute + 
          r9diabe_impute + smoker + r9BMI + hispanic + black + other + female + 
          r9age_y_int + r9cesd_elevated_sens, data = CESD_data_wide)

#summary(TTEmodel_CESD9_sens)

TTEmodel_CESD9_sens_results <- tidy(TTEmodel_CESD9_sens, 
                                    exponentiate = FALSE, conf.int = TRUE)

table_effect_ests_sens[which(table_effect_ests_sens$Exposure == "CES-D Wave 9"), 
                       c("beta", "SE", "LCI_beta", "UCI_beta")] <- 
  c(TTEmodel_CESD9_sens_results[nrow(TTEmodel_CESD9_sens_results), 
                                c("estimate", "std.error", "conf.low", "conf.high")])

#---- **Elevated Average CES-D ----
TTEmodel_elevated_avg_CESD <- 
  coxph(Surv(survtime, observed) ~ r4not_married_partnered + r4widowed + 
          ed_cat + r4drinking_cat + r4memrye_impute + r4stroke_impute + 
          r4hearte_impute + r4lunge_impute + r4cancre_impute + r4hibpe_impute + 
          r4diabe_impute + smoker + r4BMI + hispanic + black + other + female + 
          r4age_y_int + avg_cesd_elevated, data = CESD_data_wide)

TTEmodel_elevated_avg_CESD_results <- tidy(TTEmodel_elevated_avg_CESD, 
                                           exponentiate = FALSE, 
                                           conf.int = TRUE)

table_effect_ests[which(table_effect_ests$Exposure == "Elevated Average CES-D"), 
                  c("beta", "SE", "LCI_beta", "UCI_beta")] <- 
  c(TTEmodel_elevated_avg_CESD_results[nrow(TTEmodel_elevated_avg_CESD_results), 
                                       c("estimate", "std.error", "conf.low", 
                                         "conf.high")])

#---- **Elevated Average CES-D sens ----
TTEmodel_elevated_avg_CESD_sens <- 
  coxph(Surv(survtime, observed) ~ r4not_married_partnered + r4widowed + 
          ed_cat + r4drinking_cat + r4memrye_impute + r4stroke_impute + 
          r4hearte_impute + r4lunge_impute + r4cancre_impute + r4hibpe_impute + 
          r4diabe_impute + smoker + r4BMI + hispanic + black + other + female + 
          r4age_y_int + avg_cesd_elevated_sens, data = CESD_data_wide)

TTEmodel_elevated_avg_CESD_sens_results <- tidy(TTEmodel_elevated_avg_CESD_sens, 
                                                exponentiate = FALSE, 
                                                conf.int = TRUE)

table_effect_ests_sens[which(table_effect_ests_sens$Exposure == 
                               "Elevated Average CES-D"), 
                  c("beta", "SE", "LCI_beta", "UCI_beta")] <- 
  c(TTEmodel_elevated_avg_CESD_sens_results[
    nrow(TTEmodel_elevated_avg_CESD_sens_results), 
                                       c("estimate", "std.error", "conf.low", 
                                         "conf.high")])

#---- **Prop Elevated CES-D ----
TTEmodel_prop_CESD <- 
  coxph(Surv(survtime, observed) ~ r4not_married_partnered + r4widowed + 
          ed_cat + r4drinking_cat + r4memrye_impute + r4stroke_impute + 
          r4hearte_impute + r4lunge_impute + r4cancre_impute + r4hibpe_impute + 
          r4diabe_impute + smoker + r4BMI + hispanic + black + other + female + 
          r4age_y_int + prop_elevated_cesd, data = CESD_data_wide)

TTEmodel_prop_CESD_results <- tidy(TTEmodel_prop_CESD, 
                                   exponentiate = FALSE, conf.int = TRUE)

table_effect_ests[which(table_effect_ests$Exposure == "Elevated CES-D Prop"), 
                  c("beta", "SE", "LCI_beta", "UCI_beta")] <- 
  c(TTEmodel_prop_CESD_results[nrow(TTEmodel_prop_CESD_results), 
                               c("estimate", "std.error", "conf.low", 
                                 "conf.high")])

#---- **Prop Elevated CES-D sens ----
TTEmodel_prop_CESD_sens <- 
  coxph(Surv(survtime, observed) ~ r4not_married_partnered + r4widowed + 
          ed_cat + r4drinking_cat + r4memrye_impute + r4stroke_impute + 
          r4hearte_impute + r4lunge_impute + r4cancre_impute + r4hibpe_impute + 
          r4diabe_impute + smoker + r4BMI + hispanic + black + other + female + 
          r4age_y_int + prop_elevated_cesd_sens, data = CESD_data_wide)

TTEmodel_prop_CESD_sens_results <- tidy(TTEmodel_prop_CESD_sens, 
                                   exponentiate = FALSE, conf.int = TRUE)

table_effect_ests_sens[which(table_effect_ests_sens$Exposure == 
                               "Elevated CES-D Prop"), 
                  c("beta", "SE", "LCI_beta", "UCI_beta")] <- 
  c(TTEmodel_prop_CESD_sens_results[nrow(TTEmodel_prop_CESD_sens_results), 
                               c("estimate", "std.error", "conf.low", 
                                 "conf.high")])

#---- **save tables ----
write_csv(table_effect_ests, 
          paste0(path_to_dropbox, "/exposure_trajectories/data/truth.csv"))
write_csv(table_effect_ests_sens, 
          paste0(path_to_dropbox, "/exposure_trajectories/data/truth_sens.csv"))
