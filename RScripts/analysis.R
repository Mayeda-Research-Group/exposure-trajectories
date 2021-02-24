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
exposures <- c("CES-D Wave 4", "CES-D Wave 9", "Elevated Average CES-D", 
               "Elevated CES-D Count")
methods <- c("JMVN", "FCS")
mechanisms <- c("MCAR", "MAR", "NMAR")
mask_props <- c(.10, .25, .50)

table_effect_ests <- 
  data.frame("Exposure" = rep(rep(exposures, 13), 3),
             "beta" = NA, "LCI" = NA, "UCI" = NA, 
             "Method" = rep(c(rep("Truth", 4), rep(methods, each = 12)), 3), 
             "Missingness" = rep(c(rep("0%", 4), 
                                   rep(rep(paste0(mask_props*100, "%"), 
                                           each = 4), 4)), 3), 
             "Type" = rep(mechanisms, each = 52))

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

apply(all_combos[1:6, ], 1, function(x)
  mask_impute_pool(CESD_data_wide, mechanism = x[1], method = x[2],
                   mask_percent = x[4],
                   num_impute = 5, save = "yes"))

#---- get pooled effect estimates ----
for(i in 1:6){
  mechanism = as.character(all_combos[i, "mechanisms"])
  method = as.character(all_combos[i, "methods"])
  mask_percent = as.character(all_combos[i, "mask_percent"])
  
  multi_runs <- 
    replicate(2, mask_impute_pool(CESD_data_wide, mechanism = mechanism, 
                                  method = method, mask_percent = mask_percent,
                                  num_impute = 5, save = "no"), 
              simplify = FALSE)
  
  #Formatting data
  formatted <- do.call(rbind, multi_runs)
  
  #Summarizing results
  results <- formatted %>% group_by(Exposure) %>%
    summarize_at(.vars = c("beta", "LCI", "UCI"), .funs = mean)
  
  #Storing results
  table_effect_ests[which(table_effect_ests$Method == method & 
                            table_effect_ests$Missingness == mask_percent & 
                            table_effect_ests$Type == mechanism), 
                    c("Exposure", "beta", "LCI", "UCI")] <- results
}

#---- save tables ----
#Round numbers in dataframe
table_effect_ests %<>% mutate(across(where(is.numeric), ~ round(., 2)))

#Save results
table_list <- list("Table 2" = table_effect_ests)
write.xlsx(table_list, file = paste0(path_to_dropbox, 
                                     "/exposure_trajectories/manuscript/", 
                                     "tables/main_text_tables.xlsx"))









