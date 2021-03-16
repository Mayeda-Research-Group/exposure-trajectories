#---- package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "mice", "broom", "ResourceSelection", 
       "survival", "openxlsx", "lubridate")

#No scientific notation
options(scipen = 999)

set.seed(20200819)

#---- source scripts ----
source(here::here("RScripts", "mask_impute_pool.R"))

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
           col_types = cols(HHIDPN = col_character())) %>% 
  mutate_if(is.character, as.factor) 
#Will see if we want to use ordered factors
# %>% 
#   mutate_at(vars(c(paste0("r", seq(4, 9), "drinking_cat"), ed_cat)), 
#             ~factor(., ordered = TRUE)) %>%
#   mutate_at(vars(paste0("r", seq(4, 9), "drinking_cat")), 
#             ~fct_relevel(., c("No Drinking", "Moderate Drinking", 
#                               "Heavy Drinking"))) %>% 
#   mutate_at(vars("ed_cat"), 
#             ~fct_relevel(., c("Less than HS", "HS", "Some College", 
#                               "Bachelors")))

# #Check column types
# sapply(CESD_data_wide, class)

#---- Table 2 shell: Effect Estimates ----
#Number of simulation runs
num_runs <- 2
exposures <- c("CES-D Wave 4", "CES-D Wave 9", "Elevated Average CES-D", 
               "Elevated CES-D Count")
methods <- c("FCS", "JMVN", "PMM")
mechanisms <- c("MCAR")
mask_props <- c(.10, .20, .30)

table_effect_ests <- 
  data.frame(expand_grid(exposures, "Truth", mechanisms, "0%")) %>% 
  set_colnames(c("Exposure", "Method", "Type", "Missingness")) %>% 
  rbind(expand_grid(exposures, methods, mechanisms, 
                    paste0(mask_props*100, "%")) %>% 
          set_colnames(c("Exposure", "Method", "Type", "Missingness"))) %>% 
  mutate("beta" = NA, "SD" = NA, "mean_LCI" = NA, "mean_UCI" = NA, "LCI" = NA, 
         "UCI" = NA)

#---- truth ----
#---- **CES-D Wave 4 ----
TTEmodel_CESD4 <- 
  coxph(Surv(survtime, observed) ~ r4mstat_cat + ed_cat + r4drinking_cat + 
          r4memrye_impute + r4stroke_impute + r4hearte_impute + r4lunge_impute + 
          r4cancre_impute + r4hibpe_impute + r4diabe_impute + smoker + r4BMI + 
          hispanic + black + other + female + r4age_y_int + r4shlt + 
          r4cesd_elevated, data = CESD_data_wide)

#summary(TTEmodel_CESD4)

TTEmodel_CESD4_results <- tidy(TTEmodel_CESD4, 
                               exponentiate = FALSE, conf.int = TRUE)

table_effect_ests[which(table_effect_ests$Exposure == "CES-D Wave 4" & 
                          table_effect_ests$Method == "Truth"), 
                  c("beta", "SD", "mean_LCI", "mean_UCI", "LCI", "UCI")] <- 
  c(TTEmodel_CESD4_results[nrow(TTEmodel_CESD4_results), 
                           c("estimate", "std.error",  
                             rep(c("conf.low", "conf.high"), 2))])

#---- **CES-D Wave 9 ----
TTEmodel_CESD9 <- 
  coxph(Surv(survtime, observed) ~ r9mstat_cat + ed_cat + r9drinking_cat + 
          r9memrye_impute + r9stroke_impute + r9hearte_impute + r9lunge_impute + 
          r9cancre_impute + r9hibpe_impute + r9diabe_impute + smoker + r9BMI + 
          hispanic + black + other + female + r9age_y_int + r9shlt + 
          r9cesd_elevated, data = CESD_data_wide)

#summary(TTEmodel_CESD9)

TTEmodel_CESD9_results <- tidy(TTEmodel_CESD9, 
                               exponentiate = FALSE, conf.int = TRUE)

table_effect_ests[which(table_effect_ests$Exposure == "CES-D Wave 9" & 
                          table_effect_ests$Method == "Truth"), 
                  c("beta", "SD", "mean_LCI", "mean_UCI", "LCI", "UCI")] <- 
  c(TTEmodel_CESD9_results[nrow(TTEmodel_CESD9_results), 
                           c("estimate", "std.error", 
                             rep(c("conf.low", "conf.high"), 2))])

#---- **Total Count Elevated CES-D ----
TTEmodel_total_CESD <- 
  coxph(Surv(survtime, observed) ~ r4mstat_cat + ed_cat + r4drinking_cat + 
          r4memrye_impute + r4stroke_impute + r4hearte_impute + r4lunge_impute + 
          r4cancre_impute + r4hibpe_impute + r4diabe_impute + smoker + r4BMI + 
          hispanic + black + other + female + r4age_y_int + r4shlt + 
          total_elevated_cesd, data = CESD_data_wide)

TTEmodel_total_CESD_results <- tidy(TTEmodel_total_CESD, 
                                    exponentiate = FALSE, conf.int = TRUE)

table_effect_ests[which(table_effect_ests$Exposure == "Elevated CES-D Count" & 
                          table_effect_ests$Method == "Truth"), 
                  c("beta", "SD", "mean_LCI", "mean_UCI", "LCI", "UCI")] <- 
  c(TTEmodel_total_CESD_results[nrow(TTEmodel_total_CESD_results), 
                                c("estimate", "std.error", 
                                  rep(c("conf.low", "conf.high"), 2))])

#---- **Elevated Average CES-D ----
TTEmodel_elevated_avg_CESD <- 
  coxph(Surv(survtime, observed) ~ r4mstat_cat + ed_cat + r4drinking_cat + 
          r4memrye_impute + r4stroke_impute + r4hearte_impute + r4lunge_impute + 
          r4cancre_impute + r4hibpe_impute + r4diabe_impute + smoker + r4BMI + 
          hispanic + black + other + female + r4age_y_int + r4shlt + 
          avg_cesd_elevated, data = CESD_data_wide)

TTEmodel_elevated_avg_CESD_results <- tidy(TTEmodel_elevated_avg_CESD, 
                                           exponentiate = FALSE, 
                                           conf.int = TRUE)

table_effect_ests[which(table_effect_ests$Exposure == "Elevated Average CES-D" & 
                          table_effect_ests$Method == "Truth"), 
                  c("beta", "SD", "mean_LCI", "mean_UCI", "LCI", "UCI")] <- 
  c(TTEmodel_elevated_avg_CESD_results[nrow(TTEmodel_elevated_avg_CESD_results), 
                                       c("estimate", "std.error",
                                         rep(c("conf.low", "conf.high"), 2))])

#---- create one set of imputations for plot ----
all_combos <- expand_grid(mechanisms, methods, mask_props) %>%
  mutate("mask_percent" = paste0(100*mask_props, "%"))

single_run <- apply(all_combos[8:9, ], 1, function(x)
  mask_impute_pool(data_wide = CESD_data_wide, exposures = exposures, 
                   mechanism = x["mechanisms"], 
                   method = x["methods"], 
                   mask_percent = x["mask_percent"], num_impute = 5, 
                   save = "yes"))

#---- get pooled effect estimates ----
start <- Sys.time()
for(i in which(!table_effect_ests$Method == "Truth")){
  #because we fill multiple rows at a time
  if(is.na(table_effect_ests[i, "beta"])){
    mechanism = table_effect_ests[i, "Type"]
    method = table_effect_ests[i, "Method"]
    mask_percent = table_effect_ests[i, "Missingness"]
    
    multi_runs <- 
      replicate(num_runs, 
                mask_impute_pool(CESD_data_wide, exposures, 
                                 mechanism = mechanism, method = method, 
                                 mask_percent = mask_percent,
                                 num_impute = 5, save = "no"), 
                simplify = FALSE)
    
    #Formatting data
    formatted <- do.call(rbind, multi_runs)
    
    #Storing results
    table_effect_ests[which(table_effect_ests$Method == method & 
                              table_effect_ests$Missingness == mask_percent & 
                              table_effect_ests$Type == mechanism), 
                      c("Exposure", "beta", "SD", "mean_LCI", "mean_UCI")] <- 
      formatted %>% group_by(Exposure) %>%
      summarize_at(.vars = c("beta", "SD", "LCI", "UCI"), .funs = mean)
    
    table_effect_ests[which(table_effect_ests$Method == method & 
                              table_effect_ests$Missingness == mask_percent & 
                              table_effect_ests$Type == mechanism), "LCI"] <- 
      formatted %>% group_by(Exposure) %>%
      summarize_at(.vars = "beta", ~ quantile(.x, 0.025)) %>% 
      dplyr::select("beta")
    
    table_effect_ests[which(table_effect_ests$Method == method & 
                              table_effect_ests$Missingness == mask_percent & 
                              table_effect_ests$Type == mechanism), "UCI"] <- 
      formatted %>% group_by(Exposure) %>%
      summarize_at(.vars = "beta", ~ quantile(.x, 0.975)) %>% 
      dplyr::select("beta")
  }
}
end <- Sys.time() - start

#---- save tables ----
#Round numbers in dataframe
table_effect_ests %<>% mutate(across(where(is.numeric), ~ round(., 2)))

#Save results
table_list <- list("Table 2" = table_effect_ests)
write.xlsx(table_list, file = paste0(path_to_dropbox,
                                     "/exposure_trajectories/manuscript/",
                                     "tables/main_text_tables", num_runs, "_", 
                                     format(now(), "%Y%m%d_%H%M%S_"),
                                     ".xlsx"))









