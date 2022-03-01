#---- package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "broom", "ResourceSelection", 
       "survival", "openxlsx", "lubridate", "future.apply", "lme4", "devtools", 
       "miceFast")
devtools::install_github(repo = "amices/mice")
library(mice)

#No scientific notation
options(scipen = 999)

set.seed(20200819)

#---- note ----
# Since the difference between win and OS, put substituted directory here
# Yingyan's directory: C:/Users/yingyan_wu
#                      C:/Users/yingyan_wu/Dropbox
# Crystal's directory: /Users/CrystalShaw
#                     ~/Dropbox/Projects
# MRG desktop directory: C:/Users/cshaw/Dropbox/Projects

#Changing directories here will change them throughout the script
path_to_dropbox <- "~/Dropbox/Projects"

#---- source scripts ----
source(here::here("RScripts", "mask.R"))
source(here::here("RScripts", "mask_impute_pool.R"))
source(here::here("RScripts", "fast_impute.R"))

#---- read in analytical sample ----
CESD_data_wide <- 
  read_csv(paste0(path_to_dropbox, 
                  "/exposure_trajectories/data/", 
                  "CESD_data_wide.csv"), 
           col_types = cols(HHIDPN = col_character())) %>% as.data.frame() 

# #Check column types
# sapply(CESD_data_wide, class)

# #---- all combos ----
# methods <- c("JMVN")
# mechanisms <- c("MCAR", "MAR", "MNAR")
# mask_props <- c(.10, 0.20, 0.30)
# all_combos <- expand_grid(mechanisms, methods, mask_props) %>%
#   mutate("mask_percent" = paste0(100*mask_props, "%"))

# #---- create one set of imputations for plot ----
# start <- Sys.time()
# single_run <- apply(all_combos[2:nrow(all_combos), ], 1, function(x)
#   mask_impute_pool(data_wide = CESD_data_wide, exposures = exposures,
#                    mechanism = x["mechanisms"],
#                    method = x["methods"],
#                    mask_percent = x["mask_percent"], beta_0_table = beta_0_table, 
#                    beta_mat = beta_mat, truth = truth, save = "yes"))
# end <- Sys.time() - start

#---- create cluster ----
plan(multisession, gc = FALSE, workers = 10)

#---- get pooled effect estimates ----
start <- Sys.time()
for(i in which(!table_effect_ests$Method == "Truth")){
  #because we fill multiple rows at a time
  if(is.na(table_effect_ests[i, "beta"])){
    mechanism = table_effect_ests[i, "Type"]
    method = table_effect_ests[i, "Method"]
    mask_percent = table_effect_ests[i, "Missingness"]
    
    multi_runs <- 
      future_replicate(num_runs, 
                       mask_impute_pool(CESD_data_wide, exposures, 
                                        mechanism = mechanism, method = method, 
                                        mask_percent = mask_percent, 
                                        beta_mat = beta_mat, 
                                        beta_0_table = beta_0_table,
                                        truth = truth, save = "no"), 
                       simplify = FALSE)
    #Formatting data
    formatted <- do.call(rbind, multi_runs)
    
    #Storing results
    table_effect_ests[which(table_effect_ests$Method == method & 
                              table_effect_ests$Missingness == mask_percent & 
                              table_effect_ests$Type == mechanism), 
                      c("Exposure", "beta", "SE", "mean_LCI", "mean_UCI")] <- 
      formatted %>% group_by(Exposure) %>%
      summarize_at(.vars = c("beta", "SE", "LCI", "UCI"), .funs = mean)
    
    table_effect_ests[which(table_effect_ests$Method == method & 
                              table_effect_ests$Missingness == mask_percent & 
                              table_effect_ests$Type == mechanism), 
                      "LCI_beta"] <- 
      formatted %>% group_by(Exposure) %>%
      summarize_at(.vars = "beta", ~ quantile(.x, 0.025)) %>% 
      dplyr::select("beta")
    
    table_effect_ests[which(table_effect_ests$Method == method & 
                              table_effect_ests$Missingness == mask_percent & 
                              table_effect_ests$Type == mechanism), 
                      "UCI_beta"] <- 
      formatted %>% group_by(Exposure) %>%
      summarize_at(.vars = "beta", ~ quantile(.x, 0.975)) %>% 
      dplyr::select("beta")
    
    table_effect_ests[which(table_effect_ests$Method == method & 
                              table_effect_ests$Missingness == mask_percent & 
                              table_effect_ests$Type == mechanism), 
                      "truth_capture"] <- 
      formatted %>% group_by(Exposure) %>%
      summarize_at(.vars = c("capture_truth"), .funs = mean) %>% 
      dplyr::select("capture_truth")
  }
}
end <- Sys.time() - start

#---- shut down cluster ----
future::plan("sequential")

#---- save tables ----
# #Round numbers in dataframe
# table_effect_ests %<>% mutate(across(where(is.numeric), ~ round(., 3)))

#Save results 
write_csv(table_effect_ests, 
          file = paste0(path_to_dropbox,
                        "/exposure_trajectories/manuscript/",
                        "tables/results_", method, "_", num_runs, 
                        "_", format(now(), "%Y%m%d"),
                        ".csv"))
