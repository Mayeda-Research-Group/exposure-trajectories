#---- package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "broom", "survival")

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

#---- read in datasets ----
CESD_data_wide <- 
  read_csv(paste0(path_to_dropbox, 
                  "/exposure_trajectories/data/", 
                  "CESD_data_wide.csv"), 
           col_types = cols(HHIDPN = col_character())) %>% as.data.frame() 

beta_0_table <- read_csv(paste0(path_to_dropbox, "/exposure_trajectories/data/", 
                                "beta_0_table.csv"))

beta_mat <- read_csv(paste0(path_to_dropbox, "/exposure_trajectories/", 
                            "data/beta_mat.csv"))

#---- function ----
test_cc_analysis <- function(data, mechanism, mask_percent){
  #---- mask data ----
  complete_data <- 
    mask(data, mechanism, mask_percent, beta_0_table, beta_mat)
  
  complete_data %<>% 
    mutate("r4cesd_elevated" = ifelse(r4cesd >= 4, 1, 0), 
           "r9cesd_elevated" = ifelse(r9cesd >= 4, 1, 0))
  
  #indicate where to take averages
  indicator <- complete_data %>% 
    dplyr::select(paste0("r", seq(4, 9), "cesd")) %>% is.na() %>% 
    rowSums() > 2
  complete_data %<>% mutate("avg_indicator" = 1 - (1*indicator))
  
  complete_data[, "avg_cesd"] <- 
    rowMeans(complete_data %>% 
               dplyr::select(paste0("r", seq(4, 9), "cesd")), na.rm = TRUE) 
  complete_data[which(complete_data$avg_indicator == 0), "avg_cesd"] <- NA
  complete_data %<>% 
    mutate("avg_cesd_elevated" = ifelse(avg_cesd >= 4, 1, 0))
  
  complete_data[, "prop_elevated_cesd"] <- 
    rowMeans(complete_data %>% 
               dplyr::select(paste0("r", seq(4, 9), "cesd")) %>% 
               mutate_all(function(x) ifelse(x >= 4, 1, 0)), na.rm = TRUE) 
  complete_data[which(complete_data$avg_indicator == 0), 
                "prop_elevated_cesd"] <- NA
  
  #---- models ----
  model_list <- vector(mode = "list", length = length(exposures)) %>% 
    set_names(exposures)
  
  for(exposure in exposures){
    if(exposure == "CES-D Wave 4"){
      model_list[[exposure]] <- 
        with(complete_data, 
             coxph(Surv(survtime, observed) ~ r4not_married_partnered + 
                     r4widowed + ed_cat + r4drinking_cat + r4memrye_impute + 
                     r4stroke_impute + r4hearte_impute + r4lunge_impute + 
                     r4cancre_impute + r4hibpe_impute + r4diabe_impute + 
                     smoker + r4BMI + hispanic + black + other + female + 
                     r4age_y_int + r4cesd_elevated))
    } else if(exposure == "CES-D Wave 9"){
      model_list[[exposure]] <- 
        with(complete_data, 
             coxph(Surv(survtime, observed) ~ r9not_married_partnered + 
                     r9widowed + ed_cat + r9drinking_cat + r9memrye_impute + 
                     r9stroke_impute + r9hearte_impute + r9lunge_impute + 
                     r9cancre_impute + r9hibpe_impute + r9diabe_impute + 
                     smoker + r9BMI + hispanic + black + other + female + 
                     r9age_y_int + r9cesd_elevated))
    } else if(exposure == "Elevated CES-D Prop"){
      model_list[[exposure]] <- 
        with(complete_data, 
             coxph(Surv(survtime, observed) ~ r4not_married_partnered + 
                     r4widowed + ed_cat + r4drinking_cat + r4memrye_impute + 
                     r4stroke_impute + r4hearte_impute + r4lunge_impute + 
                     r4cancre_impute + r4hibpe_impute + r4diabe_impute + 
                     smoker + r4BMI + hispanic + black + other + female + 
                     r4age_y_int + prop_elevated_cesd))
      
    } else{
      model_list[[exposure]] <- 
        with(complete_data, 
             coxph(Surv(survtime, observed) ~ r4not_married_partnered + 
                     r4widowed + ed_cat + r4drinking_cat + r4memrye_impute + 
                     r4stroke_impute + r4hearte_impute + r4lunge_impute + 
                     r4cancre_impute + r4hibpe_impute + r4diabe_impute + 
                     smoker + r4BMI + hispanic + black + other + female + 
                     r4age_y_int + avg_cesd_elevated))
    }
  }
  
  #---- create shell for output ----
  pooled_effect_ests <- 
    data.frame("Exposure" = exposures, "beta" = NA, "SD" = NA, "LCI" = NA, 
               "UCI" = NA, "Missingness" = mask_percent, "Type" = mechanism)
  
  for(exposure in exposures){
    #---- store output ----
    pooled_model <- broom::tidy(model_list[[exposure]])
    
    pooled_effect_ests[which(pooled_effect_ests$Exposure == exposure), 
                       c("beta", "SD")] <- 
      pooled_model[nrow(pooled_model), c("estimate", "std.error")]
  }
  
  pooled_effect_ests[, "LCI"] <- 
    pooled_effect_ests$beta - 1.96*pooled_effect_ests$SD
  
  pooled_effect_ests[, "UCI"] <- 
    pooled_effect_ests$beta + 1.96*pooled_effect_ests$SD
  
  #---- **return ----
  return(pooled_effect_ests)
}

#---- run tests ----
exposures <- c("CES-D Wave 4", "CES-D Wave 9", "Elevated Average CES-D", 
               "Elevated CES-D Prop")

num_runs <- 1
mechanisms <- c("MCAR", "MAR", "MNAR")
percents <- c("10%", "20%", "30%")

for(mech in mechanisms){
  for(percent in percents){
    if(!exists("results")){
      results <- 
        replicate(num_runs, 
                  test_cc_analysis(data = CESD_data_wide, 
                                   mechanism = mech, mask_percent = percent), 
                  simplify = FALSE) %>% do.call(rbind, .) %>% 
        mutate("run" = rep(seq(1, num_runs), each = 4), .)
    } else{
      results <-
        rbind(results,  
              replicate(num_runs, 
                        test_cc_analysis(data = CESD_data_wide, 
                                         mechanism = mech, 
                                         mask_percent = percent), 
                        simplify = FALSE) %>% do.call(rbind, .) %>% 
                mutate("run" = rep(seq(1, num_runs), each = 4), .))
    }
  }
}


MCAR_10_avg <- MCAR_10 %>% group_by(Exposure) %>%
  summarize_at(.vars = c("beta", "SD", "LCI", "UCI"), .funs = mean)

#code for pivoting
test <- MCAR_10 %>% dplyr::select(c("run", "Exposure", "beta"))
pivot_wider(data = test, id_cols = run, names_from = Exposure, values_from = beta)



