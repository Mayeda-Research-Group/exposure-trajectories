#---- package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "broom", "survival", "gcookbook", 
       "stringr")

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

#---- Sanity check full dataset analyses ----
full_dataset_analyses <- as.data.frame(matrix(nrow = 4, ncol = 7)) %>% 
  set_colnames(c("Exposure", "beta", "SE", "stat", "p-value", "LCI", "UCI"))

#---- **CES-D Wave 4 ----
TTEmodel_CESD4 <- 
  coxph(Surv(survtime, observed) ~ r4not_married_partnered + r4widowed + 
          ed_cat + r4drinking_cat + r4memrye_impute + r4stroke_impute + 
          r4hearte_impute + r4lunge_impute + r4cancre_impute + r4hibpe_impute + 
          r4diabe_impute + smoker + r4BMI + hispanic + black + other + female + 
          r4age_y_int + r4cesd_elevated, data = CESD_data_wide)

results <- tidy(TTEmodel_CESD4, exponentiate = FALSE, conf.int = TRUE)
full_dataset_analyses[1, ] <- results[nrow(results), ]

#---- **CES-D Wave 9 ----
TTEmodel_CESD9 <- 
  coxph(Surv(survtime, observed) ~ r9not_married_partnered + r9widowed + 
          ed_cat + r9drinking_cat + r9memrye_impute + r9stroke_impute + 
          r9hearte_impute + r9lunge_impute + r9cancre_impute + r9hibpe_impute + 
          r9diabe_impute + smoker + r9BMI + hispanic + black + other + female + 
          r9age_y_int + r9cesd_elevated, data = CESD_data_wide)

results <- tidy(TTEmodel_CESD9, exponentiate = FALSE, conf.int = TRUE)
full_dataset_analyses[2, ] <- results[nrow(results), ]

#---- **Prop Elevated CES-D ----
TTEmodel_prop_CESD <- 
  coxph(Surv(survtime, observed) ~ r4not_married_partnered + r4widowed + 
          ed_cat + r4drinking_cat + r4memrye_impute + r4stroke_impute + 
          r4hearte_impute + r4lunge_impute + r4cancre_impute + r4hibpe_impute + 
          r4diabe_impute + smoker + r4BMI + hispanic + black + other + female + 
          r4age_y_int + prop_elevated_cesd, data = CESD_data_wide)

results <- tidy(TTEmodel_prop_CESD, exponentiate = FALSE, conf.int = TRUE)
full_dataset_analyses[3, ] <- results[nrow(results), ]

#---- **Elevated Average CES-D ----
TTEmodel_elevated_avg_CESD <- 
  coxph(Surv(survtime, observed) ~ r4not_married_partnered + r4widowed + 
          ed_cat + r4drinking_cat + r4memrye_impute + r4stroke_impute + 
          r4hearte_impute + r4lunge_impute + r4cancre_impute + r4hibpe_impute + 
          r4diabe_impute + smoker + r4BMI + hispanic + black + other + female + 
          r4age_y_int + avg_cesd_elevated, data = CESD_data_wide)

results <- 
  tidy(TTEmodel_elevated_avg_CESD, exponentiate = FALSE, conf.int = TRUE)
full_dataset_analyses[4, ] <- results[nrow(results), ]

#---- **check table ----
full_dataset_analyses %<>% 
  mutate(Exposure = c("CES-D Wave 4", "CES-D Wave 9", "Elevated CES-D Prop", 
                      "Elevated Average CES-D"))

full_dataset_analyses

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
    data.frame("Exposure" = exposures, "beta" = NA, "SE" = NA, "LCI" = NA, 
               "UCI" = NA, "Missingness" = mask_percent, "Type" = mechanism)
  
  for(exposure in exposures){
    #---- store output ----
    pooled_model <- broom::tidy(model_list[[exposure]])
    
    pooled_effect_ests[which(pooled_effect_ests$Exposure == exposure), 
                       c("beta", "SE")] <- 
      pooled_model[nrow(pooled_model), c("estimate", "std.error")]
  }
  
  pooled_effect_ests[, "LCI"] <- 
    pooled_effect_ests$beta - 1.96*pooled_effect_ests$SE
  
  pooled_effect_ests[, "UCI"] <- 
    pooled_effect_ests$beta + 1.96*pooled_effect_ests$SE
  
  #---- **return ----
  return(pooled_effect_ests)
}

#---- run tests ----
exposures <- c("CES-D Wave 4", "CES-D Wave 9", "Elevated Average CES-D", 
               "Elevated CES-D Prop")

num_runs <- 5
mechanisms <- c("MCAR", "MAR", "MNAR")
percents <- paste0(seq(0, 30, by = 10), "%")
for(mech in mechanisms){
  for(percent in percents){
    if(!exists("test_results")){
      test_results <- 
        replicate(num_runs, 
                  test_cc_analysis(data = CESD_data_wide, 
                                   mechanism = mech, mask_percent = percent), 
                  simplify = FALSE) %>% do.call(rbind, .) %>% 
        mutate("run" = rep(seq(1, num_runs), each = 4), .)
    } else{
      test_results <-
        rbind(test_results,  
              replicate(num_runs, 
                        test_cc_analysis(data = CESD_data_wide, 
                                         mechanism = mech, 
                                         mask_percent = percent), 
                        simplify = FALSE) %>% do.call(rbind, .) %>% 
                mutate("run" = rep(seq(1, num_runs), each = 4), .))
    }
  }
}

write_csv(test_results, 
          paste0(path_to_dropbox, 
                 "/exposure_trajectories/inducing_missingness_troubleshooting/", 
                 "CC_analysis_", num_runs, "_runs.csv"))

# SE_summary <- test_results %>% group_by(Exposure, Type) %>%
#   summarize_at(.vars = c("SE"), .funs = c(min, max)) %>% 
#   set_colnames(c("Exposure", "Type", "min_SE_full", "max_SE_full"))

average_results <- test_results %>% group_by(Exposure, Type, Missingness) %>%
  summarize_at(.vars = c("beta", "SE", "LCI", "UCI"), .funs = mean) %>% 
  left_join(., full_dataset_analyses %>% dplyr::select(c(Exposure, SE)) %>% 
              set_colnames(c("Exposure", "SE_full")), by = "Exposure") %>% 
  ungroup()

# %>% 
#   mutate("xmin" = 
#            case_when(Exposure %in% c("CES-D Wave 4", "CES-D Wave 9") ~ 0.10, 
#                      TRUE ~ 0.02),
#          "xmax" = 
#            case_when(Exposure %in% c("CES-D Wave 4", "CES-D Wave 9") ~ 0.30, 
#                      TRUE ~ 0.26)) %>% 
#   left_join(., SE_summary, by = c("Exposure", "Type")) 

average_results %<>% 
  mutate("xaxis" = 
           as.numeric(str_sub(average_results$Missingness, end = -2))/100) 

#---- **plot ----
for(exposure in unique(average_results$Exposure)){
  for(type in unique(average_results$Type)){
    plot_data <- average_results %>% filter(Exposure == exposure & Type == type)
    
    ggplot(data = plot_data, aes(x = xaxis, y = SE)) + 
      geom_point() + geom_line() + 
      geom_hline(yintercept = unique(plot_data$SE_full), linetype = "dashed") +
      theme_minimal() + 
      ggtitle(paste0(unique(plot_data$Exposure), " | ", 
                     unique(plot_data$Type), " Missingness")) +
      # annotate("rect", xmin = unique(plot_data$xmin), 
      #          xmax = unique(plot_data$xmax), 
      #          ymin = unique(plot_data$min_SE_full), 
      #          ymax = unique(plot_data$max_SE_full), 
      #          alpha = .25 , fill = "blue") + 
      scale_x_continuous("Missing Percent", breaks = c(seq(0, 0.30, by = 0.10)), 
                         labels = paste0(seq(0, 30, by = 10), "%"))
    
    ggsave(filename = 
             paste0(path_to_dropbox, "/exposure_trajectories/", 
                    "inducing_missingness_troubleshooting/", 
                    "CC_analysis_plots/CC_analysis_", type, "_", exposure, 
                    ".jpeg"), width = 6, height = 4, units = "in")
  }
}

ggplot(data = average_results, aes(x = xaxis, y = SE)) + 
  geom_point() + geom_line() + 
  geom_hline(aes(yintercept = SE_full), linetype = "dashed") +
  theme_minimal() + 
  ggtitle("SE by percent missing, exposure, missingness mechanism") +
  scale_x_continuous("Missing Percent", breaks = c(seq(0, 0.30, by = 0.10)), 
                     labels = paste0(seq(0, 30, by = 10), "%")) + 
  facet_grid(rows = vars(Type), cols = vars(Exposure))

ggsave(filename = 
         paste0(path_to_dropbox, "/exposure_trajectories/", 
                "inducing_missingness_troubleshooting/", 
                "CC_analysis_plots/CC_analysis_", type, "_", exposure, 
                ".jpeg"), width = 6, height = 4, units = "in")

#---- OLD ----
#---- check masking ----
total_indices <- nrow(CESD_data_wide) 

#---- **plot data shell ----
mask_props <- seq(0.10, 0.80, by = 0.1)
exposures <- c("CES-D Wave 4", "CES-D Wave 9", "Elevated Average CES-D", 
               "Elevated CES-D Prop")

plot_data <- expand_grid(exposures, mask_props) %>% 
  mutate("SE" = NA)

plot_data <- data.frame("mask_prop" = seq(0.10, 0.80, by = 0.1), 
                        "SE" = NA)

for(prop in plot_data$mask_prop){
  mask_index <- sample(seq(1, total_indices), size = floor(prop*total_indices), 
                       replace = FALSE)
  masked_data <- CESD_data_wide
  masked_data[mask_index, "r4cesd"] <- NA
  masked_data %<>% mutate(r4cesd_elevated = ifelse(r4cesd >= 4, 1, 0))
  
  TTEmodel_CESD4 <- 
    coxph(Surv(survtime, observed) ~ r4not_married_partnered + r4widowed + 
            ed_cat + r4drinking_cat + r4memrye_impute + r4stroke_impute + 
            r4hearte_impute + r4lunge_impute + r4cancre_impute + 
            r4hibpe_impute + r4diabe_impute + smoker + r4BMI + hispanic + 
            black + other + female + r4age_y_int + r4cesd_elevated, 
          data = masked_data)
  
  results <- 
    tidy(TTEmodel_CESD4, exponentiate = FALSE, conf.int = TRUE)
  
  plot_data[which(plot_data$mask_prop == prop), "SE"] <- 
    unlist(results[nrow(results), "std.error"]) 
}

#---- **plot results ----
ggplot(data = plot_data, aes(x = mask_prop, y = SE)) + 
  geom_point() + geom_line() + 
  geom_hline(yintercept = 
               full_dataset_analyses[which(full_dataset_analyses$Exposure == 
                                             "r4cesd_elevated"), "SE"], 
             linetype = "dashed") +
  theme_minimal() + xlab("masking proportion") + 
  ggtitle("Elevated Baseline CES-D") +
  annotate("rect", xmin = 0.10, xmax = 0.30, 
           ymin = min(plot_data$SE), ymax = max(plot_data$SE), 
           alpha = .1 ,fill = "blue")

