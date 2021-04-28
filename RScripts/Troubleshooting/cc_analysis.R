#---- package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "broom", "ResourceSelection", 
       "survival", "openxlsx", "lubridate", "lme4", "ghibli")

#No scientific notation
options(scipen = 999)

set.seed(20200819)

#---- source scripts ----
source(here::here("RScripts", "mask.R"))

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
           col_types = cols(HHIDPN = col_character())) 

# #Check column types
# sapply(CESD_data_wide, class)

#---- average missingess per wave ----
avg_miss <- function(data, mechanism, mask_percent){
  #---- mask data ----
  data_wide <- mask(data, mechanism, mask_percent)
  
  #---- percent missing per wave ----
  return(apply(data_wide[, paste0("r", seq(4, 9), "cesd")], 2, 
               function(x) mean(is.na(x))))
}

#---- **run sim ----
mechanisms <- c("MCAR", "MAR", "MNAR")
percents <- c("10%", "20%", "30%")
all_combos <- expand_grid(mechanisms, percents) 

for(combo in 1:nrow(all_combos)){
  mechanism = all_combos[[combo, "mechanisms"]]
  percent = all_combos[[combo, "percents"]]
  
  assign(paste0("results_", mechanism, percent), 
         rowMeans(replicate(100, avg_miss(CESD_data_wide, 
                                          mechanism = mechanism, 
                                          mask_percent = percent))))
}

#---- shell table ----
exposures <- c("CES-D Wave 4", "CES-D Wave 9", "Elevated Average CES-D", 
               "Elevated CES-D Count")

table_effect_ests <- 
  data.frame(expand_grid(exposures, "Truth", mechanisms, "0%")) %>% 
  set_colnames(c("Exposure", "Method", "Type", "Missingness")) %>% 
  rbind(expand_grid(exposures, "CC", mechanisms, 
                    paste0(mask_props*100, "%")) %>% 
          set_colnames(c("Exposure", "Method", "Type", "Missingness"))) %>% 
  mutate("beta" = NA, "SD" = NA, "mean_LCI" = NA, "mean_UCI" = NA, 
         "truth_capture" = NA, "people_dropped" = NA)

#---- truth ----
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

table_effect_ests[which(table_effect_ests$Exposure == "CES-D Wave 4" & 
                          table_effect_ests$Method == "Truth"), 
                  c("beta", "SD", "mean_LCI", "mean_UCI")] <- 
  c(TTEmodel_CESD4_results[nrow(TTEmodel_CESD4_results), 
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

table_effect_ests[which(table_effect_ests$Exposure == "CES-D Wave 9" & 
                          table_effect_ests$Method == "Truth"), 
                  c("beta", "SD", "mean_LCI", "mean_UCI")] <- 
  c(TTEmodel_CESD9_results[nrow(TTEmodel_CESD9_results), 
                           c("estimate", "std.error", "conf.low", "conf.high")])

#---- **Total Count Elevated CES-D ----
TTEmodel_total_CESD <- 
  coxph(Surv(survtime, observed) ~ r4not_married_partnered + r4widowed + 
          ed_cat + r4drinking_cat + r4memrye_impute + r4stroke_impute + 
          r4hearte_impute + r4lunge_impute + r4cancre_impute + r4hibpe_impute + 
          r4diabe_impute + smoker + r4BMI + hispanic + black + other + female + 
          r4age_y_int + total_elevated_cesd, data = CESD_data_wide)

TTEmodel_total_CESD_results <- tidy(TTEmodel_total_CESD, 
                                    exponentiate = FALSE, conf.int = TRUE)

table_effect_ests[which(table_effect_ests$Exposure == "Elevated CES-D Count" & 
                          table_effect_ests$Method == "Truth"), 
                  c("beta", "SD", "mean_LCI", "mean_UCI")] <- 
  c(TTEmodel_total_CESD_results[nrow(TTEmodel_total_CESD_results), 
                                c("estimate", "std.error", 
                                  "conf.low", "conf.high")])

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

table_effect_ests[which(table_effect_ests$Exposure == "Elevated Average CES-D" & 
                          table_effect_ests$Method == "Truth"), 
                  c("beta", "SD", "mean_LCI", "mean_UCI")] <- 
  c(TTEmodel_elevated_avg_CESD_results[nrow(TTEmodel_elevated_avg_CESD_results), 
                                       c("estimate", "std.error",
                                         "conf.low", "conf.high")])
                                         
truth <- table_effect_ests %>% filter(Method == "Truth", Type == "MCAR")

#---- cc analysis ----
cc <- function(data, mechanism, mask_percent, truth){
  cc_results <- 
    data.frame("Exposure" = exposures, "beta" = NA, "SD" = NA, "LCI" = NA, 
               "UCI" = NA, "Missingness" = mask_percent, 
               "Type" = mechanism, "capture_truth" = NA, "people_dropped" = NA)
  
  #---- mask data ----
  data_wide <- mask(data, mechanism, mask_percent)
  
  #---- **CES-D Wave 4 ----
  TTEmodel_CESD4 <- 
    coxph(Surv(survtime, observed) ~ r4not_married_partnered + r4widowed + 
            ed_cat + r4drinking_cat + r4memrye_impute + r4stroke_impute + 
            r4hearte_impute + r4lunge_impute + r4cancre_impute + 
            r4hibpe_impute + r4diabe_impute + smoker + r4BMI + hispanic + 
            black + other + female + r4age_y_int + r4cesd_elevated, 
          data = data_wide)
  
  cc_results[which(cc_results$Exposure == "CES-D Wave 4"), "people_dropped"] <- 
    (1 - TTEmodel_CESD4$n/nrow(data))
  
  TTEmodel_CESD4_results <- tidy(TTEmodel_CESD4, 
                                 exponentiate = FALSE, conf.int = TRUE)
  
  cc_results[which(cc_results$Exposure == "CES-D Wave 4"), 
                    c("beta", "SD", "LCI", "UCI")] <- 
    c(TTEmodel_CESD4_results[nrow(TTEmodel_CESD4_results), 
                             c("estimate", "std.error", 
                               "conf.low", "conf.high")])
  #---- **CES-D Wave 9 ----
  TTEmodel_CESD9 <- 
    coxph(Surv(survtime, observed) ~ r9not_married_partnered + r9widowed + 
            ed_cat + r9drinking_cat + r9memrye_impute + r9stroke_impute + 
            r9hearte_impute + r9lunge_impute + r9cancre_impute + r9hibpe_impute + 
            r9diabe_impute + smoker + r9BMI + hispanic + black + other + female + 
            r9age_y_int + r9cesd_elevated, data = data_wide)
  
  cc_results[which(cc_results$Exposure == "CES-D Wave 9"), "people_dropped"] <- 
    (1 - TTEmodel_CESD9$n/nrow(data))
  
  TTEmodel_CESD9_results <- tidy(TTEmodel_CESD9, 
                                 exponentiate = FALSE, conf.int = TRUE)
  
  cc_results[which(cc_results$Exposure == "CES-D Wave 9"), 
                    c("beta", "SD", "LCI", "UCI")] <- 
    c(TTEmodel_CESD9_results[nrow(TTEmodel_CESD9_results), 
                             c("estimate", "std.error", 
                               "conf.low", "conf.high")])
                               
  #---- **Total Count Elevated CES-D ----
  TTEmodel_total_CESD <- 
    coxph(Surv(survtime, observed) ~ r4not_married_partnered + r4widowed + 
            ed_cat + r4drinking_cat + r4memrye_impute + r4stroke_impute + 
            r4hearte_impute + r4lunge_impute + r4cancre_impute + r4hibpe_impute + 
            r4diabe_impute + smoker + r4BMI + hispanic + black + other + female + 
            r4age_y_int + total_elevated_cesd, data = data_wide)
  
  cc_results[which(cc_results$Exposure == "Elevated CES-D Count"), 
             "people_dropped"] <- (1 - TTEmodel_total_CESD$n/nrow(data))
  
  TTEmodel_total_CESD_results <- tidy(TTEmodel_total_CESD, 
                                      exponentiate = FALSE, conf.int = TRUE)
  
  cc_results[which(cc_results$Exposure == "Elevated CES-D Count"), 
                    c("beta", "SD", "LCI", "UCI")] <- 
    c(TTEmodel_total_CESD_results[nrow(TTEmodel_total_CESD_results), 
                                  c("estimate", "std.error", 
                                    "conf.low", "conf.high")])
                                    
  #---- **Elevated Average CES-D ----
  TTEmodel_elevated_avg_CESD <- 
    coxph(Surv(survtime, observed) ~ r4not_married_partnered + r4widowed + 
            ed_cat + r4drinking_cat + r4memrye_impute + r4stroke_impute + 
            r4hearte_impute + r4lunge_impute + r4cancre_impute + r4hibpe_impute + 
            r4diabe_impute + smoker + r4BMI + hispanic + black + other + female + 
            r4age_y_int + avg_cesd_elevated, data = data_wide)
  
  cc_results[which(cc_results$Exposure == "Elevated Average CES-D"), 
             "people_dropped"] <- (1 - TTEmodel_elevated_avg_CESD$n/nrow(data))
  
  TTEmodel_elevated_avg_CESD_results <- tidy(TTEmodel_elevated_avg_CESD, 
                                             exponentiate = FALSE, 
                                             conf.int = TRUE)
  
  cc_results[which(cc_results$Exposure == "Elevated Average CES-D"), 
                    c("beta", "SD", "LCI", "UCI")] <- 
    c(TTEmodel_elevated_avg_CESD_results[nrow(TTEmodel_elevated_avg_CESD_results), 
                                         c("estimate", "std.error", 
                                           "conf.low", "conf.high")])
  
  #---- truth capture ----
  cc_results$capture_truth <- 
    (truth$beta > cc_results$LCI)*(truth$beta < cc_results$UCI)
  
  #---- return ----
  return(cc_results)
}

#---- **run sim ----
start <- Sys.time()
all_combos <- expand_grid(mechanisms, percents)

for(combo in 1:nrow(all_combos)){
  mechanism = all_combos[[combo, "mechanisms"]]
  percent = all_combos[[combo, "percents"]]
  
  multi_runs <- replicate(100, cc(CESD_data_wide, mechanism, percent, truth), 
                          simplify = FALSE)
  #Formatting data
  formatted <- do.call(rbind, multi_runs)
  
  #Storing results
  table_effect_ests[which(table_effect_ests$Method == "CC" & 
                            table_effect_ests$Missingness == percent & 
                            table_effect_ests$Type == mechanism), 
                    c("Exposure", "beta", "SD", "mean_LCI", "mean_UCI", 
                      "truth_capture", "people_dropped")] <- 
    formatted %>% group_by(Exposure) %>%
    summarize_at(.vars = c("beta", "SD", "LCI", "UCI", "capture_truth", 
                           "people_dropped"), .funs = mean)
}
end <- Sys.time() - start

#save output
write_csv(table_effect_ests, 
          file = paste0(path_to_dropbox,
                        "/exposure_trajectories/manuscript/",
                        "tables/results_CC_100", format(now(), "%Y%m%d"), 
                        ".csv"))

#---- make figure ----
#or read in the data
results <- table_effect_ests

#---- **formatting the data ----
mask_props <- unique(results$Missingness)[-1] #Don't want 0%
results %<>% mutate_at(c("Missingness"), as.factor) 
results$Method <- factor(results$Method, levels = c("Truth", "CC"))
results$Type <- factor(results$Type, levels = c("MCAR", "MAR", "MNAR"))

#update exposure defs
results[which(results$Exposure == "CES-D Wave 4"), "Exposure"] <- 
  "Elevated CES-D Wave 4" 
results[which(results$Exposure == "CES-D Wave 9"), "Exposure"] <- 
  "Elevated CES-D Wave 9" 
results$Exposure <- 
  factor(results$Exposure, 
         levels = c("Elevated CES-D Wave 4", "Elevated CES-D Wave 9" , 
                    "Elevated Average CES-D", "Elevated CES-D Count"))

#---- **make plot ----
ggplot(results, 
       aes(x = beta, y = Missingness, color = Method, shape = Method)) +
  geom_point(size = 2.0, position = position_dodge(0.75)) + 
  scale_shape_manual(values = c(rep("square", (nrow(results))))) + 
  geom_errorbar(aes(xmin = mean_LCI, xmax = mean_UCI), width = .3, 
                position = position_dodge(0.75)) + theme_minimal() + 
  theme(legend.position = "bottom", legend.direction = "horizontal") + 
  scale_color_ghibli_d("LaputaMedium", direction = -1) + 
  scale_y_discrete(limits = rev(levels(results$Missingness))) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
  facet_grid(rows = vars(Type), cols = vars(Exposure)) + 
  ggtitle(paste0("Mean 95% CI of beta across 100 runs"))

ggsave(paste0(path_to_dropbox, "/exposure_trajectories/",
              "manuscript/figures/effect_ests_CC.jpeg"), 
       device = "jpeg", dpi = 300, width = 9, height = 7, units = "in")


