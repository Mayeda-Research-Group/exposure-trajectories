#---- package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "data.table")

#No scientific notation
options(scipen = 999)

#---- note ----
# Since the difference between win and OS, put substituted directory here
# Yingyan's directory: C:/Users/yingyan_wu
#                      C:/Users/yingyan_wu/Dropbox
# Crystal's directory: /Users/CrystalShaw
#                     ~/Dropbox/Projects

#Changing directories here will change them throughout the script
path_to_box <- "/Users/CrystalShaw"
path_to_dropbox <- "~/Dropbox/Projects"

#---- color palette ----
# The palette with grey:
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#FFD700", "#0072B2", 
               "#D55E00", "#CC79A7")

#---- Figure 2: results ----
#---- **get filepaths ----
all_paths <- 
  list.files(path = paste0(path_to_dropbox,
                           "/exposure_trajectories/data/hoffman_transfer/",
                           "results"), full.names = TRUE, pattern = "*.csv")

main_paths <- all_paths[!str_detect(all_paths, "sens")]
sens_paths <- all_paths[str_detect(all_paths, "sens")]

#---- **read in data ----
read_results <- function(paths){
  data.table::fread(paths, fill = TRUE) %>% na.omit() %>%
    set_colnames(c("Exposure", "Beta", "SE", "LCI", "UCI", "Method",
                   "Percent", "Mechanism", "Truth Capture", "Time"))
}

# test <- data.table::fread(main_paths[2], fill = TRUE) %>%
#   set_colnames(c("Exposure", "Beta", "SE", "LCI", "UCI", "Method",
#                  "Percent", "Mechanism", "Truth Capture", "Time")) %>% na.omit()

main_results <- do.call(rbind, lapply(main_paths, read_results)) %>% na.omit()
sens_analyses <- do.call(rbind, lapply(sens_paths, read_results)) %>% na.omit()

#test for invalid rows
colSums(is.na(main_results))
colSums(is.na(sens_analyses))

#---- **take first 1000 runs (in case of extra) ----
main_results %<>% 
  group_by(Method, Mechanism, Percent, Exposure) %>% slice_head(n = 1000) %>% 
  na.omit()

sens_analyses %<>% 
  group_by(Method, Mechanism, Percent, Exposure) %>% slice_head(n = 1000) %>% 
  na.omit()

#---- **check scenario counts ----
#should be num_runs*num_exposures = 1000*4 = 4000 in each cell
table(main_results$Mechanism, main_results$Percent, main_results$Method)
table(sens_analyses$Mechanism, sens_analyses$Percent, sens_analyses$Method)

#---- **average run time ----
main_run_times <- main_results %>% 
  group_by(Method) %>% summarize_at(.vars = c("Time"), ~mean(., na.rm = TRUE)) 


#---- NEED TO EDIT DATAFRAME NAMES ----
#---- **summarize data ----
results_summary <- main_results %>%
  summarize_at(.vars = c("Beta", "SE", "LCI", "UCI", "Truth Capture"), 
               ~mean(., na.rm = TRUE)) 

# Sanity check
# results_sum_test <- 
#   results %>% group_by(Method, Mechanism, Percent, Exposure) %>%
#   summarise_all(list(mean))
# 
# diffdf::diffdf(results_summary, results_sum_test)

#---- I've remade the truth table with fewer columns ----
#---- **read in truth table ----
truth <- read_csv(paste0(path_to_dropbox, 
                         "/exposure_trajectories/data/", "truth.csv")) %>%
  select(-c(LCI_beta, UCI_beta)) %>%
  dplyr::rename(
    "Mechanism" = "Type",    
    "Percent" = "Missingness",
    "Beta" = "beta",
    "LCI" = "mean_LCI",
    "UCI" = "mean_UCI",
    "Truth Capture" = "truth_capture")

truth_multiple <- do.call("rbind", replicate(
  3, truth, simplify = FALSE)) %>%
  mutate(Mechanism = c(rep("MCAR", length(unique(truth$Exposure))), 
                       rep("MAR", length(unique(truth$Exposure))), 
                       rep("MNAR", length(unique(truth$Exposure)))))
results_summary %<>% rbind(truth_multiple)

#---- **format data ----
results_summary$Method <- 
  factor(results_summary$Method, levels = c("Truth", methods))
results_summary$Percent <- factor(results_summary$Percent)
results_summary$Mechanism <- 
  factor(results_summary$Mechanism, levels = c("MCAR", "MAR", "MNAR"))

#---- **plot ----
ggplot(results_summary, 
       aes(x = Beta, y = Percent, color = Method, shape = Method)) +
  geom_point(size = 2.0, position = position_dodge(0.75)) + 
  scale_shape_manual(values = c(rep("square", (nrow(results_summary))))) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = .3,
                position = position_dodge(0.75)) +
  theme_minimal() + 
  theme(legend.position = "bottom", legend.direction = "horizontal") + 
  scale_color_manual(values = cbPalette) + 
  scale_y_discrete(limits = rev(levels(results_summary$Percent))) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
  facet_grid(rows = vars(Mechanism), cols = vars(Exposure)) + 
  ggtitle(paste0("Mean 95% CI of beta across 100 runs"))

ggsave(paste0(path_to_dropbox, "/exposure_trajectories/",
              "manuscript/figures/effect_ests_mean_CI.jpeg"), 
       device = "jpeg", dpi = 300, width = 9, height = 7, units = "in")

#---- Supplement Figure 1: traceplots ----
#---- **read in data ----
test <- readRDS(here::here("MI datasets", "jmvn_mcar10"))
jpeg(paste0(path_to_dropbox, "/exposure_trajectories/",
            "manuscript/figures/trace_jmvn_mcar10.jpeg"), 
     width = 8, height = 10, units = "in", res = 300)
plot(test, r4cesd + r5cesd + r6cesd + r7cesd + r8cesd + r9cesd ~ 
       .it | .ms, layout = c(2, 6))
dev.off()

test <- readRDS(here::here("MI datasets", "jmvn_mnar30"))
jpeg(paste0(path_to_dropbox, "/exposure_trajectories/",
            "manuscript/figures/trace_jmvn_mnar30.jpeg"), 
     width = 8, height = 10, units = "in", res = 300)
plot(test, r4cesd + r5cesd + r6cesd + r7cesd + r8cesd + r9cesd ~ 
       .it | .ms, layout = c(2, 6))
dev.off()

plot(test, r4married_partnered + r5married_partnered + r6married_partnered + 
       r7married_partnered + r8married_partnered + r9married_partnered ~ 
       .it | .ms, layout = c(2, 6))



#---- OLD ----
# for(method in methods){
#   file_paths <-
#     list.files(path = paste0(path_to_dropbox,
#                              "/exposure_trajectories/data/hoffman_transfer/",
#                              "results/main/", method), full.names = TRUE,
#                pattern = "*.csv")
#   
#   if(!exists("main_results")){
#     main_results <- 
#       vroom(file_paths, col_names = FALSE) 
#   } else{
#     main_results %<>% rbind(vroom(file_paths, col_names = FALSE))
#   }
# }

# main_results %<>% 
#   set_colnames(c("Exposure", "Beta", "SE", "LCI", "UCI", "Method", "Percent", 
#                  "Mechanism", "Truth Capture", "Time"))

# for(method in methods){
#   file_paths <-
#     list.files(path = paste0(path_to_dropbox,
#                              "/exposure_trajectories/data/hoffman_transfer/",
#                              "results/main/", method), full.names = TRUE,
#                pattern = "*.csv")
#   
#   if(!exists("main_results")){
#     main_results <- do.call(rbind.data.frame, lapply(file_paths, read_results))
#     
#   } else{
#     main_results %<>% rbind(
#       do.call(rbind.data.frame, lapply(file_paths, read_results)))
#   }
# }

#---- read in data ----
#---- **CESD data ----
CESD_data_wide <- 
  read_csv(paste0(path_to_dropbox, 
                  "/exposure_trajectories/data/", 
                  "CESD_data_wide.csv"), 
           col_types = cols(HHIDPN = col_character())) %>% 
  mutate_if(is.character, as.factor) 

#---- **results ----
methods <- c("JMVN", "PMM", "FCS")
num_runs <- 10

#complete case
results <- 
  read_csv(paste0(path_to_dropbox, 
                  "/exposure_trajectories/manuscript/tables/", 
                  "complete_case_analyses/results_JMVN_100_20210411_cc.csv"))
results[which(results$Method != "Truth"), "Method"] <- "Complete Case"

for(method in methods){
  results <- 
    rbind(results, read_csv(Sys.glob(
      paste0(path_to_dropbox, "/exposure_trajectories/manuscript/tables/", 
             "results_", method, "_", num_runs, "*.csv"))) %>% 
        filter(Method != "Truth"))
}

#---- visualizations ----
#---- **effect estimates ----
mask_props <- unique(results$Missingness)[-1] #Don't want 0%
results %<>% mutate_at(c("Missingness"), as.factor) 
results$Method <- factor(results$Method, 
                         levels = c("Truth", "Complete Case", methods))
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

# #---- ****Distribution of beta ----
# ggplot(results, 
#        aes(x = beta, y = Missingness, color = Method, shape = Method)) +
#   geom_point(size = 2, position = position_dodge(0.60)) + 
#   scale_shape_manual(values = c(rep("square", (nrow(results))))) + 
#   geom_errorbar(aes(xmin = LCI_beta, xmax = UCI_beta), width = .2, 
#                 position = position_dodge(0.60)) + theme_minimal() + 
#   theme(legend.position = "bottom", legend.direction = "horizontal") + 
#   scale_color_ghibli_d("LaputaMedium", direction = -1) + 
#   scale_y_discrete(limits = rev(levels(results$Missingness))) + 
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
#   facet_grid(rows = vars(Type), cols = vars(Exposure)) + 
#   ggtitle(paste0("95% CI of beta across ", num_runs, " runs"))
# 
# ggsave(paste0(path_to_dropbox, "/exposure_trajectories/",
#               "manuscript/figures/effect_ests_dist_beta.jpeg"), device = "jpeg", 
#        dpi = 300, width = 9, height = 5, units = "in")

#---- ****mean LCI and mean UCI ----
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
  ggtitle(paste0("Mean 95% CI of beta across ", num_runs, " runs"))

ggsave(paste0(path_to_dropbox, "/exposure_trajectories/",
              "manuscript/figures/effect_ests_mean_LCI_UCI.jpeg"), 
       device = "jpeg", dpi = 300, width = 9, height = 7, units = "in")

#---- ****CI with mean SD ----
ggplot(results, 
       aes(x = beta, y = Missingness, color = Method, shape = Method)) +
  geom_point(size = 2, position = position_dodge(0.50)) + 
  scale_shape_manual(values = c(rep("square", (nrow(results))))) + 
  geom_errorbar(aes(xmin = beta - 1.96*SD, xmax = beta + 1.96*SD), width = .2, 
                position = position_dodge(0.50)) + theme_minimal() + 
  theme(legend.position = "bottom", legend.direction = "horizontal") + 
  scale_color_ghibli_d("LaputaMedium", direction = -1) + 
  scale_y_discrete(limits = rev(levels(results$Missingness))) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
  facet_grid(rows = vars(Type), cols = vars(Exposure)) + 
  ggtitle(paste0("95% CI of beta using mean SD across ", num_runs, " runs"))

ggsave(paste0(path_to_dropbox, "/exposure_trajectories/",
              "manuscript/figures/effect_ests_CI_mean_SD.jpeg"), 
       device = "jpeg", dpi = 300, width = 9, height = 5, units = "in")

#---- **individual imputations ----
#Read in data
methods <- c("jmvn_", "pmm_", "fcs_")
type <- c("mcar", "mar", "mnar")
mask_percent <- c("10", "20", "30")

filenames <- expand_grid(methods, type, mask_percent) %>% 
  unite("filenames", sep = "", remove = FALSE) 

for(name in filenames$filenames){
  assign(noquote(name), readRDS(here::here("MI datasets", name)))
}

#---- ****obs vs. pred CESD ----
for(i in 1:nrow(filenames)){
  for(wave in 4:9){
    imputed <- get(unlist(filenames[i, "filenames"]))[["imp"]][[
      paste0("r", wave, "cesd")]] %>% rowMeans()
    observed <- CESD_data_wide[as.integer(names(imputed)), 
                               paste0("r", wave, "cesd")]
    if(i == 1){
      obs_vs_pred_data <- cbind(observed, imputed) %>% as.data.frame() %>% 
        set_colnames(c("Observed", "Imputed")) %>% 
        mutate("Missingness" = paste0(unlist(filenames[i, "mask_percent"]), "%"), 
               "Type" = toupper(unlist(filenames[i, "type"])), 
               "Method" = toupper(str_sub(unlist(filenames[i, "methods"]), 
                                          end = -2)))
    } else{
      obs_vs_pred_data %<>% 
        rbind(., 
              cbind(observed, imputed) %>% as.data.frame() %>% 
                set_colnames(c("Observed", "Imputed")) %>% 
                mutate("Missingness" = 
                         paste0(unlist(filenames[i, "mask_percent"]), "%"), 
                       "Type" = toupper(unlist(filenames[i, "type"])), 
                       "Method" = 
                         toupper(str_sub(unlist(filenames[i, "methods"]), 
                                         end = -2))))
    }
  }
}

obs_vs_pred_CESD_plot <- 
  ggplot(data = obs_vs_pred_data, 
         aes(x = Observed, y = Imputed, color = Missingness)) + 
  geom_point(alpha = 0.50) + theme_bw() + geom_smooth(method = lm, se = FALSE) +
  geom_abline(slope = 1, intercept = 0, color = "black", lty = "dashed", 
              size = 1) + 
  facet_grid(rows = vars(Method), cols = vars(Type)) + 
  scale_color_manual(values = rev(ghibli_palette("LaputaMedium"))) + 
  ggtitle(paste0("Observed vs. Imputed CES-D scores (wave-aggregated)"))

ggsave(paste0("/Users/CrystalShaw/Dropbox/Projects/exposure_trajectories/",
              "manuscript/figures/obs_vs_pred_CESD.jpeg"), 
       obs_vs_pred_CESD_plot, device = "jpeg", width = 7, height = 4.5, 
       units = "in", dpi = 300)

#---- ****obs vs. pred exposures ----

#---- ****diagnostics: trace plots ----
pmm_mcar20_cesd_chainMean_data_raw <- 
  pmm_mcar20$chainMean[paste0("r", seq(4, 9), "cesd"), , ] 

for(chain in 1:dim(pmm_mcar20_cesd_chainMean_data_raw)[3]){
  if(chain == 1){
    pmm_mcar20_cesd_chainMean_data <- 
      pmm_mcar20_cesd_chainMean_data_raw[, , paste0("Chain ", chain)]
  } else{
    pmm_mcar20_cesd_chainMean_data <-
      rbind(pmm_mcar20_cesd_chainMean_data, 
            pmm_mcar20_cesd_chainMean_data_raw[, , paste0("Chain ", chain)])
  }
}

pmm_mcar20_cesd_chainMean_data %<>% as.data.frame() %>% 
  mutate("Measure" = rep(paste0("r", seq(4, 9), "cesd"), 
                         dim(pmm_mcar20_cesd_chainMean_data_raw)[3])) %>%
  mutate("Imputation" = rep(seq(1, dim(pmm_mcar20_cesd_chainMean_data_raw)[3]), 
                            each = 6))

#plot
gamma_chain_plot <- 
  ggplot(data = gamma_plot_data, 
         aes(x = reorder(run, sort(as.numeric(run))), y = gamma, 
             colour = Predictor)) +       
  geom_line(aes(group = Predictor)) + 
  facet_grid(rows = vars(factor(Group, 
                                levels = c("Unimpaired", "MCI", "Other")))) + 
  theme_bw() + xlab("Run") + 
  scale_color_manual(values = rev(extended_pallette14))




#---- ***fcs ----
fcs_mean_imputation_10 <- vector(mode = "list", length = 6)
for(i in 1:length(fcs_mean_imputation_10)){
  wave = i + 3
  fcs_mean_imputation_10[[i]] = 
    rowMeans((as.data.frame(fcs_mcar10$imp[[c(paste0("r", wave, "cesd"))]])) %>% 
               mutate_all(as.numeric))
}

#Do we want to stratify by wave?-- aggregate for now
num_missing = 0
for(i in 1:length(fcs_mean_imputation_10)){
  num_missing = num_missing + length(fcs_mean_imputation_10[[i]])
}

plot_data <- as.data.frame(matrix(nrow = num_missing, ncol = 2)) %>% 
  set_colnames(c("Observed", "Imputed"))
plot_data[, "Imputed"] <- unlist(fcs_mean_imputation_10)

for(i in 4:9){
  index <- as.numeric(rownames(fcs_mcar10[["imp"]][[paste0("r", i, "cesd")]]))
  start <- min(which(is.na(plot_data$Observed)))
  
  plot_data[start:(start + length(index) - 1), "Observed"] <- 
    CESD_data_wide[index, paste0("r", i, "cesd")]
}

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
              "manuscript/figures/mcar10_fcs_obs_pred.jpeg"), 
       device = "jpeg", width = 7, height = 4.5, units = "in", dpi = 300)

#---- **derived exposures ----





mean_imputation_mcar <- 
  lapply(mean_imputation_mcar <- vector(mode = "list", length(mask_props)),
         function(x) x <- lapply(x <- vector(mode = "list", 4),
                                 function(x) x <- vector(mode = "list", 
                                                         length = 6)))
#naming layers of list
names(mean_imputation_mcar) <- mask_props*100
for(i in 1:length(mean_imputation_mcar)){
  names(mean_imputation_mcar[[i]]) <- c("jmvn", "fcs", "jmvn_long", "fcs_long")
}

for(i in 1:length(mean_imputation_mcar)){
  data <- get(paste0("jmvn_mcar", mask_props[i]*100))
  for(j in 1:6){
    wave = j + 3
    mean_imputation_mcar[[i]][["jmvn"]][[j]] <- 
      rowMeans(data$imp[[c(paste0("logr", wave, "cesd"))]])
  }
}

plot_data <- as.data.frame(matrix(nrow = num_missing, ncol = 2)) %>% 
  set_colnames(c("Observed", "Imputed"))
plot_data[, "Imputed"] <- unlist(mean_imputation)

observed_data <- 
  cbind(CESD_data_wide %>% 
          dplyr::select(paste0("r", seq(4, 9), "cesd")) %>% 
          pivot_longer(everything(), names_to = "rwave", 
                       values_to = "complete"), 
        mcar10 %>% 
          dplyr::select(paste0("logr", seq(4, 9), "cesd")) %>% 
          pivot_longer(everything(), names_to = "logrwave", 
                       values_to = "masked")) %>% 
  arrange(rwave) %>% filter(is.na(masked))

plot_data[, "Observed"] <- observed_data$complete
plot_data %<>% mutate("log_observed" = log(1 + Observed), 
                      "trans_imputed" = round(exp(Imputed) - 1))

# #Sanity check
# head(CESD_data_wide$r4cesd[as.numeric(names(mean_imputation[[1]]))])
# head(plot_data)
# tail(CESD_data_wide$r9cesd[as.numeric(names(mean_imputation[[6]]))])
# tail(plot_data)

#plot
ggplot(data = plot_data, aes(x = log_observed, y = Imputed)) + 
  geom_point(color = "#B4DAE5FF") + 
  geom_smooth(method = lm, se = FALSE, color = "#F0D77BFF") +
  geom_abline(slope = 1, intercept = 0, color = "#5C5992FF", lty = "dashed", 
              size = 1) + labs(x = "Log(1 + Observed)") + 
  ggtitle(paste0("Missingness Pattern: MCAR 10% \n", 
                 "Imputation Strategy: Joint Multivariate Normal")) + 
  theme_minimal()  

ggsave(paste0("/Users/CrystalShaw/Dropbox/Projects/exposure_trajectories/",
              "manuscript/figures/mcar10_jmvn_obs_pred_log.jpeg"), 
       device = "jpeg", width = 7, height = 4.5, units = "in", dpi = 300)

ggplot(data = plot_data, aes(x = Observed, y = trans_imputed)) + 
  geom_point(color = "#B4DAE5FF") + 
  geom_smooth(method = lm, se = FALSE, color = "#F0D77BFF") +
  geom_abline(slope = 1, intercept = 0, color = "#5C5992FF", lty = "dashed", 
              size = 1) + labs(y = "Exp(Imputed) - 1") + 
  ggtitle(paste0("Missingness Pattern: MCAR 10% \n", 
                 "Imputation Strategy: Joint Multivariate Normal")) + 
  theme_minimal()  

ggsave(paste0("/Users/CrystalShaw/Dropbox/Projects/exposure_trajectories/",
              "manuscript/figures/mcar10_jmvn_obs_pred.jpeg"), 
       device = "jpeg", width = 7, height = 4.5, units = "in", dpi = 300)



