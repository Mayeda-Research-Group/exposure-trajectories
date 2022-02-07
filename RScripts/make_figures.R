#---- package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "data.table", "stringr", "openxlsx")

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

#---- load scripts ----
#put read_results function here

#---- color palette ----
# The palette with grey:
cbPalette <- c(
  # "#000000", # Remove the black color
  "#E69F00", "#56B4E9", "#009E73", "#FFD700", "#0072B2", 
  "#D55E00", "#CC79A7")

#---- count scenarios ----
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
                   "Percent", "Mechanism", "Truth Capture", "Time", "Seed"))
}

main_results <- do.call(rbind, lapply(main_paths, read_results)) %>% 
  #making sure only one copy of each seed
  na.omit() %>% group_by(Method, Exposure, Seed) %>% slice_head(n = 1) %>% 
  group_by(Method, Mechanism, Percent, Exposure)

sens_analyses <- do.call(rbind, lapply(sens_paths, read_results)) %>% 
  #making sure only one copy of each seed
  na.omit() %>% group_by(Method, Exposure, Seed) %>% slice_head(n = 1) %>% 
  group_by(Method, Mechanism, Percent, Exposure)

#---- **check scenario counts ----
#should be 1000 in each cell (divide by 4 for number of exposures)
table(main_results$Mechanism, main_results$Percent, main_results$Method)/4
table(sens_analyses$Mechanism, sens_analyses$Percent, sens_analyses$Method)/4

#---- **check max seeds for paused jobs ----
main_results %>% group_by(Method) %>% 
  summarise_at(.vars = c("Seed"), .funs = max)

sens_analyses %>% group_by(Method) %>% 
  summarise_at(.vars = c("Seed"), .funs = max)

#---- **check seeds overall ----
seeds <- seq(1, 9000, by = 1)

CC_main_missing_seeds <- main_results %>% 
  filter(Method == "CC") %>% ungroup() %>% dplyr::select("Seed") %>% 
  unique() %>% unlist() %>% setdiff(seeds, .) %>% as.data.frame() %>% 
  set_colnames("Seed") %>% mutate("Diff" = Seed - lag(Seed)) %>% 
  write.xlsx(paste0(path_to_dropbox,
                   "/exposure_trajectories/data/hoffman_transfer/missing_seeds/", 
                   "CC_main_missing.xlsx"))

PMM_main_missing_seeds <- main_results %>% 
  filter(Method == "PMM") %>% ungroup() %>% dplyr::select("Seed") %>% 
  unique() %>% unlist() %>% setdiff(seeds, .) %>% as.data.frame() %>% 
  set_colnames("Seed") %>% mutate("Diff" = Seed - lag(Seed)) %>% 
  write.xlsx(paste0(path_to_dropbox,
                   "/exposure_trajectories/data/hoffman_transfer/missing_seeds/", 
                   "PMM_main_missing.xlsx"))

FCS_main_missing_seeds <- main_results %>% 
  filter(Method == "FCS") %>% ungroup() %>% dplyr::select("Seed") %>% 
  unique() %>% unlist() %>% setdiff(seeds, .) %>% as.data.frame() %>% 
  set_colnames("Seed") %>% mutate("Diff" = Seed - lag(Seed)) %>% 
  write.xlsx(paste0(path_to_dropbox,
                   "/exposure_trajectories/data/hoffman_transfer/missing_seeds/", 
                   "FCS_main_missing.xlsx"))

CC_sens_missing_seeds <- sens_analyses %>% 
  filter(Method == "CC") %>% ungroup() %>% dplyr::select("Seed") %>% 
  unique() %>% unlist() %>% setdiff(seeds, .) %>% as.data.frame() %>% 
  set_colnames("Seed") %>% mutate("Diff" = Seed - lag(Seed)) %>% 
  write.xlsx(paste0(path_to_dropbox,
                   "/exposure_trajectories/data/hoffman_transfer/missing_seeds/", 
                   "CC_sens_missing.xlsx"))

JMVN_sens_missing_seeds <- sens_analyses %>% 
  filter(Method == "JMVN") %>% ungroup() %>% dplyr::select("Seed") %>% 
  unique() %>% unlist() %>% setdiff(seeds, .) %>% as.data.frame() %>% 
  set_colnames("Seed") %>% mutate("Diff" = Seed - lag(Seed)) %>% 
  write.xlsx(paste0(path_to_dropbox,
                   "/exposure_trajectories/data/hoffman_transfer/missing_seeds/", 
                   "JMVN_sens_missing.xlsx"))

PMM_sens_missing_seeds <- sens_analyses %>% 
  filter(Method == "PMM") %>% ungroup() %>% dplyr::select("Seed") %>% 
  unique() %>% unlist() %>% setdiff(seeds, .) %>% as.data.frame() %>% 
  set_colnames("Seed") %>% mutate("Diff" = Seed - lag(Seed)) %>% 
  write.xlsx(paste0(path_to_dropbox,
                    "/exposure_trajectories/data/hoffman_transfer/missing_seeds/", 
                    "PMM_sens_missing.xlsx"))

#---- **check seeds for each scenario ----
mechanisms <- c("mcar", "mar", "mnar")
percents <- c(10, 20, 30)

seed_vecs <- expand_grid(mechanisms, percents) %>% mutate("seed" = "seeds") %>% 
  unite("names", everything(), sep = "_") %>% unlist()

for(i in 1:length(seed_vecs)){
  assign(seed_vecs[i], seq(i, 9000, by = 9))
}

#---- **check for missing seeds ----
FCS_MCAR_10_missing_seeds <- main_results %>% 
  filter(Method == "FCS" & Mechanism == "MCAR" & Percent == "10%") %>% 
  ungroup() %>% dplyr::select("Seed") %>% unique() %>% unlist() %>% 
  setdiff(mcar_10_seeds, .) %>% as.data.frame() %>% 
  set_colnames("Seed") %>% mutate("Diff" = Seed - lag(Seed))

#---- Figure 2: results ----
#---- **get filepaths ----
all_paths <- 
  list.files(path = paste0(path_to_dropbox,
                           "/exposure_trajectories/data/hoffman_transfer/",
                           "results"), full.names = TRUE, pattern = "*.csv")

main_paths <- all_paths[!str_detect(all_paths, "sens")]

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

# #test for invalid rows
# colSums(is.na(main_results))
# colSums(is.na(sens_analyses))

#---- **limit runs for figure (for now) ----
main_results %<>% 
  group_by(Method, Mechanism, Percent, Exposure) %>% slice_head(n = 700) %>% 
  na.omit()

#double-checking
table(main_results$Mechanism, main_results$Percent, main_results$Method)/4

#---- **summarize data ----
results_summary <- main_results %>% 
  group_by(Method, Mechanism, Percent, Exposure) %>%
  summarize_at(.vars = c("Beta", "SE", "LCI", "UCI", "Truth Capture"), 
               ~mean(., na.rm = TRUE)) 

# #Sanity check
# results_sum_test <-
#   main_results %>% group_by(Method, Mechanism, Percent, Exposure) %>%
#   summarise_all(list(mean))
# 
# diffdf::diffdf(results_summary, results_sum_test)

#---- **read in truth table ----
truth <- read_csv(paste0(path_to_dropbox, 
                         "/exposure_trajectories/data/", "truth.csv")) %>%
  dplyr::rename("LCI" = "LCI_beta", 
                "UCI" = "UCI_beta",
                "Beta" = "beta") %>%
  # mutate("Percent" = "0%", 
  #        "Truth Capture" = 1) 
  mutate(Exposure = 
           case_when(
             Exposure == "CES-D Wave 4" ~ "Elevated Baseline CES-D",
             Exposure == "CES-D Wave 9" ~ "Elevated End of Follow-up CES-D",
             Exposure == "Elevated CES-D Prop" ~ "Proportion Elevated CES-D",
             TRUE ~ Exposure))

truth$Exposure <- 
  factor(truth$Exposure, 
         levels = c("Elevated Baseline CES-D", "Elevated End of Follow-up CES-D", 
                    "Elevated Average CES-D", "Proportion Elevated CES-D")) 

# truth_multiple <- do.call("rbind", replicate(
#   3, truth, simplify = FALSE)) %>%
#   mutate(Mechanism = c(rep("MCAR", length(unique(truth$Exposure))), 
#                        rep("MAR", length(unique(truth$Exposure))), 
#                        rep("MNAR", length(unique(truth$Exposure)))))

# results_summary %<>% rbind(truth_multiple)

#---- **format data ----
methods <- c("CC", "JMVN", "PMM", "FCS", "LMM")
results_summary$Method <- 
  factor(results_summary$Method, levels = c("Truth", methods))
results_summary$Percent <- factor(results_summary$Percent)
results_summary$Mechanism <- 
  factor(results_summary$Mechanism, levels = c("MCAR", "MAR", "MNAR")) 

results_summary[which(results_summary$Exposure == "CES-D Wave 4"), 
                "Exposure"] <- "Elevated Baseline CES-D"
results_summary[which(results_summary$Exposure == "CES-D Wave 9"), 
                "Exposure"] <- "Elevated End of Follow-up CES-D"
results_summary[which(results_summary$Exposure == "Elevated CES-D Prop"), 
                "Exposure"] <- "Proportion Elevated CES-D"
results_summary$Exposure <- 
  factor(results_summary$Exposure, 
         levels = c("Elevated Baseline CES-D", "Elevated End of Follow-up CES-D", 
                    "Elevated Average CES-D", "Proportion Elevated CES-D")) 

#---- **plot ----
ggplot(results_summary %>% filter(Method != "LMM"), 
       aes(x = Beta, y = Percent, color = Method, shape = Method)) +
  geom_point(size = 2.0, position = position_dodge(-0.75)) + 
  scale_shape_manual(values = c(rep("square", (nrow(results_summary))))) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = .3,
                position = position_dodge(-0.75)) +
  theme_minimal() + ylab("Percent Missing Data") +
  theme(legend.position = "bottom", legend.direction = "horizontal") + 
  scale_color_manual(values = cbPalette) + 
  scale_y_discrete(limits = rev(levels(results_summary$Percent))) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "dark grey") + 
  facet_grid(rows = vars(Mechanism), cols = vars(Exposure)) + 
  geom_vline(data = truth, aes(xintercept = Beta)) +
  ggtitle(paste0("Mean 95% CI of beta across 700 runs"))

ggsave(paste0(path_to_dropbox, "/exposure_trajectories/",
              "manuscript/figures/figure2/effect_ests_mean_CI.jpeg"), 
       device = "jpeg", dpi = 300, width = 9, height = 7, units = "in")

#---- Figure 3: coverage probabilities ----
#---- **get filepaths ----
all_paths <- 
  list.files(path = paste0(path_to_dropbox,
                           "/exposure_trajectories/data/hoffman_transfer/",
                           "results"), full.names = TRUE, pattern = "*.csv")

main_paths <- all_paths[!str_detect(all_paths, "sens")]

#---- **read in data ----
read_results <- function(paths){
  data.table::fread(paths, fill = TRUE) %>% na.omit() %>%
    set_colnames(c("Exposure", "Beta", "SE", "LCI", "UCI", "Method",
                   "Percent", "Mechanism", "Truth Capture", "Time"))
}

main_results <- do.call(rbind, lapply(main_paths, read_results)) %>% na.omit()

#---- **limit runs for figure (for now) ----
main_results %<>% 
  group_by(Method, Mechanism, Percent, Exposure) %>% slice_head(n = 700) %>% 
  na.omit()

#double-checking
table(main_results$Mechanism, main_results$Percent, main_results$Method)/4

#---- **summarize data ----
results_summary <- main_results %>% 
  group_by(Method, Mechanism, Percent, Exposure) %>%
  summarize_at(.vars = c("Truth Capture"), ~mean(., na.rm = TRUE))

#---- **format data ----
# Somehow this way, we got a plot with the order: "Truth, CC, JMVN in the plot
# but not in the legend
methods <- c("CC", "JMVN", "PMM", "FCS")
results_summary$Method <- factor(results_summary$Method, 
                                 levels = c("Truth", methods))
results_summary$Mechanism <- 
  factor(results_summary$Mechanism, levels = c("MCAR", "MAR", "MNAR"))
results_summary$Percent <- factor(results_summary$Percent)

results_summary[which(results_summary$Exposure == "CES-D Wave 4"), 
                "Exposure"] <- "Elevated Baseline CES-D"
results_summary[which(results_summary$Exposure == "CES-D Wave 9"), 
                "Exposure"] <- "Elevated End of Follow-up CES-D"
results_summary[which(results_summary$Exposure == "Elevated CES-D Prop"), 
                "Exposure"] <- "Proportion Elevated CES-D"
results_summary$Exposure <- 
  factor(results_summary$Exposure, 
         levels = c("Elevated Baseline CES-D", "Elevated End of Follow-up CES-D", 
                    "Elevated Average CES-D", "Proportion Elevated CES-D")) 

#---- **plot ----
ggplot(results_summary %>% filter(!Method == "Truth"), 
       mapping = aes(x = Percent, y = `Truth Capture`, 
                     color = Method)) +
  geom_point(alpha = 0.75) + geom_line(aes(group = Method), alpha = 0.75) + 
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal") + 
  scale_color_manual(values = cbPalette[-1]) + ylab("Coverage Probability") + 
  facet_grid(rows = vars(Mechanism), cols = vars(Exposure), scales = "free_y")

ggsave(paste0(path_to_dropbox, "/exposure_trajectories/",
              "manuscript/figures/figure3/coverage_prob.jpeg"), 
       device = "jpeg", dpi = 300, width = 9, height = 7, units = "in")

#---- Figure 4: RMSE ----
#---- **read in data ----
rmse_table <- read_csv(paste0(path_to_dropbox, "/exposure_trajectories/",
                              "manuscript/tables/table2/rmse.csv"))
rmse_table %<>% 
  pivot_longer(cols = colnames(rmse_table)[grep("CES-D", 
                                                colnames(rmse_table))]) %>% 
  filter(Method != "LMM")

#---- **format data ----
methods <- c("CC", "JMVN", "PMM", "FCS")
rmse_table$Method <- factor(rmse_table$Method, levels = methods)
rmse_table$Mechanism <- factor(rmse_table$Mechanism, 
                               levels = c("MCAR", "MAR", "MNAR"))
rmse_table$`Missing Percent` <- factor(rmse_table$`Missing Percent`)

rmse_table[which(rmse_table$name == "CES-D Wave 4"), "name"] <- 
  "Elevated Baseline CES-D"
rmse_table[which(rmse_table$name == "CES-D Wave 9"), "name"] <- 
  "Elevated End of Follow-up CES-D"
rmse_table[which(rmse_table$name == "Elevated CES-D Prop"), "name"] <- 
  "Proportion Elevated CES-D"
rmse_table$name <- 
  factor(rmse_table$name, 
         levels = c("Elevated Baseline CES-D", "Elevated End of Follow-up CES-D", 
                    "Elevated Average CES-D", "Proportion Elevated CES-D"))

#---- **plot ----
ggplot(rmse_table, 
       mapping = aes(x = `Missing Percent`, y = value, 
                     color = Method)) +
  geom_point(alpha = 0.75) + geom_line(aes(group = Method), alpha = 0.75) + 
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal") + 
  scale_color_manual(values = cbPalette[-1]) + ylab("RMSE") + 
  facet_grid(rows = vars(Mechanism), cols = vars(name), scales = "free_y")

ggsave(paste0(path_to_dropbox, "/exposure_trajectories/",
              "manuscript/figures/figure4/rmse.jpeg"), 
       device = "jpeg", dpi = 300, width = 9, height = 7, units = "in")

#---- Figure 5: runtimes ----
#---- **get filepaths ----
all_paths <- 
  list.files(path = paste0(path_to_dropbox,
                           "/exposure_trajectories/data/hoffman_transfer/",
                           "results"), full.names = TRUE, pattern = "*.csv")

main_paths <- all_paths[!str_detect(all_paths, "sens")]

#---- **read in data ----
read_results <- function(paths){
  data.table::fread(paths, fill = TRUE) %>% na.omit() %>%
    set_colnames(c("Exposure", "Beta", "SE", "LCI", "UCI", "Method",
                   "Percent", "Mechanism", "Truth Capture", "Time"))
}

main_results <- do.call(rbind, lapply(main_paths, read_results)) %>% na.omit()

#---- **limit runs for figure (for now) ----
main_results %<>% 
  group_by(Method, Mechanism, Percent, Exposure) %>% slice_head(n = 700) %>% 
  na.omit()

#double-checking
table(main_results$Mechanism, main_results$Percent, main_results$Method)/4

#---- **summarize data ----
main_run_times <- main_results %>% 
  group_by(Method) %>% summarize_at(.vars = c("Time"), ~mean(., na.rm = TRUE)) 

#---- **format data ----
methods <- c("CC", "JMVN", "PMM", "FCS")
main_results$Method <- factor(main_results$Method, 
                              levels = c("Truth", methods))
main_results$Mechanism <- 
  factor(main_results$Mechanism, levels = c("MCAR", "MAR", "MNAR"))
main_results$Percent <- factor(main_results$Percent)

main_results %<>% mutate("time_hours" = Time/60)

#---- **plot ----
ggplot(data = na.omit(main_results), 
       aes(x = Percent, y = time_hours, color = Method)) + 
  geom_boxplot() + ylab("Computational Time (Hours)") + 
  xlab("Percent Missing Data") + theme_bw() + 
  theme(legend.position = "bottom", legend.direction = "horizontal") + 
  scale_color_manual(values = cbPalette[-1])

ggsave(paste0(path_to_dropbox, "/exposure_trajectories/",
              "manuscript/figures/figure5/run_times.jpeg"), 
       device = "jpeg", dpi = 300, width = 7, height = 5, units = "in")

#---- eFigure 3: sensitivity analyses ----
#---- **get filepaths ----
all_paths <- 
  list.files(path = paste0(path_to_dropbox,
                           "/exposure_trajectories/data/hoffman_transfer/",
                           "results"), full.names = TRUE, pattern = "*.csv")

sens_paths <- all_paths[str_detect(all_paths, "sens")]

#---- **read in data ----
read_results <- function(paths){
  data.table::fread(paths, fill = TRUE) %>% na.omit() %>%
    set_colnames(c("Exposure", "Beta", "SE", "LCI", "UCI", "Method",
                   "Percent", "Mechanism", "Truth Capture", "Time"))
}

sens_analyses <- do.call(rbind, lapply(sens_paths, read_results)) %>% na.omit()

#---- **limit runs for figure (for now) ----
sens_analyses %<>% 
  group_by(Method, Mechanism, Percent, Exposure) %>% slice_head(n = 100) %>% 
  na.omit()

#double-checking
table(sens_analyses$Mechanism, sens_analyses$Percent, sens_analyses$Method)/4

#---- **summarize data ----
results_summary <- sens_analyses %>% 
  group_by(Method, Mechanism, Percent, Exposure) %>%
  summarize_at(.vars = c("Beta", "SE", "LCI", "UCI", "Truth Capture"), 
               ~mean(., na.rm = TRUE)) 

#---- **read in truth table ----
truth_sens <- read_csv(paste0(path_to_dropbox, 
                              "/exposure_trajectories/data/", "truth_sens.csv")) %>%
  dplyr::rename("LCI" = "LCI_beta", 
                "UCI" = "UCI_beta",
                "Beta" = "beta") %>% 
  mutate("Percent" = "0%", 
         "Truth Capture" = 1)

truth_multiple <- do.call("rbind", replicate(
  3, truth_sens, simplify = FALSE)) %>%
  mutate(Mechanism = c(rep("MCAR", length(unique(truth_sens$Exposure))), 
                       rep("MAR", length(unique(truth_sens$Exposure))), 
                       rep("MNAR", length(unique(truth_sens$Exposure)))))

results_summary %<>% rbind(truth_multiple)

#---- **format data ----
methods <- c("CC", "JMVN", "PMM", "FCS")
results_summary$Method <- 
  factor(results_summary$Method, levels = c("Truth", methods))
results_summary$Percent <- factor(results_summary$Percent)
results_summary$Mechanism <- 
  factor(results_summary$Mechanism, levels = c("MCAR", "MAR", "MNAR")) 

results_summary[which(results_summary$Exposure == "CES-D Wave 4"), 
                "Exposure"] <- "Elevated Baseline CES-D"
results_summary[which(results_summary$Exposure == "CES-D Wave 9"), 
                "Exposure"] <- "Elevated End of Follow-up CES-D"
results_summary[which(results_summary$Exposure == "Elevated CES-D Prop"), 
                "Exposure"] <- "Proportion Elevated CES-D"
results_summary$Exposure <- 
  factor(results_summary$Exposure, 
         levels = c("Elevated Baseline CES-D", "Elevated End of Follow-up CES-D", 
                    "Elevated Average CES-D", "Proportion Elevated CES-D")) 

#---- **plot ----
ggplot(results_summary, 
       aes(x = Beta, y = Percent, color = Method, shape = Method)) +
  geom_point(size = 2.0, position = position_dodge(-0.75)) + 
  scale_shape_manual(values = c(rep("square", (nrow(results_summary))))) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = .3,
                position = position_dodge(-0.75)) +
  theme_minimal() + ylab("Percent Missing Data") +
  theme(legend.position = "bottom", legend.direction = "horizontal") + 
  scale_color_manual(values = cbPalette) + 
  scale_y_discrete(limits = rev(levels(results_summary$Percent))) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
  facet_grid(rows = vars(Mechanism), cols = vars(Exposure)) + 
  ggtitle(paste0("Mean 95% CI of beta across 100 runs"))

ggsave(paste0(path_to_dropbox, "/exposure_trajectories/",
              "manuscript/figures/efigure3/effect_ests_mean_CI_sens.jpeg"), 
       device = "jpeg", dpi = 300, width = 9, height = 7, units = "in")

#---- eFigure 4: sensitivity coverage probabilities ----
#---- **get filepaths ----
all_paths <- 
  list.files(path = paste0(path_to_dropbox,
                           "/exposure_trajectories/data/hoffman_transfer/",
                           "results"), full.names = TRUE, pattern = "*.csv")

sens_paths <- all_paths[str_detect(all_paths, "sens")]

#---- **read in data ----
read_results <- function(paths){
  data.table::fread(paths, fill = TRUE) %>% na.omit() %>%
    set_colnames(c("Exposure", "Beta", "SE", "LCI", "UCI", "Method",
                   "Percent", "Mechanism", "Truth Capture", "Time"))
}

sens_analyses <- do.call(rbind, lapply(sens_paths, read_results)) %>% na.omit()

#---- **limit runs for figure (for now) ----
sens_analyses %<>% 
  group_by(Method, Mechanism, Percent, Exposure) %>% slice_head(n = 100) %>% 
  na.omit()

#double-checking
table(sens_analyses$Mechanism, sens_analyses$Percent, sens_analyses$Method)/4

#---- **summarize data ----
results_summary <- sens_analyses %>% 
  group_by(Method, Mechanism, Percent, Exposure) %>%
  summarize_at(.vars = c("Truth Capture"), ~mean(., na.rm = TRUE))

#---- **format data ----
# Somehow this way, we got a plot with the order: "Truth, CC, JMVN in the plot
# but not in the legend
methods <- c("CC", "JMVN", "PMM", "FCS")
results_summary$Method <- factor(results_summary$Method, 
                                 levels = c("Truth", methods))
results_summary$Mechanism <- 
  factor(results_summary$Mechanism, levels = c("MCAR", "MAR", "MNAR"))
results_summary$Percent <- factor(results_summary$Percent)

results_summary[which(results_summary$Exposure == "CES-D Wave 4"), 
                "Exposure"] <- "Elevated Baseline CES-D"
results_summary[which(results_summary$Exposure == "CES-D Wave 9"), 
                "Exposure"] <- "Elevated End of Follow-up CES-D"
results_summary[which(results_summary$Exposure == "Elevated CES-D Prop"), 
                "Exposure"] <- "Proportion Elevated CES-D"
results_summary$Exposure <- 
  factor(results_summary$Exposure, 
         levels = c("Elevated Baseline CES-D", "Elevated End of Follow-up CES-D", 
                    "Elevated Average CES-D", "Proportion Elevated CES-D")) 

#---- **plot ----
ggplot(results_summary %>% filter(!Method == "Truth"), 
       mapping = aes(x = Percent, y = `Truth Capture`, 
                     color = Method)) +
  geom_point(alpha = 0.75) + geom_line(aes(group = Method), alpha = 0.75) + 
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal") + 
  scale_color_manual(values = cbPalette[-1]) + ylab("Coverage Probability") + 
  facet_grid(rows = vars(Mechanism), cols = vars(Exposure), scales = "free_y")

ggsave(paste0(path_to_dropbox, "/exposure_trajectories/",
              "manuscript/figures/efigure4/coverage_prob_sens.jpeg"), 
       device = "jpeg", dpi = 300, width = 9, height = 7, units = "in")

#---- eFigure 5: sensitivity RMSE ----
#---- **read in data ----
rmse_table <- read_csv(paste0(path_to_dropbox, "/exposure_trajectories/",
                              "manuscript/tables/etable2/rmse_sens.csv"))
rmse_table %<>% 
  pivot_longer(cols = colnames(rmse_table)[grep("CES-D", 
                                                colnames(rmse_table))]) %>% 
  filter(Method != "LMM")

#---- **format data ----
methods <- c("CC", "JMVN", "PMM", "FCS")
rmse_table$Method <- factor(rmse_table$Method, levels = methods)
rmse_table$Mechanism <- factor(rmse_table$Mechanism, 
                               levels = c("MCAR", "MAR", "MNAR"))
rmse_table$`Missing Percent` <- factor(rmse_table$`Missing Percent`)

rmse_table[which(rmse_table$name == "CES-D Wave 4"), "name"] <- 
  "Elevated Baseline CES-D"
rmse_table[which(rmse_table$name == "CES-D Wave 9"), "name"] <- 
  "Elevated End of Follow-up CES-D"
rmse_table[which(rmse_table$name == "Elevated CES-D Prop"), "name"] <- 
  "Proportion Elevated CES-D"
rmse_table$name <- 
  factor(rmse_table$name, 
         levels = c("Elevated Baseline CES-D", "Elevated End of Follow-up CES-D", 
                    "Elevated Average CES-D", "Proportion Elevated CES-D"))

#---- **plot ----
ggplot(rmse_table, 
       mapping = aes(x = `Missing Percent`, y = value, 
                     color = Method)) +
  geom_point(alpha = 0.75) + geom_line(aes(group = Method), alpha = 0.75) + 
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal") + 
  scale_color_manual(values = cbPalette[-1]) + ylab("RMSE") + 
  facet_grid(rows = vars(Mechanism), cols = vars(name), scales = "free_y")

ggsave(paste0(path_to_dropbox, "/exposure_trajectories/",
              "manuscript/figures/efigure5/rmse_sens.jpeg"), 
       device = "jpeg", dpi = 300, width = 9, height = 7, units = "in")

#---- eFigure 6: runtimes ----
#---- **get filepaths ----
all_paths <- 
  list.files(path = paste0(path_to_dropbox,
                           "/exposure_trajectories/data/hoffman_transfer/",
                           "results"), full.names = TRUE, pattern = "*.csv")

sens_paths <- all_paths[str_detect(all_paths, "sens")]

#---- **read in data ----
read_results <- function(paths){
  data.table::fread(paths, fill = TRUE) %>% na.omit() %>%
    set_colnames(c("Exposure", "Beta", "SE", "LCI", "UCI", "Method",
                   "Percent", "Mechanism", "Truth Capture", "Time"))
}

sens_analyses <- do.call(rbind, lapply(sens_paths, read_results)) %>% na.omit()

#---- **limit runs for figure (for now) ----
sens_analyses %<>% 
  group_by(Method, Mechanism, Percent, Exposure) %>% slice_head(n = 100) %>% 
  na.omit()

#double-checking
table(sens_analyses$Mechanism, sens_analyses$Percent, sens_analyses$Method)/4

#---- **summarize data ----
sens_run_times <- sens_analyses %>% 
  group_by(Method) %>% summarize_at(.vars = c("Time"), ~mean(., na.rm = TRUE)) 

#---- **format data ----
methods <- c("CC", "JMVN", "PMM", "FCS")
sens_analyses$Method <- factor(sens_analyses$Method, 
                               levels = c("Truth", methods))
sens_analyses$Mechanism <- 
  factor(sens_analyses$Mechanism, levels = c("MCAR", "MAR", "MNAR"))
sens_analyses$Percent <- factor(sens_analyses$Percent)

sens_analyses %<>% mutate("time_hours" = Time/60)

#---- **plot ----
ggplot(data = na.omit(sens_analyses), 
       aes(x = Percent, y = time_hours, color = Method)) + 
  geom_boxplot() + ylab("Computational Time (Hours)") + 
  xlab("Percent Missing Data") + theme_bw() + 
  theme(legend.position = "bottom", legend.direction = "horizontal") + 
  scale_color_manual(values = cbPalette[-1])

ggsave(paste0(path_to_dropbox, "/exposure_trajectories/",
              "manuscript/figures/efigure6/run_times_sens.jpeg"), 
       device = "jpeg", dpi = 300, width = 7, height = 5, units = "in")

#---- OLD ----
#---- eFigure 2: traceplots ----
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



