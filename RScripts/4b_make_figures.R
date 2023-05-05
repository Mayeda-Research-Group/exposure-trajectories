#This script creates all the manuscript figures

#---- package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "data.table", "stringr", "openxlsx", 
       "lemon", "cowplot")

#No scientific notation
options(scipen = 999)

#---- note ----
# Since the difference between win and OS, put substituted directory here
# Yingyan's directory: "C:/Users/Yingyan Wu/Box"
#                      
# Crystal's directory: /Users/crystalshaw/Library/CloudStorage/Box-Box/Projects
#                     

#Changing directories here will change them throughout the script
path_to_box <- "C:/Users/Yingyan Wu/Box"

#---- load scripts ----
source(here::here("RScripts", "functions", "read_results.R"))

#---- color palette ----
# The palette with grey:
cbPalette <- c(
  # "#000000", # Remove the black color
  "#E69F00", "#56B4E9", "#009E73", "#FFD700", "#0072B2", 
  "#D55E00", "#CC79A7")

#---- count scenarios ----
#---- **get filepaths ----
all_paths <- 
  list.files(path = paste0(path_to_box,
                           "/exposure_trajectories/data/hoffman_transfer/",
                           "results"), full.names = TRUE, pattern = "*.csv")

main_paths <- all_paths[!str_detect(all_paths, "sens")]
sens_paths <- all_paths[str_detect(all_paths, "sens")]

#---- **read in data ----
main_results <- do.call(rbind, lapply(main_paths, read_results)) %>% 
  #making sure only one copy of each seed
  na.omit() %>% group_by(Method, Mechanism, Percent, Exposure, Seed) %>% 
  slice_head(n = 1) %>% group_by(Method, Mechanism, Percent, Exposure)

sens_analyses <- do.call(rbind, lapply(sens_paths, read_results)) %>% 
  #making sure only one copy of each seed
  na.omit() %>% group_by(Method, Mechanism, Percent, Exposure, Seed) %>% 
  slice_head(n = 1) %>% group_by(Method, Mechanism, Percent, Exposure)

#---- **check scenario counts ----
#should be 1000 in each cell (divide by 4 for number of exposures)
table(main_results$Mechanism, main_results$Percent, main_results$Method)/4
table(sens_analyses$Mechanism, sens_analyses$Percent, sens_analyses$Method)/4

# #---- **check max seeds for paused jobs ----
# main_results %>% group_by(Method) %>% 
#   summarise_at(.vars = c("Seed"), .funs = max)
# 
# sens_analyses %>% group_by(Method) %>% 
#   summarise_at(.vars = c("Seed"), .funs = max)

#---- **check seeds overall ----
seeds <- seq(1, 9000, by = 1)

#---- ****CC main ----
main_results %>% filter(Method == "CC") %>% ungroup() %>% 
  dplyr::select("Seed") %>% unique() %>% unlist() %>% setdiff(seeds, .) %>% 
  as.data.frame() %>% set_colnames("Seed") %>% 
  mutate("Diff" = Seed - lag(Seed)) %>% 
  write.xlsx(paste0(path_to_box,
                    "/exposure_trajectories/data/hoffman_transfer/missing_seeds/", 
                    "CC_main_missing.xlsx"), overwrite = TRUE)

#---- ****JMVN main ----
main_results %>% filter(Method == "JMVN") %>% ungroup() %>% 
  dplyr::select("Seed") %>% unique() %>% unlist() %>% setdiff(seeds, .) %>% 
  as.data.frame() %>% set_colnames("Seed") %>% 
  mutate("Diff" = Seed - lag(Seed)) %>% 
  write.xlsx(paste0(path_to_box,
                    "/exposure_trajectories/data/hoffman_transfer/missing_seeds/", 
                    "JMVN_main_missing.xlsx"), overwrite = TRUE)

#---- ****PMM main ----
main_results %>% filter(Method == "PMM") %>% ungroup() %>% 
  dplyr::select("Seed") %>% unique() %>% unlist() %>% setdiff(seeds, .) %>% 
  as.data.frame() %>% set_colnames("Seed") %>% 
  mutate("Diff" = Seed - lag(Seed)) %>% 
  write.xlsx(paste0(path_to_box,
                    "/exposure_trajectories/data/hoffman_transfer/missing_seeds/", 
                    "PMM_main_missing.xlsx"), overwrite = TRUE)

#---- ****FCS main ----
main_results %>% filter(Method == "FCS") %>% ungroup() %>% 
  dplyr::select("Seed") %>% unique() %>% unlist() %>% setdiff(seeds, .) %>% 
  as.data.frame() %>% set_colnames("Seed") %>% 
  mutate("Diff" = Seed - lag(Seed)) %>% 
  write.xlsx(paste0(path_to_box,
                    "/exposure_trajectories/data/hoffman_transfer/missing_seeds/", 
                    "FCS_main_missing.xlsx"), overwrite = TRUE)

#---- ****LMM main (MNAR 30%) ----
LMM_seeds <- seq(1, 1000, by = 1)
main_results %>% filter(Method == "LMM" & Mechanism == "MNAR") %>% ungroup() %>% 
  dplyr::select("Seed") %>% unique() %>% unlist() %>% setdiff(LMM_seeds, .) %>% 
  as.data.frame() %>% set_colnames("Seed") %>% 
  mutate("Diff" = Seed - lag(Seed)) %>% 
  write.xlsx(paste0(path_to_box,
                    "/exposure_trajectories/data/hoffman_transfer/missing_seeds/", 
                    "LMM_main_MNAR_30_missing.xlsx"), overwrite = TRUE)

#---- ****CC sens ----
sens_analyses %>% filter(Method == "CC") %>% ungroup() %>% 
  dplyr::select("Seed") %>% unique() %>% unlist() %>% setdiff(seeds, .) %>% 
  as.data.frame() %>% set_colnames("Seed") %>% 
  mutate("Diff" = Seed - lag(Seed)) %>% 
  write.xlsx(paste0(path_to_box,
                    "/exposure_trajectories/data/hoffman_transfer/missing_seeds/", 
                    "CC_sens_missing.xlsx"), overwrite = TRUE)

#---- ****JVMN sens ----
sens_analyses %>% filter(Method == "JMVN") %>% ungroup() %>% 
  dplyr::select("Seed") %>% unique() %>% unlist() %>% setdiff(seeds, .) %>% 
  as.data.frame() %>% set_colnames("Seed") %>% 
  mutate("Diff" = Seed - lag(Seed)) %>% 
  write.xlsx(paste0(path_to_box,
                    "/exposure_trajectories/data/hoffman_transfer/missing_seeds/", 
                    "JMVN_sens_missing.xlsx"), overwrite = TRUE)

#---- ****PMM sens ----
sens_analyses %>% filter(Method == "PMM") %>% ungroup() %>% 
  dplyr::select("Seed") %>% unique() %>% unlist() %>% setdiff(seeds, .) %>% 
  as.data.frame() %>% set_colnames("Seed") %>% 
  mutate("Diff" = Seed - lag(Seed)) %>% 
  write.xlsx(paste0(path_to_box,
                    "/exposure_trajectories/data/hoffman_transfer/missing_seeds/", 
                    "PMM_sens_missing.xlsx"), overwrite = TRUE)

#---- ****FCS sens ----
sens_analyses %>% filter(Method == "FCS") %>% ungroup() %>% 
  dplyr::select("Seed") %>% unique() %>% unlist() %>% setdiff(seeds, .) %>% 
  as.data.frame() %>% set_colnames("Seed") %>% 
  mutate("Diff" = Seed - lag(Seed)) %>% 
  write.xlsx(paste0(path_to_box,
                    "/exposure_trajectories/data/hoffman_transfer/missing_seeds/", 
                    "FCS_sens_missing.xlsx"), overwrite = TRUE)

#---- **check seeds by scenario ----
MCAR_10_seeds <- seq(1, 9000, by = 9)
MCAR_20_seeds <- seq(2, 9000, by = 9)
MCAR_30_seeds <- seq(3, 9000, by = 9)

MAR_10_seeds <- seq(4, 9000, by = 9)
MAR_20_seeds <- seq(5, 9000, by = 9)
MAR_30_seeds <- seq(6, 9000, by = 9)

MNAR_10_seeds <- seq(7, 9000, by = 9)
MNAR_20_seeds <- seq(8, 9000, by = 9)
MNAR_30_seeds <- seq(9, 9000, by = 9)

#---- ****FCS main ----
for(mechanism in c("MCAR", "MAR", "MNAR")){
  for(percent in c(10, 20, 30)){
    main_results %>% 
      filter(Method == "FCS" & Mechanism == mechanism & 
               Percent == paste0(percent, "%")) %>% ungroup() %>% 
      dplyr::select("Seed") %>% unique() %>% unlist() %>% 
      setdiff(get(paste0(mechanism, "_", percent, "_seeds")), .) %>% 
      as.data.frame() %>% set_colnames("Seed") %>% 
      mutate("Diff" = Seed - lag(Seed)) %>% 
      write.xlsx(paste0(path_to_box,
                        "/exposure_trajectories/data/hoffman_transfer/", 
                        "missing_seeds/FCS_main_missing_", mechanism, "_", 
                        percent, ".xlsx"), 
                 overwrite = TRUE)
  }
}

#---- ****FCS sens ----
for(mechanism in c("MCAR", "MAR", "MNAR")){
  for(percent in c(10, 20, 30)){
    sens_analyses %>% 
      filter(Method == "FCS" & Mechanism == mechanism & 
               Percent == paste0(percent, "%")) %>% ungroup() %>% 
      dplyr::select("Seed") %>% unique() %>% unlist() %>% 
      setdiff(get(paste0(mechanism, "_", percent, "_seeds")), .) %>% 
      as.data.frame() %>% set_colnames("Seed") %>% 
      mutate("Diff" = Seed - lag(Seed)) %>% 
      write.xlsx(paste0(path_to_box,
                        "/exposure_trajectories/data/hoffman_transfer/", 
                        "missing_seeds/FCS_sens_missing_", mechanism, "_", 
                        percent, ".xlsx"), 
                 overwrite = TRUE)
  }
}

#---- Figure 1 + eFigure 7: results ----
#---- **get filepaths ----
all_paths <- 
  list.files(path = paste0(path_to_box,
                           "/exposure_trajectories/data/hoffman_transfer/",
                           "results"), full.names = TRUE, pattern = "*.csv")

main_paths <- all_paths[!str_detect(all_paths, "sens")]

#---- **read in data ----
main_results <- do.call(rbind, lapply(main_paths, read_results)) %>% 
  #making sure only one copy of each seed
  na.omit() %>% group_by(Method, Mechanism, Percent, Exposure, Seed) %>% 
  slice_head(n = 1) %>% 
  group_by(Method, Mechanism, Percent, Exposure)

#rename FCS --> VTS
main_results %<>% mutate("Method" = ifelse(Method == "FCS", "VTS", Method))
#rename JMVN --> NORM
main_results %<>% mutate("Method" = ifelse(Method == "JMVN", "NORM", Method))

# #test for invalid rows
# colSums(is.na(main_results))
# colSums(is.na(sens_analyses))

#---- **double-checking ----
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
truth <- read_csv(paste0(path_to_box, 
                         "/exposure_trajectories/data/", "truth.csv")) %>%
  dplyr::rename("LCI" = "LCI_beta", 
                "UCI" = "UCI_beta",
                "Beta" = "beta") %>%
  # mutate("Percent" = "0%", 
  #        "Truth Capture" = 1) 
  mutate(Exposure = 
           case_when(
             Exposure == "CES-D Wave 4" ~ "Elevated Baseline CES-D,\n1998",
             Exposure == "CES-D Wave 9" ~ "Elevated End of Exposure CES-D, \n2008",
             Exposure == 
               "Elevated CES-D Prop" ~ "Proportion Elevated CES-D,\n1998-2008",
             TRUE ~ "Elevated Average CES-D,\n1998-2008"))

truth$Exposure <- 
  factor(truth$Exposure, 
         levels = c("Elevated Baseline CES-D,\n1998", 
                    "Elevated End of Exposure CES-D, \n2008", 
                    "Elevated Average CES-D,\n1998-2008", 
                    "Proportion Elevated CES-D,\n1998-2008")) 

#---- **format data ----
methods <- c("CC", "NORM", "PMM", "VTS", "LMM")
results_summary$Method <- 
  factor(results_summary$Method, levels = c("Truth", methods))
results_summary$Percent <- factor(results_summary$Percent)
results_summary$Mechanism <- 
  factor(results_summary$Mechanism, levels = c("MCAR", "MAR", "MNAR")) 

results_summary[which(results_summary$Exposure == "CES-D Wave 4"), 
                "Exposure"] <- "Elevated Baseline CES-D,\n1998"
results_summary[which(results_summary$Exposure == "CES-D Wave 9"), 
                "Exposure"] <- "Elevated End of Exposure CES-D, \n2008"
results_summary[which(results_summary$Exposure == "Elevated Average CES-D"), 
                "Exposure"] <- "Elevated Average CES-D,\n1998-2008"
results_summary[which(results_summary$Exposure == "Elevated CES-D Prop"), 
                "Exposure"] <- "Proportion Elevated CES-D,\n1998-2008"
results_summary$Exposure <- 
  factor(results_summary$Exposure, 
         levels = c("Elevated Baseline CES-D,\n1998", 
                    "Elevated End of Exposure CES-D, \n2008", 
                    "Elevated Average CES-D,\n1998-2008", 
                    "Proportion Elevated CES-D,\n1998-2008")) 

#---- **figure 1 plot ----
#cowplot
plot_vars <- 
  expand.grid(unique(results_summary$Exposure), 
              unique(results_summary$Mechanism)) %>% 
  set_colnames(c("Exposure", "Mechanism")) %>% arrange(Mechanism) %>% 
  mutate("Label" = paste0(LETTERS[1:12], ")"))

plot_data <- results_summary %>% filter(!Method == "LMM")
plot_data$Percent <- str_remove(plot_data$Percent, "%")
plot_data$Percent <- factor(plot_data$Percent)

figure1_plot_list <- list()

for(row in 1:nrow(plot_vars)){
  mech = plot_vars[row, "Mechanism"]
  exp = plot_vars[row, "Exposure"]
  label = plot_vars[row, "Label"]
  
  figure1_plot_list[[row]] <- 
    ggplot(plot_data %>% filter(Mechanism == mech & Exposure == exp), 
           aes(x = Beta, y = Percent, color = Method, shape = Method)) +
    geom_point(position = position_dodge(-0.8)) +
    scale_shape_manual(values = 
                         c(rep("square", (nrow(results_summary))))) + 
    geom_errorbar(aes(xmin = LCI, xmax = UCI), 
                  position = position_dodge(-0.8)) +
    theme_bw() + ylab("Missing Data, %") +
    theme(legend.position = "none") +
    scale_color_manual(values = cbPalette) + 
    scale_y_discrete(limits = rev(levels(plot_data$Percent))) + 
    scale_x_continuous(limits = c(min(results_summary$LCI), 
                                  max(results_summary$UCI))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "dark grey") + 
    geom_vline(data = truth %>% filter(Exposure == exp), 
               linetype = "dashed", aes(xintercept = Beta)) + 
    xlab(expression(beta)) + 
    theme(text = element_text(size = 8, color = "black"), 
          axis.text.x = element_text(size = 8, color = "black"), 
          axis.text.y = element_text(size = 8, color = "black")) + 
    theme(panel.border = element_blank(), axis.line = element_line(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(t = 13, r = 1, b = 1, l = 13), unit = "pt"))
  
  if (row %in% c(4, 8)){
    # for plot D and H, label the far right X axis.
    figure1_plot_list[[row]] <- figure1_plot_list[[row]] +
      scale_x_continuous(limits = c(min(results_summary$LCI), 
                                    max(results_summary$UCI)),
                         breaks = c(-1.0, -0.5, 0.0, 0.5, 0.7))
  }
  
  if (row %in% c(11, 12)){
    # for plot K and L, label the far left X axis.
    figure1_plot_list[[row]] <- figure1_plot_list[[row]] +
      scale_x_continuous(limits = c(min(results_summary$LCI), 
                                    max(results_summary$UCI)),
                         breaks = c(-1.4, -1.0, -0.5, 0.0, 0.5))
  }
  
  figure1_plot_list[[row]] <- 
    plot_grid(figure1_plot_list[[row]], labels = plot_vars$Label[[row]], 
              align = "vh", label_size = 8, hjust = 0, vjust = 1, 
              label_fontface = "plain")
  
  ggsave(plot = figure1_plot_list[[row]],
         filename = paste0(path_to_box, "/exposure_trajectories/",
                           "manuscript/figures/figure1/figure1_panel_", 
                           substr(plot_vars$Label, 1, 1)[[row]], ".pdf"),
         device = "pdf", dpi = 300, width = 7/4, height = 5/3, units = "in")
  
  ggsave(plot = figure1_plot_list[[row]],
         filename = paste0(path_to_box, "/exposure_trajectories/",
                           "manuscript/figures/figure1/figure1_panel_", 
                           substr(plot_vars$Label, 1, 1)[[row]], ".eps"),
         device = "eps", dpi = 300, width = 7/4, height = 5/3, units = "in")
  
}

figure1_plot_forlegend <- 
  ggplot(plot_data, 
         aes(x = Beta, y = Percent, color = Method, shape = Method)) +
  geom_point(size = 2.25, position = position_dodge(-0.8)) + 
  scale_shape_manual(values = 
                       c(rep("square", (nrow(results_summary))))) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = .5, size = 0.75,
                position = position_dodge(-0.8)) +
  theme_bw() + ylab("Missing Data, %") +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  scale_color_manual(values = cbPalette) + 
  theme(text = element_text(size = 8, color = "black"), 
        axis.text.x = element_text(size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"),
        legend.background = 
          element_rect(fill = "white", linetype = "solid", colour ="black"), 
        legend.title.align = 0.5,
        legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8)) + 
  guides(shape = guide_legend(title = expression(underline(Method)), 
                              title.position = "top",
                              title.vjust = -1.3,
                              title.theme = element_text(size = 8, color = "black")),
         color = guide_legend(title = expression(underline(Method)), 
                              title.position = "top",
                              title.vjust = -1.3,
                              title.theme = element_text(size = 8, color = "black")))

legend_b <- ggpubr::get_legend(figure1_plot_forlegend)
# plot(legend_b)

figure1_panel <- 
  plot_grid(plotlist = figure1_plot_list, align = "vh", 
            ncol = 4) +
  theme(plot.margin = unit(c(t = 0, r = 0, b = 8, l = 0), unit = "pt"))

figure1_panel_final <- plot_grid(figure1_panel, 
                                 legend_b,
                                 ncol = 1, rel_heights = c(1, .1)) +
  theme(plot.margin = unit(c(t = 21, r = 10, b = 10, l = 0), unit = "pt")) +
  geom_text(data = data.frame(
    x = seq(0.13, 0.91, by = 0.26), y = rep(1, 4),
    label = paste0(levels(plot_vars$Exposure), "\n\n")),
    # This is stupid and I didn't figure out how to solve this without putting two \ns
    mapping = aes(x = x, y = y, label = label),
    size = 5/14*8, inherit.aes = FALSE) + 
  geom_text(data = data.frame(x = rep(1, 3), y = c(0.9, 0.6, 0.3),
                              label = paste0(levels(plot_vars$Mechanism), "\n")),
            mapping = aes(x = x, y = y, label = label),
            size = 5/14*8, angle = -90L, inherit.aes = FALSE)
figure1_panel_final

ggsave(plot = figure1_panel_final,
       filename = paste0(path_to_box, "/exposure_trajectories/",
                         "manuscript/figures/figure1/figure1_panel.pdf"),
       device = "pdf", dpi = 300, width = 7, height = 5, units = "in")
ggsave(plot = figure1_panel_final,
       paste0(path_to_box, "/exposure_trajectories/",
              "manuscript/figures/figure1/figure1_panel.eps"),
       device = "eps", dpi = 300, width = 7, height = 5, units = "in")
ggsave(plot = figure1_panel_final,
       paste0(path_to_box, "/exposure_trajectories/",
              "manuscript/figures/figure1/figure1_panel.jpeg"),
       device = "jpeg", dpi = 300, width = 7, height = 5, units = "in")

#---- **figure 1 plot OLD ----
#using lemon package
ggplot(results_summary %>% filter(!Method == "LMM"), 
       aes(x = Beta, y = Percent, color = Method, shape = Method)) +
  geom_point(size = 2.25, position = position_dodge(-0.8)) + 
  scale_shape_manual(values = c(rep("square", (nrow(results_summary))))) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = .5, size = 0.75,
                position = position_dodge(-0.8)) +
  theme_minimal() + ylab("Percent Missing Data") +
  theme(legend.position = "bottom", legend.direction = "horizontal") + 
  scale_color_manual(values = cbPalette) + 
  scale_y_discrete(limits = rev(levels(results_summary$Percent))) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "dark grey", 
             size = 0.75) + 
  facet_rep_grid(rows = vars(Mechanism), cols = vars(Exposure), 
                 repeat.tick.labels = TRUE) + 
  geom_vline(data = truth, linetype = "dashed", aes(xintercept = Beta), 
             size = 0.75) + 
  xlab("\u03B2 (ln(hazard ratio))") + 
  theme(text = element_text(size = 9, color = "black"), 
        axis.text.x = element_text(color = "black"), 
        axis.text.y = element_text(color = "black")) + 
  theme(panel.border = element_blank(), axis.line = element_line(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(paste0(path_to_box, "/exposure_trajectories/",
              "manuscript/figures/figure1/effect_ests_mean_CI.jpeg"), 
       device = "jpeg", dpi = 300, width = 7, height = 5, units = "in")

ggsave(paste0(path_to_box, "/exposure_trajectories/",
              "submission/AJE/figures/figure1_effect_ests_mean_CI.eps"), 
       device = "eps", dpi = 300, width = 7, height = 5, units = "in")

#---- **efigure 7 plot ----
ggplot(results_summary, 
       aes(x = Beta, y = Percent, color = Method, shape = Method)) +
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = .5, size = 0.75,
                position = position_dodge(-0.8)) +
  geom_point(size = 2.25, position = position_dodge(-0.8)) + 
  scale_shape_manual(values = c(rep("square", 4), rep("asterisk", 1))) + 
  theme_minimal() + ylab("Percent Missing Data") +
  theme(legend.position = "bottom", legend.direction = "horizontal") + 
  scale_color_manual(values = cbPalette) + 
  scale_y_discrete(limits = rev(levels(results_summary$Percent))) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "dark grey", 
             size = 0.75) + 
  facet_grid(rows = vars(Mechanism), cols = vars(Exposure)) + 
  geom_vline(data = truth, aes(xintercept = Beta), size = 0.75) + 
  xlab("Beta (ln(hazard ratio))")

ggsave(paste0(path_to_box, "/exposure_trajectories/",
              "manuscript/figures/efigure7/effect_ests_mean_CI.jpeg"), 
       device = "jpeg", dpi = 300, width = 9, height = 7, units = "in")

#---- Figure 2 + eFigure 8: bias ----
#---- **read in data ----
bias_table <- 
  read_csv(paste0(path_to_box, "/exposure_trajectories/",
                  "manuscript/tables/table2/bias_plot_table.csv")) %>%
  set_colnames(c("Exposure", "Method", "Missing Percent", "Mechanism", "Bias"))

#rename FCS --> VTS
bias_table %<>% mutate("Method" = ifelse(Method == "FCS", "VTS", Method))
#rename JMVN --> NORM
bias_table %<>% mutate("Method" = ifelse(Method == "JMVN", "NORM", Method))

#---- **format data ----
methods <- c("CC", "NORM", "PMM", "VTS", "LMM")
bias_table$Method <- factor(bias_table$Method, levels = methods)
bias_table$Mechanism <- factor(bias_table$Mechanism, 
                               levels = c("MCAR", "MAR", "MNAR"))
bias_table$`Missing Percent` <- factor(bias_table$`Missing Percent`)

bias_table[which(bias_table$Exposure == "CES-D Wave 4"), "Exposure"] <- 
  "Elevated Baseline CES-D,\n1998"
bias_table[which(bias_table$Exposure == "CES-D Wave 9"), "Exposure"] <- 
  "Elevated End of Exposure CES-D, \n2008"
bias_table[which(bias_table$Exposure == "Elevated Average CES-D"), "Exposure"] <- 
  "Elevated Average CES-D,\n1998-2008"
bias_table[which(bias_table$Exposure == "Elevated CES-D Prop"), "Exposure"] <- 
  "Proportion Elevated CES-D,\n1998-2008"
bias_table$Exposure <- 
  factor(bias_table$Exposure, 
         levels = c("Elevated Baseline CES-D,\n1998", 
                    "Elevated End of Exposure CES-D, \n2008", 
                    "Elevated Average CES-D,\n1998-2008", 
                    "Proportion Elevated CES-D,\n1998-2008"))

#---- **figure 2 plot ----
#cowplot
plot_vars <- 
  expand.grid(unique(bias_table$Exposure), 
              unique(bias_table$Mechanism)) %>% 
  set_colnames(c("Exposure", "Mechanism")) %>% arrange(Mechanism) %>% 
  mutate("Label" = paste0(LETTERS[1:12], ")"))

plot_data <- bias_table %>% filter(!Method == "LMM")
plot_data$Percent <- str_remove(plot_data$`Missing Percent`, "%")
plot_data$Percent <- factor(plot_data$Percent)

figure2_plot_list <- list()

for(row in 1:nrow(plot_vars)){
  mech = plot_vars[row, "Mechanism"]
  exp = plot_vars[row, "Exposure"]
  label = plot_vars[row, "Label"]
  
  figure2_plot_list[[row]] <- 
    ggplot(plot_data %>% filter(Mechanism == mech & Exposure == exp), 
           aes(x = Percent, y = Bias, color = Method)) +
    geom_line(aes(group = Method), alpha = 0.75) + geom_point(alpha = 0.75) + 
    theme_bw() + 
    scale_y_continuous(limits = c(min(plot_data$Bias), max(plot_data$Bias))) +
    # ylim(c(min(plot_data$Bias), max(plot_data$Bias))) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "dark grey") +
    theme(text = element_text(size = 8, color = "black"), 
          axis.text = element_text(size = 8, color = "black"),
          panel.border = element_blank(), axis.line = element_line(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.position = "none") +
    scale_color_manual(values = cbPalette) + 
    ylab("Bias") + xlab("Missing Data, %")  + 
    theme(plot.margin = unit(c(t = 13, r = 1, b = 1, l = 13), unit = "pt"))
  
  if (row %in% c(11, 12)){
    # for plot K and L, label the far left X axis.
    figure2_plot_list[[row]] <- figure2_plot_list[[row]] +
      scale_y_continuous(limits = c(-1.5, 
                                    max(plot_data$Bias)),
                         breaks = c(-1.5, -1.0, -0.5, 0.0))
  }
  
  figure2_plot_list[[row]] <- 
    plot_grid(figure2_plot_list[[row]], labels = plot_vars$Label[[row]], 
              align = "vh", label_size = 8, hjust = 0, vjust = 1, 
              label_fontface = "plain")
  
  ggsave(plot = figure2_plot_list[[row]],
         filename = paste0(path_to_box, "/exposure_trajectories/",
                           "manuscript/figures/figure2/figure2_panel_",
                           substr(plot_vars$Label, 1, 1)[[row]], ".pdf"),
         device = "pdf", dpi = 300, width = 7/4, height = 5/3, units = "in")
  
  ggsave(plot = figure2_plot_list[[row]],
         filename = paste0(path_to_box, "/exposure_trajectories/",
                           "manuscript/figures/figure2/figure2_panel_",
                           substr(plot_vars$Label, 1, 1)[[row]], ".tiff"),
         device = "tiff", dpi = 300, width = 7/4, height = 5/3, units = "in")
}

figure2_plot_forlegend <- 
  ggplot(plot_data, 
         aes(x = Percent, y = Bias, color = Method)) +
  geom_line(aes(group = Method), alpha = 0.75) + geom_point(alpha = 0.75) + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) +
  theme(text = element_text(size = 8, color = "black"), 
        axis.text.x = element_text(color = "black"), 
        axis.text.y = element_text(color = "black"),
        legend.background = 
          element_rect(fill = "white", linetype = "solid", colour ="black"), 
        legend.title.align = 0.5,
        legend.position = "bottom", legend.direction = "horizontal") + 
  guides(shape = guide_legend(title = expression(underline(Method)), 
                              title.position = "top",
                              title.vjust = -1.3),
         color = guide_legend(title = expression(underline(Method)), 
                              title.position = "top",
                              title.vjust = -1.3)) +
  ylab("Bias") + xlab("Missing Data, %")

legend_b <- ggpubr::get_legend(figure2_plot_forlegend)
# plot(legend_b)

figure2_panel <- 
  plot_grid(plotlist = figure2_plot_list, align = "vh", 
            ncol = 4) + 
  theme(plot.margin = unit(c(t = 0, r = 0, b = 8, l = 0), unit = "pt"))
# figure2_panel

figure2_panel_final <- plot_grid(figure2_panel, 
                                 legend_b,
                                 ncol = 1, rel_heights = c(1, .1)) +
  theme(plot.margin = unit(c(t = 21, r = 10, b = 8, l = 0), unit = "pt")) +
  geom_text(data = data.frame(
    x = seq(0.13, 0.91, by = 0.26), y = rep(1, 4),
    label = paste0(levels(plot_vars$Exposure), "\n\n")),
    # This is stupid and I didn't figure out how to solve this without putting two \ns
    mapping = aes(x = x, y = y, label = label),
    size = 5/14*8, inherit.aes = FALSE) + 
  geom_text(data = data.frame(x = rep(1, 3), y = c(0.9, 0.6, 0.3),
                              label = paste0(levels(plot_vars$Mechanism), "\n")),
            mapping = aes(x = x, y = y, label = label),
            size = 5/14*8, angle = -90L, inherit.aes = FALSE)
figure2_panel_final

ggsave(plot = figure2_panel_final, 
       filename = paste0(path_to_box, "/exposure_trajectories/",
                         "manuscript/figures/figure2/figure2_panel.jpeg"), 
       device = "jpeg", dpi = 300, width = 7, height = 5, units = "in")
ggsave(plot = figure2_panel_final, 
       filename = paste0(path_to_box, "/exposure_trajectories/",
                         "manuscript/figures/figure2/figure2_panel.pdf"), 
       device = "pdf", dpi = 300, width = 7, height = 5, units = "in")
# Somehow can't save to EPS
ggsave(plot = figure2_panel_final, 
       filename = paste0(path_to_box, "/exposure_trajectories/",
                         "manuscript/figures/figure2/figure2_panel.tiff"), 
       device = "tiff", dpi = 300, width = 7, height = 5, units = "in")

#---- figure 2 OLD ----
ggplot(bias_table %>% filter(!Method == "LMM"), 
       mapping = aes(x = `Missing Percent`, y = Bias, 
                     color = Method)) +
  geom_line(aes(group = Method), alpha = 0.75) + geom_point(alpha = 0.75) + 
  theme_bw() + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "dark grey") +
  theme(text = element_text(size = 12, color = "black"), 
        axis.text = element_text(size = 12, color = "black"),
        panel.border = element_blank(), axis.line = element_line(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "bottom", legend.direction = "horizontal",
        legend.background = 
          element_rect(fill = "white", linetype = "solid", colour ="black"), 
        legend.title.align = 0.5) +
  scale_color_manual(values = cbPalette) + ylab("Bias") + 
  facetscales::facet_grid_sc(rows = vars(Mechanism), cols = vars(Exposure)) 

ggsave(paste0(path_to_box, "/exposure_trajectories/",
              "manuscript/figures/figure2/bias.jpeg"), 
       device = "jpeg", dpi = 300, width = 7, height = 5, units = "in")
# Somehow can't save to EPS
ggsave(paste0(path_to_box, "/exposure_trajectories/",
              "submission/AJE/figures/figure2_bias.tiff"), 
       device = "tiff", dpi = 300, width = 7, height = 5, units = "in")

#---- **efigure 8 plot ----
p_load("devtools")
devtools::install_github("zeehio/facetscales")

ggplot(bias_table, 
       mapping = aes(x = `Missing Percent`, y = Bias, color = Method, 
                     shape = Method)) +
  geom_line(aes(group = Method), alpha = 0.75, size = 0.75) +
  geom_point(alpha = 0.75, size = 2.25) + theme_bw() + 
  scale_shape_manual(values = c(rep("circle", 4), rep("asterisk", 1))) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "dark grey", 
             size = 0.75) +
  theme(legend.position = "bottom", legend.direction = "horizontal") + 
  scale_color_manual(values = cbPalette) + ylab("Bias") + 
  facetscales::facet_grid_sc(rows = vars(Mechanism), cols = vars(Exposure)) 

ggsave(paste0(path_to_box, "/exposure_trajectories/",
              "manuscript/figures/efigure8/bias.jpeg"), 
       device = "jpeg", dpi = 300, width = 9, height = 7, units = "in")

#---- Figure 3 + eFigure 9: RMSE ----
#---- **read in data ----
rmse_table <- read_csv(paste0(path_to_box, "/exposure_trajectories/",
                              "manuscript/tables/table3/rmse.csv"))
rmse_table %<>% 
  pivot_longer(cols = colnames(rmse_table)[grep("CES-D", 
                                                colnames(rmse_table))]) 
#rename FCS --> VTS
rmse_table %<>% mutate("Method" = ifelse(Method == "FCS", "VTS", Method))
#rename JMVN --> NORM
rmse_table %<>% mutate("Method" = ifelse(Method == "JMVN", "NORM", Method))

#---- **format data ----
methods <- c("CC", "NORM", "PMM", "VTS", "LMM")
rmse_table$Method <- factor(rmse_table$Method, levels = methods)
rmse_table$Mechanism <- factor(rmse_table$Mechanism, 
                               levels = c("MCAR", "MAR", "MNAR"))
rmse_table$`Missing Percent` <- factor(rmse_table$`Missing Percent`)

rmse_table[which(rmse_table$name == "CES-D Wave 4"), "name"] <- 
  "Elevated Baseline CES-D,\n1998"
rmse_table[which(rmse_table$name == "CES-D Wave 9"), "name"] <- 
  "Elevated End of Exposure CES-D, \n2008"
rmse_table[which(rmse_table$name == "Elevated Average CES-D"), "name"] <- 
  "Elevated Average CES-D,\n1998-2008"
rmse_table[which(rmse_table$name == "Elevated CES-D Prop"), "name"] <- 
  "Proportion Elevated CES-D,\n1998-2008"
rmse_table$name <- 
  factor(rmse_table$name, 
         levels = c("Elevated Baseline CES-D,\n1998", 
                    "Elevated End of Exposure CES-D, \n2008", 
                    "Elevated Average CES-D,\n1998-2008", 
                    "Proportion Elevated CES-D,\n1998-2008"))

#---- **figure 3 plot ----
#cowplot
plot_vars <- 
  expand.grid(unique(rmse_table$name),
              unique(rmse_table$Mechanism)) %>%
  set_colnames(c("name", "Mechanism")) %>%
  arrange("Mechanism") %>%
  mutate("Label" = paste0(LETTERS[1:12], ")"))

plot_data <- rmse_table %>% filter(!Method == "LMM")
plot_data$Percent <- str_remove(plot_data$`Missing Percent`, "%")
plot_data$Percent <- factor(plot_data$Percent)

figure3_plot_list <- list()

for(row in 1:nrow(plot_vars)){
  mech = plot_vars[row, "Mechanism"]
  exp = plot_vars[row, "name"]
  label = plot_vars[row, "Label"]
  
  figure3_plot_list[[row]] <-
    ggplot(plot_data %>% filter(Mechanism == mech & name == exp), 
           aes(x = Percent, y = value, color = Method)) +
    geom_line(aes(group = Method), alpha = 0.75) + geom_point(alpha = 0.75) + 
    theme_bw() + 
    scale_y_continuous(
      limits = case_when(mech %in% c("MCAR", "MAR") ~ c(0.00, 0.10),
                         mech == "MNAR" ~ c(0, 1.5)),
      breaks = case_when(mech %in% c("MCAR", "MAR") ~ seq(0.00, 0.10, 0.02),
                         mech == "MNAR" ~ seq(0, 1.5, 0.3))) +
    theme(text = element_text(size = 8, color = "black"), 
          axis.text = element_text(size = 8, color = "black"),
          panel.border = element_blank(), axis.line = element_line(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.position = "none") +
    scale_color_manual(values = cbPalette) + 
    ylab("RMSE") + xlab("Missing Data, %")  + 
    theme(plot.margin = unit(c(t = 13, r = 1, b = 1, l = 13), unit = "pt"))
  
  if(mech == "MNAR"){
    figure3_plot_list[[row]] <- figure3_plot_list[[row]] + 
      scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.5))
  }
  
  figure3_plot_list[[row]] <- 
    plot_grid(figure3_plot_list[[row]], labels = plot_vars$Label[[row]], 
              align = "vh", label_size = 8, hjust = 0, vjust = 1, 
              label_fontface = "plain")
  
  ggsave(plot = figure3_plot_list[[row]],
         filename = paste0(path_to_box, "/exposure_trajectories/",
                           "manuscript/figures/figure3/figure3_panel_",
                           substr(plot_vars$Label, 1, 1)[[row]], ".pdf"),
         device = "pdf", dpi = 300, width = 7/4, height = 5/3, units = "in")
  
  ggsave(plot = figure3_plot_list[[row]],
         filename = paste0(path_to_box, "/exposure_trajectories/",
                           "manuscript/figures/figure3/figure3_panel_",
                           substr(plot_vars$Label, 1, 1)[[row]], ".tiff"),
         device = "tiff", dpi = 300, width = 7/4, height = 5/3, units = "in")
}

figure3_plot_forlegend <- 
  ggplot(plot_data, 
         aes(x = Percent, y = value, color = Method)) +
  geom_line(aes(group = Method), alpha = 0.75) + geom_point(alpha = 0.75) + 
  theme_bw() + 
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  scale_color_manual(values = cbPalette) +
  theme(text = element_text(size = 8, color = "black"), 
        axis.text.x = element_text(color = "black"), 
        axis.text.y = element_text(color = "black"),
        legend.background = 
          element_rect(fill = "white", linetype = "solid", colour ="black"), 
        legend.title.align = 0.5) + 
  guides(shape = guide_legend(title = expression(underline(Method)), 
                              title.position = "top",
                              title.vjust = -1.3),
         color = guide_legend(title = expression(underline(Method)), 
                              title.position = "top",
                              title.vjust = -1.3)) +
  ylab("RMSE") + xlab("Missing Data, %")
# figure3_plot_forlegend

legend_b <- ggpubr::get_legend(figure3_plot_forlegend)
# plot(legend_b)

figure3_panel <- 
  plot_grid(plotlist = figure3_plot_list, align = "vh", 
            ncol = 4) + 
  theme(plot.margin = unit(c(t = 0, r = 0, b = 8, l = 0), unit = "pt"))
# figure3_panel

figure3_panel_final <- plot_grid(figure3_panel, 
                                 legend_b,
                                 ncol = 1, rel_heights = c(1, .1)) +
  theme(plot.margin = unit(c(t = 21, r = 10, b = 10, l = 0), unit = "pt")) +
  geom_text(data = data.frame(
    x = seq(0.13, 0.91, by = 0.26), y = rep(1, 4),
    label = paste0(levels(plot_vars$name), "\n\n")),
    # This is stupid and I didn't figure out how to solve this without putting two \ns
    mapping = aes(x = x, y = y, label = label),
    size = 5/14*8, inherit.aes = FALSE) + 
  geom_text(data = data.frame(x = rep(1, 3), y = c(0.9, 0.6, 0.3),
                              label = paste0(levels(plot_vars$Mechanism), "\n")),
            mapping = aes(x = x, y = y, label = label),
            size = 5/14*8, angle = -90L, inherit.aes = FALSE)
# figure3_panel_final

ggsave(plot = figure3_panel_final, 
       filename = paste0(path_to_box, "/exposure_trajectories/",
                         "manuscript/figures/figure3/figure3_panel.jpeg"), 
       device = "jpeg", dpi = 300, width = 7, height = 5, units = "in")
ggsave(plot = figure3_panel_final, 
       filename = paste0(path_to_box, "/exposure_trajectories/",
                         "manuscript/figures/figure3/figure3_panel.pdf"), 
       device = "pdf", dpi = 300, width = 7, height = 5, units = "in")
# Somehow can't save to EPS
ggsave(plot = figure3_panel_final, 
       filename = paste0(path_to_box, "/exposure_trajectories/",
                         "manuscript/figures/figure3/figure3_panel.tiff"), 
       device = "tiff", dpi = 300, width = 7, height = 5, units = "in")

#---- **figure 3 plot OLD ----
p_load("devtools")
devtools::install_github("zeehio/facetscales")

ggplot(rmse_table %>% filter(!Method == "LMM"), 
       mapping = aes(x = `Missing Percent`, y = value, 
                     color = Method)) +
  geom_line(aes(group = Method), alpha = 0.75, size = 0.75) +
  geom_point(alpha = 0.75, size = 2.25) + theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal") + 
  scale_color_manual(values = cbPalette) + ylab("RMSE") + 
  facetscales::facet_grid_sc(
    rows = vars(Mechanism), cols = vars(name), 
    scales = list(y = list(
      `MCAR` = scale_y_continuous(limits = c(0.00, 0.10), 
                                  breaks = seq(0.00, 0.10, 0.02)),
      `MAR` = scale_y_continuous(limits = c(0.00, 0.10), 
                                 breaks = seq(0.00, 0.10, 0.02)),
      `MNAR` = scale_y_continuous(limits = c(0, 1.5), 
                                  breaks = seq(0, 1.5, 0.5)))))

ggsave(paste0(path_to_box, "/exposure_trajectories/",
              "manuscript/figures/figure3/rmse.jpeg"), 
       device = "jpeg", dpi = 300, width = 9, height = 7, units = "in")
# Somehow can't save to EPS
ggsave(paste0(path_to_box, "/exposure_trajectories/",
              "submission/AJE/figures/figure3_rmse.tiff"), 
       device = "tiff", dpi = 300, width = 9, height = 7, units = "in")

#---- **efigure 9 plot ----
p_load("devtools")
devtools::install_github("zeehio/facetscales")

ggplot(rmse_table, 
       mapping = aes(x = `Missing Percent`, y = value, color = Method, 
                     shape = Method)) +
  geom_line(aes(group = Method), alpha = 0.75, size = 0.75) + 
  geom_point(alpha = 0.75, size = 2.25) + 
  scale_shape_manual(values = c(rep("circle", 4), rep("asterisk", 1))) + 
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal") + 
  scale_color_manual(values = cbPalette) + ylab("RMSE") + 
  facetscales::facet_grid_sc(
    rows = vars(Mechanism), cols = vars(name), 
    scales = list(y = list(
      `MCAR` = scale_y_continuous(limits = c(0.00, 0.10), 
                                  breaks = seq(0.00, 0.10, 0.02)),
      `MAR` = scale_y_continuous(limits = c(0.00, 0.10), 
                                 breaks = seq(0.00, 0.10, 0.02)),
      `MNAR` = scale_y_continuous(limits = c(0, 1.75), 
                                  breaks = seq(0, 1.75, 0.5)))))

ggsave(paste0(path_to_box, "/exposure_trajectories/",
              "manuscript/figures/efigure9/rmse.jpeg"), 
       device = "jpeg", dpi = 300, width = 9, height = 7, units = "in")

#---- Figure 4 + eFigure 10: runtimes ----
#---- **get filepaths ----
all_paths <- 
  list.files(path = paste0(path_to_box,
                           "/exposure_trajectories/data/hoffman_transfer/",
                           "results"), full.names = TRUE, pattern = "*.csv")

main_paths <- all_paths[!str_detect(all_paths, "sens")]

#---- **read in data ----
main_results <- do.call(rbind, lapply(main_paths, read_results)) %>% 
  #making sure only one copy of each seed
  na.omit() %>% group_by(Method, Exposure, Seed, Mechanism, Percent) %>% 
  slice_head(n = 1) %>% group_by(Method, Mechanism, Percent, Exposure)

#Rename FCS --> VTS
main_results %<>% mutate("Method" = ifelse(Method == "FCS", "VTS", Method))
#rename JMVN --> NORM
main_results %<>% mutate("Method" = ifelse(Method == "JMVN", "NORM", Method))

#double-checking
table(main_results$Mechanism, main_results$Percent, main_results$Method)/4

#---- **summarize data ----
main_run_times <- main_results %>% 
  group_by(Method) %>% summarize_at(.vars = c("Time"), ~mean(., na.rm = TRUE)) 

#---- **format data ----
methods <- c("CC", "NORM", "PMM", "VTS", "LMM")
main_results$Method <- factor(main_results$Method, 
                              levels = c("Truth", methods))
main_results$Mechanism <- 
  factor(main_results$Mechanism, levels = c("MCAR", "MAR", "MNAR"))
main_results$Percent <- factor(main_results$Percent)

main_results %<>% mutate("time_hours" = Time/60)

#---- **summary stats ----
main_results %>% group_by(Method, Percent) %>% 
  summarize_at(.vars = "time_hours", ~mean(.))

main_results$Percent <- str_remove(main_results$Percent, "%")
main_results$Percent <- factor(main_results$Percent)

#---- **figure 4 plot ----
ggplot(data = na.omit(main_results) %>% filter(!Method == "LMM"), 
       aes(x = Percent, y = time_hours, color = Method)) + 
  geom_boxplot() + ylim(c(0, 20)) +
  ylab("Computation Time, hours") + 
  xlab("Missing Data, %") + theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(text = element_text(size = 14, color = "black"), 
        axis.text = element_text(size = 14, color = "black"),
        panel.border = element_blank(), axis.line = element_line(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = c(0.10, 0.75), 
        legend.background = 
          element_rect(fill = "white", linetype = "solid", colour ="black"), 
        legend.title.align = 0.5,
        legend.text = element_text(size = 14)) + 
  guides(color = guide_legend(title = expression(underline(Method))))

ggsave(paste0(path_to_box, "/exposure_trajectories/",
              "manuscript/figures/figure4/figure4_run_times.jpeg"), 
       device = "jpeg", dpi = 300, width = 7, height = 5, units = "in")

ggsave(paste0(path_to_box, "/exposure_trajectories/",
              "manuscript/figures/figure4/figure4_run_times.pdf"), 
       device = "pdf", dpi = 300, width = 7, height = 5, units = "in")

ggsave(paste0(path_to_box, "/exposure_trajectories/",
              "manuscript/figures/figure4/figure4_run_times.eps"), 
       device = "eps", dpi = 300, width = 7, height = 5, units = "in")

# ggsave(paste0(path_to_box, "/exposure_trajectories/",
#               "submission/AJE/figures/figure4_run_times.eps"), 
#        device = "eps", dpi = 300, width = 7, height = 5, units = "in")

#---- **efigure 10 plot ----
ggplot(data = na.omit(main_results), 
       aes(x = Percent, y = time_hours, color = Method)) + 
  geom_boxplot(size = 0.75) + ylab("Computation Time (Hours)") + 
  xlab("Percent Missing Data") + theme_bw() + 
  theme(legend.position = "bottom", legend.direction = "horizontal") + 
  scale_color_manual(values = cbPalette)

ggsave(paste0(path_to_box, "/exposure_trajectories/",
              "manuscript/figures/efigure10/run_times.jpeg"), 
       device = "jpeg", dpi = 300, width = 9, height = 7, units = "in")

#---- eFigure 3: sensitivity results ----
#---- **get filepaths ----
all_paths <- 
  list.files(path = paste0(path_to_box,
                           "/exposure_trajectories/data/hoffman_transfer/",
                           "results"), full.names = TRUE, pattern = "*.csv")

sens_paths <- all_paths[str_detect(all_paths, "sens")]

#---- **read in data ----
sens_analyses <- do.call(rbind, lapply(sens_paths, read_results)) %>% 
  #making sure only one copy of each seed
  na.omit() %>% group_by(Method, Exposure, Seed, Mechanism, Percent) %>% 
  slice_head(n = 1) %>% group_by(Method, Mechanism, Percent, Exposure)

#Rename FCS --> VTS
sens_analyses %<>% mutate("Method" = ifelse(Method == "FCS", "VTS", Method))
#rename JMVN --> NORM
sens_analyses %<>% mutate("Method" = ifelse(Method == "JMVN", "NORM", Method))

#double-checking
table(sens_analyses$Mechanism, sens_analyses$Percent, sens_analyses$Method)/4

#---- **summarize data ----
results_summary <- sens_analyses %>% 
  group_by(Method, Mechanism, Percent, Exposure) %>%
  summarize_at(.vars = c("Beta", "SE", "LCI", "UCI", "Truth Capture"), 
               ~mean(., na.rm = TRUE)) 

#---- **read in truth table ----
truth_sens <- read_csv(paste0(path_to_box, 
                              "/exposure_trajectories/data/truth_sens.csv")) %>%
  dplyr::rename("LCI" = "LCI_beta", 
                "UCI" = "UCI_beta",
                "Beta" = "beta") %>%
  # mutate("Percent" = "0%", 
  #        "Truth Capture" = 1) 
  mutate(Exposure = 
           case_when(
             Exposure == "CES-D Wave 4" ~ "Elevated Baseline CES-D,\n1998",
             Exposure == "CES-D Wave 9" ~ "Elevated End of Exposure CES-D, \n2008",
             Exposure == "Elevated CES-D Prop" ~ 
               "Proportion Elevated CES-D,\n1998-2008",
             TRUE ~ "Elevated Average CES-D,\n1998-2008"))

truth_sens$Exposure <- 
  factor(truth_sens$Exposure, 
         levels = c("Elevated Baseline CES-D,\n1998", 
                    "Elevated End of Exposure CES-D, \n2008", 
                    "Elevated Average CES-D,\n1998-2008", 
                    "Proportion Elevated CES-D,\n1998-2008")) 

#---- **format data ----
methods <- c("CC", "NORM", "PMM", "VTS")
results_summary$Method <- 
  factor(results_summary$Method, levels = c("Truth", methods))
results_summary$Percent <- factor(results_summary$Percent)
results_summary$Mechanism <- 
  factor(results_summary$Mechanism, levels = c("MCAR", "MAR", "MNAR")) 

results_summary[which(results_summary$Exposure == "CES-D Wave 4"), 
                "Exposure"] <- "Elevated Baseline CES-D,\n1998"
results_summary[which(results_summary$Exposure == "CES-D Wave 9"), 
                "Exposure"] <- "Elevated End of Exposure CES-D, \n2008"
results_summary[which(results_summary$Exposure == "Elevated Average CES-D"), 
                "Exposure"] <- "Elevated Average CES-D,\n1998-2008"
results_summary[which(results_summary$Exposure == "Elevated CES-D Prop"), 
                "Exposure"] <- "Proportion Elevated CES-D,\n1998-2008"
results_summary$Exposure <- 
  factor(results_summary$Exposure, 
         levels = c("Elevated Baseline CES-D,\n1998", 
                    "Elevated End of Exposure CES-D, \n2008", 
                    "Elevated Average CES-D,\n1998-2008", 
                    "Proportion Elevated CES-D,\n1998-2008")) 

#---- **plot ----
ggplot(results_summary, 
       aes(x = Beta, y = Percent, color = Method, shape = Method)) +
  geom_point(size = 2.25, position = position_dodge(-0.8)) + 
  scale_shape_manual(values = c(rep("square", (nrow(results_summary))))) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = .5, size = 0.75,
                position = position_dodge(-0.8)) +
  theme_minimal() + ylab("Percent Missing Data") +
  theme(legend.position = "bottom", legend.direction = "horizontal") + 
  scale_color_manual(values = cbPalette) + 
  scale_y_discrete(limits = rev(levels(results_summary$Percent))) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "dark grey", 
             size = 0.75) + 
  facet_grid(rows = vars(Mechanism), cols = vars(Exposure)) + 
  geom_vline(data = truth_sens, aes(xintercept = Beta), size = 0.75) + 
  xlab("Beta (ln(hazard ratio))")

ggsave(paste0(path_to_box, "/exposure_trajectories/",
              "manuscript/figures/efigure3/effect_ests_mean_CI_sens.jpeg"), 
       device = "jpeg", dpi = 300, width = 9, height = 7, units = "in")

#---- eFigure 4: sensitivity bias ----
#---- **read in data ----
bias_table <- 
  read_csv(paste0(path_to_box, "/exposure_trajectories/",
                  "manuscript/tables/etable4/sens_bias_plot_table.csv")) %>%
  set_colnames(c("Exposure", "Method", "Missing Percent", "Mechanism", "Bias"))

#Rename FCS --> VTS
bias_table %<>% mutate("Method" = ifelse(Method == "FCS", "VTS", Method))
#rename JMVN --> NORM
bias_table %<>% mutate("Method" = ifelse(Method == "JMVN", "NORM", Method))

#---- **format data ----
methods <- c("CC", "NORM", "PMM", "VTS")
bias_table$Method <- factor(bias_table$Method, levels = methods)
bias_table$Mechanism <- factor(bias_table$Mechanism, 
                               levels = c("MCAR", "MAR", "MNAR"))
bias_table$`Missing Percent` <- factor(bias_table$`Missing Percent`)

bias_table[which(bias_table$Exposure == "CES-D Wave 4"), "Exposure"] <- 
  "Elevated Baseline CES-D,\n1998"
bias_table[which(bias_table$Exposure == "CES-D Wave 9"), "Exposure"] <- 
  "Elevated End of Exposure CES-D, \n2008"
bias_table[which(bias_table$Exposure == "Elevated Average CES-D"), "Exposure"] <- 
  "Elevated Average CES-D,\n1998-2008"
bias_table[which(bias_table$Exposure == "Elevated CES-D Prop"), "Exposure"] <- 
  "Proportion Elevated CES-D,\n1998-2008"
bias_table$Exposure <- 
  factor(bias_table$Exposure, 
         levels = c("Elevated Baseline CES-D,\n1998", 
                    "Elevated End of Exposure CES-D, \n2008", 
                    "Elevated Average CES-D,\n1998-2008", 
                    "Proportion Elevated CES-D,\n1998-2008"))

#---- **plot ----
p_load("devtools")
devtools::install_github("zeehio/facetscales")

ggplot(bias_table, 
       mapping = aes(x = `Missing Percent`, y = Bias, 
                     color = Method)) +
  geom_line(aes(group = Method), alpha = 0.75, size = 0.75) +
  geom_point(alpha = 0.75, size = 2.25) + theme_bw() + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "dark grey", 
             size = 0.75) +
  theme(legend.position = "bottom", legend.direction = "horizontal") + 
  scale_color_manual(values = cbPalette) + ylab("Bias") + 
  facetscales::facet_grid_sc(rows = vars(Mechanism), cols = vars(Exposure)) 

ggsave(paste0(path_to_box, "/exposure_trajectories/",
              "manuscript/figures/efigure4/bias_sens.jpeg"), 
       device = "jpeg", dpi = 300, width = 9, height = 7, units = "in")

#---- eFigure 5: sensitivity RMSE ----
#---- **read in data ----
rmse_table <- read_csv(paste0(path_to_box, "/exposure_trajectories/",
                              "manuscript/tables/etable3/rmse_sens.csv"))
rmse_table %<>% 
  pivot_longer(cols = colnames(rmse_table)[grep("CES-D", 
                                                colnames(rmse_table))]) 
#Rename FCS --> VTS
rmse_table %<>% mutate("Method" = ifelse(Method == "FCS", "VTS", Method))
#rename JMVN --> NORM
rmse_table %<>% mutate("Method" = ifelse(Method == "JMVN", "NORM", Method))

#---- **format data ----
methods <- c("CC", "NORM", "PMM", "VTS")
rmse_table$Method <- factor(rmse_table$Method, levels = methods)
rmse_table$Mechanism <- factor(rmse_table$Mechanism, 
                               levels = c("MCAR", "MAR", "MNAR"))
rmse_table$`Missing Percent` <- factor(rmse_table$`Missing Percent`)

rmse_table[which(rmse_table$name == "CES-D Wave 4"), "name"] <- 
  "Elevated Baseline CES-D,\n1998"
rmse_table[which(rmse_table$name == "CES-D Wave 9"), "name"] <- 
  "Elevated End of Exposure CES-D, \n2008"
rmse_table[which(rmse_table$name == "Elevated Average CES-D"), "name"] <- 
  "Elevated Average CES-D,\n1998-2008"
rmse_table[which(rmse_table$name == "Elevated CES-D Prop"), "name"] <- 
  "Proportion Elevated CES-D,\n1998-2008"
rmse_table$name <- 
  factor(rmse_table$name, 
         levels = c("Elevated Baseline CES-D,\n1998", 
                    "Elevated End of Exposure CES-D, \n2008", 
                    "Elevated Average CES-D,\n1998-2008", 
                    "Proportion Elevated CES-D,\n1998-2008"))

#---- **plot ----
p_load("devtools")
devtools::install_github("zeehio/facetscales")

ggplot(rmse_table, 
       mapping = aes(x = `Missing Percent`, y = value, 
                     color = Method)) +
  geom_line(aes(group = Method), alpha = 0.75, size = 0.75) + 
  geom_point(alpha = 0.75, size = 2.25) + theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal") + 
  scale_color_manual(values = cbPalette) + ylab("RMSE") + 
  facetscales::facet_grid_sc(
    rows = vars(Mechanism), cols = vars(name), 
    scales = list(y = list(
      `MCAR` = scale_y_continuous(limits = c(0, 0.10), 
                                  breaks = seq(0, 0.10, 0.02)),
      `MAR` = scale_y_continuous(limits = c(0, 0.10), 
                                 breaks = seq(0, 0.10, 0.02)),
      `MNAR` = scale_y_continuous(limits = c(0, 0.6), 
                                  breaks = seq(0, 0.6, 0.1)))))

ggsave(paste0(path_to_box, "/exposure_trajectories/",
              "manuscript/figures/efigure5/rmse_sens.jpeg"), 
       device = "jpeg", dpi = 300, width = 9, height = 7, units = "in")

#---- eFigure 6: sensitivity runtimes ----
#---- **get filepaths ----
all_paths <- 
  list.files(path = paste0(path_to_box,
                           "/exposure_trajectories/data/hoffman_transfer/",
                           "results"), full.names = TRUE, pattern = "*.csv")

sens_paths <- all_paths[str_detect(all_paths, "sens")]

#---- **read in data ----
sens_analyses <- do.call(rbind, lapply(sens_paths, read_results)) %>% 
  #making sure only one copy of each seed
  na.omit() %>% group_by(Method, Exposure, Seed, Mechanism, Percent) %>% 
  slice_head(n = 1) %>% group_by(Method, Mechanism, Percent, Exposure)

#Rename FCS --> VTS
sens_analyses %<>% mutate("Method" = ifelse(Method == "FCS", "VTS", Method))
#rename JMVN --> NORM
sens_analyses %<>% mutate("Method" = ifelse(Method == "JMVN", "NORM", Method))

#double-checking
table(sens_analyses$Mechanism, sens_analyses$Percent, sens_analyses$Method)/4

#---- **summarize data ----
sens_run_times <- sens_analyses %>% 
  group_by(Method) %>% summarize_at(.vars = c("Time"), ~mean(., na.rm = TRUE)) 

#---- **format data ----
methods <- c("CC", "NORM", "PMM", "VTS")
sens_analyses$Method <- factor(sens_analyses$Method, 
                               levels = c("Truth", methods))
sens_analyses$Mechanism <- 
  factor(sens_analyses$Mechanism, levels = c("MCAR", "MAR", "MNAR"))
sens_analyses$Percent <- factor(sens_analyses$Percent)

sens_analyses %<>% mutate("time_hours" = Time/60)

#---- **plot ----
ggplot(data = na.omit(sens_analyses), 
       aes(x = Percent, y = time_hours, color = Method)) + 
  geom_boxplot(size = 0.75) + ylab("Computation Time (Hours)") + 
  xlab("Percent Missing Data") + theme_bw() + 
  theme(legend.position = "bottom", legend.direction = "horizontal") + 
  scale_color_manual(values = cbPalette)

ggsave(paste0(path_to_box, "/exposure_trajectories/",
              "manuscript/figures/efigure6/run_times_sens.jpeg"), 
       device = "jpeg", dpi = 300, width = 7, height = 5, units = "in")
