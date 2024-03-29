# Aim: Make tables for results
# Authors: Crystal Shaw, Yingyan Wu
# 
# Input: CESD_data_wide.csv, truth(_sens).csv, all sims results
# Output: Table 1, Table 2 (RMSE), eTable 2 (sensitivity RMSE)

#---- package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "labelled", "gtsummary", "writexl")

#No scientific notation
options(scipen = 999)

#---- note ----
# Since the difference between win and OS, put substituted directory here
# Yingyan's directory: C:/Users/yingyan_wu/Dropbox
# Crystal's directory: ~/Dropbox/Projects

#Changing directories here will change them throughout the script
path_to_dropbox <- "C:/Users/yingyan_wu/Dropbox"

#---- source scripts ----
source(here::here("RScripts", "functions", "read_results.R"))

#---- manuscript calcs ----
#---- **read in analytic sample ----
CESD_data_wide <- 
  read_csv(paste0(path_to_dropbox, 
                  "/exposure_trajectories/data/", 
                  "CESD_data_wide.csv"), 
           col_types = cols(.default = col_double()))

#---- **baseline CES-D ----
sum(CESD_data_wide$r4cesd_elevated)
mean(CESD_data_wide$r4cesd_elevated)
mean(CESD_data_wide$r4cesd_elevated_sens)

#---- Table 1: characteristics description ----
# colnames(CESD_data_wide)

#---- **table shell ----
# Showing 2 digits for percentages for categorical variables)
options(
  gtsummary.tbl_summary.percent_fun = function(x) sprintf(x * 100, fmt='%#.1f'))

table_1 <- CESD_data_wide %>% 
  mutate("raceeth" = case_when(hispanic == 0 & white == 1 ~ "Non-Hispanic White",
                               hispanic == 0 & black == 1 ~ "Non-Hispanic Black",
                               hispanic == 1 ~ "Hispanic",
                               other == 1 ~ "Other"),
         "r4mstat_cat" = case_when(
           r4married_partnered == 1 ~ "Married/Partnered",
           r4not_married_partnered == 1 ~ "Not married/partnered",
           r4widowed == 1 ~ "Widowed")) %>%
  select("r4age_y_int","female", "raceeth", "ed_cat", "r4mstat_cat",
         "r4hibpe_impute", "r4diabe_impute", "r4hearte_impute",
         "r4stroke_impute", "r4cancre_impute", "r4lunge_impute",
         "r4memrye_impute", "r4conde_impute", "r4BMI", "r4drinking_cat",
         "smoker", "r4cesd_elevated") %>%
  labelled::set_variable_labels(
    r4age_y_int = "Baseline age in years",
    female = "Female (%)",
    raceeth = "Race/Ethnicity (%)",
    ed_cat = "Education level (%)",
    r4mstat_cat = "Marital status (%)",
    r4hibpe_impute = "Hypertension(Diagnosed) (%)",
    r4diabe_impute = "Diabetes (%)",
    r4hearte_impute = "Heart disease (%)",
    r4stroke_impute = "Stroke (%)",
    r4cancre_impute = "Cancer (%)",
    r4lunge_impute = "Lung disease (%)",
    r4memrye_impute = "Self-reported memory problems (%)",
    r4conde_impute = "Chronic condition count",
    r4BMI = "BMI",
    r4drinking_cat = "Drinking behavior (%)",
    smoker = "Smoking status (%)",
    r4cesd_elevated = "Baseline Elevated CES-D"
  ) %>%
  labelled::drop_unused_value_labels() %>%
  labelled::set_value_labels(female = c("Yes" = 1, "No" = 0),
                             ed_cat = c("Less than high school" = 1, "High school" = 2, 
                                        "Some college" = 3, "Bachelor's degree" = 4, 
                                        "Graduate studies" = 5),
                             r4hibpe_impute = c("Yes" = 1, "No" = 0),
                             r4diabe_impute = c("Yes" = 1, "No" = 0),
                             r4hearte_impute = c("Yes" = 1, "No" = 0),
                             r4stroke_impute = c("Yes" = 1, "No" = 0),
                             r4cancre_impute = c("Yes" = 1, "No" = 0),
                             r4lunge_impute = c("Yes" = 1, "No" = 0),
                             r4memrye_impute = c("Yes" = 1, "No" = 0),
                             r4drinking_cat = c("Heavy alcohol use" = 2,
                                                "Moderate alcohol use" = 1,
                                                "No alcohol use" = 0),
                             smoker = c("Ever smoked" = 1, "Never smoked" = 0),
                             r4cesd_elevated = c("Elevated CES-D" = 1, 
                                                 "Not Elevated CES-D" = 0)) %>%
  modify_if(is.labelled, to_factor) %>%
  tbl_summary(statistic = list(all_categorical() ~ "{n} ({p})", 
                               all_continuous() ~ "{mean} ({sd})"),
              type = c("r4conde_impute") ~ "continuous",
              by = r4cesd_elevated,
              missing = "no") %>%
  add_overall %>%
  modify_header(label = "") %>%
  modify_spanning_header(starts_with("stat_") ~ "**Baseline CES-D**") %>%
  bold_labels() %>%
  modify_footnote(
    update = list(
      starts_with("stat_") ~
        "Mean (SD) for continuous variables; n (%) for categorical variables

      Abbreviations: BMI, body mass index; CES-D, Center for Epidemiologic Studies-Depression")
  )

table_1 %>% as_flex_table()

#---- **save the table ----
# The table will be formatted in a linked excel
writexl::write_xlsx(as_tibble(table_1), 
                     paste0(path_to_dropbox,
                            "/exposure_trajectories/manuscript/tables/table1/", 
                            "table_1_temp.xlsx"))

#---- Table (Figure 4): RMSE ----
#---- **read in truth table ----
truth <- read_csv(paste0(path_to_dropbox, 
                         "/exposure_trajectories/data/", "truth.csv")) %>%
  dplyr::rename("LCI" = "LCI_beta", 
                "UCI" = "UCI_beta",
                "Beta" = "beta") %>% 
  mutate("Percent" = "0%", 
         "Truth Capture" = 1)

#---- **get filepaths ----
all_paths <- 
  list.files(path = paste0(path_to_dropbox,
                           "/exposure_trajectories/data/hoffman_transfer/",
                           "results"), full.names = TRUE, pattern = "*.csv")

main_paths <- all_paths[!str_detect(all_paths, "sens")]

#---- **read in data ----
main_results <- do.call(rbind, lapply(main_paths, read_results)) %>% 
  #making sure only one copy of each seed
  na.omit() %>% group_by(Method, Exposure, Seed, Mechanism, Percent) %>% 
  slice_head(n = 1) %>% group_by(Method, Mechanism, Percent, Exposure)

#double-checking
table(main_results$Mechanism, main_results$Percent, main_results$Method)/4

#---- **table shell ----
rmse_table <- 
  data.frame("Method" = rep(unique(main_results$Method), 
                            each = length(unique(main_results$Mechanism))*
                              length(unique(main_results$Percent))), 
             "Mechanism" = rep(c("MCAR", "MAR", "MNAR"), 
                               each = length(unique(main_results$Percent))), 
             "Missing Percent" = rep(unique(main_results$Percent), 
                                     length(unique(main_results$Mechanism))))

rmse_table <- cbind(rmse_table, matrix(nrow = nrow(rmse_table), ncol = 4)) %>% 
  set_colnames(c("Method", "Mechanism", "Missing Percent", 
                 unique(main_results$Exposure)))

#---- **calculate RMSE ----
for(i in 1:nrow(rmse_table)){
  method = rmse_table[i, "Method"]
  mechanism = rmse_table[i, "Mechanism"]
  percent = rmse_table[i, "Missing Percent"]
  
  for(exposure in unique(main_results$Exposure)){
    subset <- main_results %>% 
      filter(Exposure == exposure, Method == method, Mechanism == mechanism, 
             Percent == percent)
    
    rmse_table[i, exposure] <- 
      round(sqrt(mean((subset$Beta - 
                         truth[[which(truth$Exposure == exposure), 
                                "Beta"]])^2)), 3)
  }
}

#---- **save results ----
write_csv(rmse_table, paste0(path_to_dropbox, "/exposure_trajectories/",
                             "manuscript/tables/table3/rmse.csv"))

#---- eTable 1: percent agreement ----
#---- **read in data ----
complete_sample <- read_csv(paste0(path_to_dropbox,
                                   "/exposure_trajectories/data/",
                                   "CESD_data_wide.csv"))

#---- **compute agreement ----
sum(diag(table(complete_sample$r4cesd_elevated, 
               complete_sample$r9cesd_elevated)))/nrow(dataset)*100

sum(diag(table(complete_sample$r4cesd_elevated, 
               complete_sample$avg_cesd_elevated)))/nrow(dataset)*100

sum(diag(table(complete_sample$r9cesd_elevated, 
               complete_sample$avg_cesd_elevated)))/nrow(dataset)*100

#---- eTable 3: bias ----
#---- **read in truth table ----
truth <- read_csv(paste0(path_to_dropbox, 
                         "/exposure_trajectories/data/", "truth.csv")) %>%
  dplyr::rename("LCI" = "LCI_beta", 
                "UCI" = "UCI_beta",
                "Beta" = "beta") %>% 
  mutate("Percent" = "0%", 
         "Truth Capture" = 1)

#---- **get filepaths ----
all_paths <- 
  list.files(path = paste0(path_to_dropbox,
                           "/exposure_trajectories/data/hoffman_transfer/",
                           "results"), full.names = TRUE, pattern = "*.csv")

main_paths <- all_paths[!str_detect(all_paths, "sens")]

#---- **read in data ----
main_results <- do.call(rbind, lapply(main_paths, read_results)) %>% 
  #making sure only one copy of each seed
  na.omit() %>% group_by(Method, Exposure, Seed, Mechanism, Percent) %>% 
  slice_head(n = 1) %>% group_by(Method, Mechanism, Percent, Exposure)

#double-checking
table(main_results$Mechanism, main_results$Percent, main_results$Method)/4

#---- **join data ----
bias_plot_table <- left_join(main_results, truth, by = c("Exposure")) %>% 
  mutate("bias" = `Beta.x` - `Beta.y`) %>% 
  group_by(Exposure, `Method.x`, `Percent.x`, Mechanism) %>% 
  summarize_at("bias", mean) 

bias_table <- bias_plot_table %>% 
  pivot_wider(names_from = Exposure, values_from = bias) %>% 
  set_colnames(c("Method", "Percent", "Mechanism", "Baseline CES-D", 
                 "End of Exposure CES-D", "Elevated Average CES-D", 
                 "Proportion Elevated CES-D")) %>% 
  dplyr::select("Mechanism", "Method", "Percent", everything()) 

bias_table$Mechanism <- 
  factor(bias_table$Mechanism, levels = c("MCAR", "MAR", "MNAR"))
bias_table$Method <- 
  factor(bias_table$Method, levels = c("CC", "JMVN", "PMM", "FCS"))

bias_table %<>% arrange(Mechanism, Method)

#---- **save results ----
write_csv(bias_plot_table, 
          paste0(path_to_dropbox, "/exposure_trajectories/manuscript/tables/", 
                 "table2/bias_plot_table.csv"))

write_csv(bias_table, 
          paste0(path_to_dropbox, "/exposure_trajectories/manuscript/tables/", 
                 "table2/bias_table.csv"))

#---- eTable 4: sensitivity bias ----
#---- **read in truth table ----
truth_sens <- 
  read_csv(paste0(path_to_dropbox, 
                  "/exposure_trajectories/data/", "truth_sens.csv")) %>%
  dplyr::rename("LCI" = "LCI_beta", 
                "UCI" = "UCI_beta",
                "Beta" = "beta") %>% 
  mutate("Percent" = "0%", 
         "Truth Capture" = 1)

#---- **get filepaths ----
all_paths <- 
  list.files(path = paste0(path_to_dropbox,
                           "/exposure_trajectories/data/hoffman_transfer/",
                           "results"), full.names = TRUE, pattern = "*.csv")

sens_paths <- all_paths[str_detect(all_paths, "sens")]

#---- **read in data ----
sens_analysis <- do.call(rbind, lapply(sens_paths, read_results)) %>% 
  #making sure only one copy of each seed
  na.omit() %>% group_by(Method, Exposure, Seed, Mechanism, Percent) %>% 
  slice_head(n = 1) %>% group_by(Method, Mechanism, Percent, Exposure)

#double-checking
table(sens_analysis$Mechanism, sens_analysis$Percent, sens_analysis$Method)/4

#---- **join data ----
bias_plot_table <- left_join(sens_analysis, truth_sens, by = c("Exposure")) %>% 
  mutate("bias" = `Beta.x` - `Beta.y`) %>% 
  group_by(Exposure, `Method.x`, `Percent.x`, Mechanism) %>% 
  summarize_at("bias", mean) 

bias_table <- bias_plot_table %>% 
  pivot_wider(names_from = Exposure, values_from = bias) %>% 
  set_colnames(c("Method", "Percent", "Mechanism", "Baseline CES-D", 
                 "End of Exposure CES-D", "Elevated Average CES-D", 
                 "Proportion Elevated CES-D")) %>% 
  dplyr::select("Mechanism", "Method", "Percent", everything()) 

bias_table$Mechanism <- 
  factor(bias_table$Mechanism, levels = c("MCAR", "MAR", "MNAR"))
bias_table$Method <- 
  factor(bias_table$Method, levels = c("CC", "JMVN", "PMM", "FCS"))

bias_table %<>% arrange(Mechanism, Method)

#---- **save results ----
write_csv(bias_plot_table, 
          paste0(path_to_dropbox, "/exposure_trajectories/manuscript/tables/", 
                 "etable4/sens_bias_plot_table.csv"))

write_csv(bias_table, 
          paste0(path_to_dropbox, "/exposure_trajectories/manuscript/tables/", 
                 "etable4/sens_bias_table.csv"))

#---- eTable 3: sensitivity RMSE ----
#---- **read in truth table ----
truth_sens <- 
  read_csv(paste0(path_to_dropbox, 
                  "/exposure_trajectories/data/", "truth_sens.csv")) %>%
  dplyr::rename("LCI" = "LCI_beta", 
                "UCI" = "UCI_beta",
                "Beta" = "beta") %>% 
  mutate("Percent" = "0%", 
         "Truth Capture" = 1)

#---- **get filepaths ----
all_paths <- 
  list.files(path = paste0(path_to_dropbox,
                           "/exposure_trajectories/data/hoffman_transfer/",
                           "results"), full.names = TRUE, pattern = "*.csv")

sens_paths <- all_paths[str_detect(all_paths, "sens")]

#---- **read in data ----
sens_analyses <- do.call(rbind, lapply(sens_paths, read_results)) %>% 
  #making sure only one copy of each seed
  na.omit() %>% group_by(Method, Exposure, Seed) %>% slice_head(n = 1) %>% 
  group_by(Method, Mechanism, Percent, Exposure)

#double-checking
table(sens_analyses$Mechanism, sens_analyses$Percent, sens_analyses$Method)/4

#---- **table shell ----
rmse_table <- 
  data.frame("Method" = rep(unique(sens_analyses$Method), 
                            each = length(unique(sens_analyses$Mechanism))*
                              length(unique(sens_analyses$Percent))), 
             "Mechanism" = rep(c("MCAR", "MAR", "MNAR"), 
                               each = length(unique(sens_analyses$Percent))), 
             "Missing Percent" = rep(unique(sens_analyses$Percent), 
                                     length(unique(sens_analyses$Mechanism))))

rmse_table <- cbind(rmse_table, matrix(nrow = nrow(rmse_table), ncol = 4)) %>% 
  set_colnames(c("Method", "Mechanism", "Missing Percent", 
                 unique(sens_analyses$Exposure)))

#---- **calculate RMSE ----
for(i in 1:nrow(rmse_table)){
  method = rmse_table[i, "Method"]
  mechanism = rmse_table[i, "Mechanism"]
  percent = rmse_table[i, "Missing Percent"]
  
  for(exposure in unique(sens_analyses$Exposure)){
    subset <- sens_analyses %>% 
      filter(Exposure == exposure, Method == method, Mechanism == mechanism, 
             Percent == percent)
    
    rmse_table[i, exposure] <- 
      round(sqrt(mean((subset$Beta - 
                         truth_sens[[which(truth_sens$Exposure == exposure), 
                                     "Beta"]])^2)), 3)
  }
}

#---- **save results ----
write_csv(rmse_table, paste0(path_to_dropbox, "/exposure_trajectories/",
                             "manuscript/tables/etable3/rmse_sens.csv"))



