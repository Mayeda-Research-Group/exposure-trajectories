# Aim: Make tables for results
# Authors: Crystal Shaw, Yingyan Wu

#---- package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "labelled", "gtsummary", "writexl")

#No scientific notation
options(scipen = 999)

#---- note ----
# Since the difference between win and OS, put substituted directory here
# Yingyan's directory: C:/Users/yingyan_wu
#                      C:/Users/yingyan_wu/Dropbox
# Crystal's directory: /Users/CrystalShaw
#                     ~/Dropbox/Projects

#Changing directories here will change them throughout the script
path_to_box <- "C:/Users/Yingyan Wu"
path_to_dropbox <- "C:/Users/Yingyan Wu/Dropbox"

#---- Table 1: characteristics description ----
#---- **read in analytical sample ----
CESD_data_wide <- 
  read_csv(paste0(path_to_dropbox, 
                  "/exposure_trajectories/data/", 
                  "CESD_data_wide.csv"), 
           col_types = cols(.default = col_double()))

# colnames(CESD_data_wide)

#---- **table shell ----
# Showing 2 digits for percentages for categorical variables)
options(
  gtsummary.tbl_summary.percent_fun = function(x) sprintf(x * 100, fmt='%#.2f'))

table_1 <- CESD_data_wide %>% 
  mutate("raceeth" = case_when(hispanic == 0 & white == 1 ~ "Non-hispanic White",
                               hispanic == 0 & black == 1 ~ "Non-hispanic Black",
                               hispanic == 1 ~ "Hispanic",
                               other == 1 ~ "Other"),
         "r4mstat_cat" = case_when(
           r4married_partnered == 1 ~ "Married/Partnered",
           r4not_married_partnered == 1 ~ "Not Married/Partnered",
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
    r4memrye_impute = "Memory problems (%)",
    r4conde_impute = "Self-reported chronic condition count",
    r4BMI = "BMI",
    r4drinking_cat = "Drinking behavior (%)",
    smoker = "Smoking status (%)",
    r4cesd_elevated = "Baseline Elevated CES-D"
  ) %>%
  labelled::drop_unused_value_labels() %>%
  labelled::set_value_labels(female = c("Yes" = 1, "No" = 0),
                             ed_cat = c("Less than High School" = 1, "High School" = 2, 
                                        "Some college" = 3, "Bachelor's degree" = 4, 
                                        "Graduate studies" = 5),
                             r4hibpe_impute = c("Yes" = 1, "No" = 0),
                             r4diabe_impute = c("Yes" = 1, "No" = 0),
                             r4hearte_impute = c("Yes" = 1, "No" = 0),
                             r4stroke_impute = c("Yes" = 1, "No" = 0),
                             r4cancre_impute = c("Yes" = 1, "No" = 0),
                             r4lunge_impute = c("Yes" = 1, "No" = 0),
                             r4memrye_impute = c("Yes" = 1, "No" = 0),
                             r4drinking_cat = c("Heavy Drinking" = 2,
                                                "Moderate Drinking" = 1,
                                                "No Drinking" = 0),
                             smoker = c("Ever smoke" = 1, "Never smoke" = 0),
                             r4cesd_elevated = c("Elevated CES-D" = 1, 
                                                 "Not Elevated CES-D" = 0)) %>%
  modify_if(is.labelled, to_factor) %>%
  tbl_summary(statistic = list(all_categorical() ~ "{n} ({p})", 
                               all_continuous() ~ "{mean} ({sd})"),
              type = list(c("r4conde_impute") ~ "continuous"),
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
write_xl::write_xlsx(as.tibble(table_1), 
           paste0(path_to_dropbox,
                  "/exposure_trajectories/manuscript/tables/", 
                  "table_1_temp.xlsx"))

#---- Table 2: RMSE ----
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
read_results <- function(paths){
  data.table::fread(paths, fill = TRUE) %>% na.omit() %>%
    set_colnames(c("Exposure", "Beta", "SE", "LCI", "UCI", "Method",
                   "Percent", "Mechanism", "Truth Capture", "Time", "Seed"))
}

# test <- data.table::fread(main_paths[2], fill = TRUE) %>%
#   set_colnames(c("Exposure", "Beta", "SE", "LCI", "UCI", "Method",
#                  "Percent", "Mechanism", "Truth Capture", "Time")) %>% na.omit()

main_results <- do.call(rbind, lapply(main_paths, read_results)) %>% na.omit()

#---- **limit runs for table (for now) ----
main_results %<>% 
  group_by(Method, Mechanism, Percent, Exposure) %>% slice_head(n = 1000) %>% 
  na.omit()

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
                             "manuscript/tables/table2/rmse.csv"))

#---- eTable 2: sensitivity RMSE ----
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
read_results <- function(paths){
  data.table::fread(paths, fill = TRUE) %>% na.omit() %>%
    set_colnames(c("Exposure", "Beta", "SE", "LCI", "UCI", "Method",
                   "Percent", "Mechanism", "Truth Capture", "Time", "Seed"))
}

sens_analyses <- do.call(rbind, lapply(sens_paths, read_results)) %>% na.omit()

#---- **limit runs for table (for now) ----
sens_analyses %<>% 
  group_by(Method, Mechanism, Percent, Exposure) %>% slice_head(n = 1000) %>% 
  na.omit()

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
                             "manuscript/tables/etable2/rmse_sens.csv"))



