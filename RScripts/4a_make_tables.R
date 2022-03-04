#---- package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse")

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


