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

#---- Table 2 ----
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
                   "Percent", "Mechanism", "Truth Capture", "Time"))
}

main_results <- do.call(rbind, lapply(main_paths, read_results)) %>% na.omit()

#---- **table shell ----
methods <- c("CC", "JMVN", "PMM", "FCS")

rmse_table <- 
  data.frame("Method" = rep(unique(results$Method), 
                            each = length(unique(results$Type))*
                              length(unique(results$Missingness))), 
             "Mechanism" = rep(c("MCAR", "MAR", "MNAR"), 
                               each = length(unique(results$Missingness))), 
             "Missing Percent" = rep(unique(results$Missingness), 
                                     length(unique(results$Type))))
rmse_table <- cbind(rmse_table, matrix(nrow = nrow(rmse_table), ncol = 4)) %>% 
  set_colnames(c("Method", "Mechanism", "Missing Percent", 
                 unique(results$Exposure)))

#---- **calculate RMSE ----
for(i in 1:nrow(rmse_table)){
  method = rmse_table[i, "Method"]
  mechanism = rmse_table[i, "Mechanism"]
  percent = rmse_table[i, "Missing Percent"]
  
  for(exposure in unique(results$Exposure)){
    subset <- results %>% filter(Exposure == exposure, Method == method, 
                                 Type == mechanism, Missingness == percent)
    
    rmse_table[i, exposure] <- round(sqrt(mean((subset$beta - 
                                                  truth[[which(truth$Exposure == exposure), "beta"]])^2)), 3)
  }
}

#---- **save results ----
write_csv(rmse_table, paste0(path_to_dropbox, "/exposure_trajectories/",
                             "manuscript/tables/rmse.csv"))



