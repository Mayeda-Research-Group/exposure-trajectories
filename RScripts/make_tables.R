#---- package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "vroom")

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
#---- **read in true data ----
truth <- read_csv(paste0(path_to_dropbox, 
                         "/exposure_trajectories/manuscript/tables/", 
                         "results_CC_1000.csv")) %>% 
  dplyr::select(-one_of("people_dropped")) %>% 
  filter(Method == "Truth") %>% group_by(Exposure) %>% slice(1)

#---- **read in Hoffman data ----
methods <- c("JMVN", "PMM", "FCS")

for(method in methods){
  if(!exists("results")){
    results <- 
      vroom(paste0(path_to_dropbox, "/exposure_trajectories/data/",
                   "hoffman_transfer/", method, "/", 
                   list.files(path = 
                                paste0(path_to_dropbox, 
                                       "/exposure_trajectories/data/", 
                                       "hoffman_transfer/", method))), 
            col_names = FALSE) %>% 
      set_colnames(c("Exposure", "beta", "SD", "LCI", "UCI", "Method", 
                     "Missingness", "Type", "truth_capture")) %>% 
      group_by(Method, Type, Missingness, Exposure) %>% slice_head(n = 100)
  } else{
    results <- 
      rbind(results, 
            vroom(paste0(path_to_dropbox, 
                         "/exposure_trajectories/data/", 
                         "hoffman_transfer/", method, "/", 
                         list.files(path = paste0(path_to_dropbox, 
                                                  "/exposure_trajectories/data/", 
                                                  "hoffman_transfer/", method))), 
                  col_names = FALSE) %>% 
              set_colnames(c("Exposure", "beta", "SD", "LCI", "UCI", "Method", 
                             "Missingness", "Type", "truth_capture")) %>% 
              group_by(Method, Type, Missingness, Exposure) %>% 
              slice_head(n = 100))
  }
}

#---- RMSE table ----
#---- **table shell ----
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



