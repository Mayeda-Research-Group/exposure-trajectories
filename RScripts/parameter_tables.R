#This script creates all the parameter tables needed for the simulations:
#-- optimized betas for masking function
#-- 

#---- package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "broom", "ResourceSelection", 
       "survival")

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
           col_types = cols(HHIDPN = col_character())) %>% as.data.frame() 

#---- masking: optimized betas table ----
#don't need MCAR
mechanisms <- c("MAR", "MNAR")
percents <- c(10, 20, 30)

beta_0_table <- expand_grid(mechanisms, percents) %>% 
  mutate("beta0" = 0)

for(mechanism in mechanisms){
  for(percent in percents){
    optimized <- 
      read_rds(file = paste0(path_to_dropbox, "/exposure_trajectories/data/", 
                             "optimized_masking_intercepts/optim_", mechanism, 
                             percent, ".RDS"))
    
    beta_0_table[which(beta_0_table$mechanisms == mechanism & 
                         beta_0_table$percents == percent), "beta0"] <- 
      optimized$minimum
  }
}

write_csv(beta_0_table, paste0(path_to_dropbox, "/exposure_trajectories/data/", 
                               "beta_0_table.csv"))

#---- masking: beta matrix ----
beta_mat <- #effect sizes
  matrix(c(log(1.10), log(1.05), log(1.05), log(1.25), log(1.10), log(1.25)), 
         nrow = 1) %>% 
  #MAR
  set_colnames(c("cesdpre", "condepre", "cesdpre_condepre",
                 #MNAR
                 "death2018", "cesdcurrent", "death2018_cesdcurrent")) %>% 
  set_rownames(c("beta"))

write_csv(as.data.frame(beta_mat), 
          paste0(path_to_dropbox, "/exposure_trajectories/data/beta_mat.csv"))