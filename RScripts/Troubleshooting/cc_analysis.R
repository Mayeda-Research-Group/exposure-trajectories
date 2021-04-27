#---- package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "broom", "ResourceSelection", 
       "survival", "openxlsx", "lubridate", "lme4")

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



