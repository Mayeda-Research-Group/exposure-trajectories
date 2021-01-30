#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "knitr", "fmsb")

#No scientific notation
options(scipen = 999)

set.seed(20200819)

#---- color palette ----
light_blue <- "#B4DAE5FF"

#---- note ----
# Since the difference between win and OS, put substituted directory here
# Yingyan's directory: C:/Users/yingyan_wu
#                      C:/Users/yingyan_wu/Dropbox
# Crystal's directory: /Users/CrystalShaw
#                     ~/Dropbox/Projects

#Changing directories here will change them throughout the script
path_to_box <- "C:/Users/yingyan_wu"
path_to_dropbox <- "C:/Users/yingyan_wu/Dropbox"

#---- read in dataset ----
CESD_data_wide <- read_csv(paste0(path_to_dropbox,
                                  "/exposure_trajectories/data/",
                                  "CESD_data_wide.csv"))

#---- KAPPA test ----

fmsb::Kappa.test(
  table(CESD_data_wide$r4cesd_elevated, CESD_data_wide$r9cesd_elevated))$
  Result$estimate

fmsb::Kappa.test(
  table(CESD_data_wide$r4cesd_elevated, CESD_data_wide$avg_cesd_elevated))$
  Result$estimate

fmsb::Kappa.test(
  table(CESD_data_wide$r9cesd_elevated, CESD_data_wide$avg_cesd_elevated))$
  Result$estimate

