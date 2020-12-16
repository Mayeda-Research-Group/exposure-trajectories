#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "mice", "broom", "ghibli", 
       "ResourceSelection", "survival", "openxlsx")

#No scientific notation
options(scipen = 999)

set.seed(20200819)

#---- note ----
# Since the difference between win and OS, put substituted directory here
# Yingyan's directory: C:/Users/yingyan_wu
#                      C:/Users/yingyan_wu/Dropbox
# Crystal's directory: /Users/CrystalShaw
#                     ~/Dropbox/Projects

#Changing directories here will change them throughout the script
path_to_box <- "/Users/CrystalShaw"
path_to_dropbox <- "~/Dropbox/Projects"

#---- read in analytical sample ----
CESD_data_wide <- 
  read_csv(paste0(path_to_dropbox, 
                  "/exposure_trajectories/data/", 
                  "CESD_data_wide.csv"), 
           col_types = cols(.default = col_double(), HHIDPN = col_character(), 
                            death2018 = col_integer(), DOD = col_character(), 
                            Bday = col_character(), ed_cat = col_factor(), 
                            drop = col_logical(), r4mstat_cat = col_factor(), 
                            r9mstat_cat = col_factor(),
                            drinking4_cat_impute = col_factor(),
                            drinking9_cat_impute = col_factor(),
                            female = col_factor(), hispanic = col_factor(), 
                            black = col_factor(), other = col_factor(), 
                            smoker = col_integer()))

#---- select variables ----
CESD_vars <- c("HHIDPN", "4age_y", "9age_y", "female", "hispanic", "black", 
               "other", "smoker", "drinking4_cat_impute", 
               "drinking9_cat_impute", "r4mstat_cat", "r9mstat_cat", 
               "r4cesd_elevated", "r9cesd_elevated", "avg_cesd_elevated", 
               "total_elevated_cesd", "survtime", "observed")

CESD_subset <- CESD_data_wide %>% dplyr::select(all_of(CESD_vars))

#---- Check missingness ----
colSums(is.na(CESD_subset))

#---- induce missingness ----
CESD_only <- CESD_data_wide %>% 
  dplyr::select("HHIDPN", paste0("r", seq(4, 9), "cesd"))

mcar10 <- CESD_only
for(i in 2:ncol(mcar10)){
  mask <- sample(seq(1, nrow(mcar10)), size = ceiling(0.10*nrow(mcar10)), 
                 replace = FALSE)
  mcar10[mask, i] <- NA
}

# #Sanity check
# sum(is.na(mcar10))
# floor(0.10*nrow(mcar10)*(ncol(mcar10) - 1))
# table(rowSums(is.na(mcar10)))

#---- MI ----

#---- imputed exposure defs ----


