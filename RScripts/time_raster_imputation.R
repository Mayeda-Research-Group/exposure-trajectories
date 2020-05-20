#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here")

set.seed(20200520)

#---- Read in data ----
analytic_df <- read_csv(here::here("Data", "analytic_df.csv"))

