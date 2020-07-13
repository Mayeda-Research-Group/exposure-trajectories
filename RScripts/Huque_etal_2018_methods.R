#---- Package loading + seed setting ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here")

set.seed(20200517)
