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
table(CESD_data_wide$total_elevated_cesd)

#---- KAPPA test ----

kappa_4_9 <- fmsb::Kappa.test(
  table(CESD_data_wide$r4cesd_elevated, CESD_data_wide$r9cesd_elevated))$
  Result$estimate

kappa_4_avg <- fmsb::Kappa.test(
  table(CESD_data_wide$r4cesd_elevated, CESD_data_wide$avg_cesd_elevated))$
  Result$estimate

kappa_9_avg <- fmsb::Kappa.test(
  table(CESD_data_wide$r9cesd_elevated, CESD_data_wide$avg_cesd_elevated))$
  Result$estimate

#---- KAPPA matrix ----
kappa_value <- c(NA, kappa_4_9,kappa_4_avg, 
                kappa_4_9, NA, kappa_9_avg,
                kappa_4_avg, kappa_9_avg, NA)
rowname <- c("r4cesd_elevated", "r9cesd_elevated", "avg_cesd_elevated")
colname <- c("r4cesd_elevated", "r9cesd_elevated", "avg_cesd_elevated")

KAPPA_mat <- matrix(kappa_value, nrow = 3, ncol = 3, dimnames = list(rowname, colname))
# heatmap(KAPPA_mat)

KAPPA_dat <- reshape2::melt(KAPPA_mat)

#---- heat map ----
ggplot(data = KAPPA_dat,
       aes(x = Var1, y = Var2, fill = value)) +
  theme_minimal()+
  geom_tile(color = "white") + 
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.6, limit = c(0,1))+
  geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 6)+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank())


