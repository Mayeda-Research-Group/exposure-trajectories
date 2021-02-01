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

#---- KAPPA dat_function----
KAPPA_dat_func <- function(occasion){
  
  # subset based on # of occasions
  CESD_data_subset <- CESD_data_wide %>%
    select(contains("elevated")) %>%
    filter(total_elevated_cesd == occasion)
  
  # vector to fix the positivity issue
  `0` <- c(0,0)
  
  # KAPPA test
  `4_9_table` <- 
    table(CESD_data_subset$r4cesd_elevated, CESD_data_subset$r9cesd_elevated)
  if(ncol( `4_9_table`) == 2){
    kappa_4_9 <- fmsb::Kappa.test(`4_9_table`)$ Result$estimate
  } else{
    kappa_4_9 <- NA
  }

  avg_4_table <- 
    table(CESD_data_subset$r4cesd_elevated, CESD_data_subset$avg_cesd_elevated)
  if(ncol(avg_4_table) == 2){
  kappa_4_avg <- fmsb::Kappa.test(avg_4_table)$ Result$estimate
  } else{
    kappa_4_avg <- NA
  }
  
  avg_9_table <- 
    table(CESD_data_subset$r9cesd_elevated, CESD_data_subset$avg_cesd_elevated)
  if(ncol(avg_9_table) == 2){
    kappa_9_avg <- fmsb::Kappa.test(avg_9_table)$ Result$estimate
  } else{
    kappa_9_avg <- NA
  }
  
  # KAPPA matrix
  kappa_value <- c(1, kappa_4_9,kappa_4_avg, 
                   kappa_4_9, 1, kappa_9_avg,
                   kappa_4_avg, kappa_9_avg, 1)
  rowname <- c("r4cesd_elevated", "r9cesd_elevated", "avg_cesd_elevated")
  colname <- c("r4cesd_elevated", "r9cesd_elevated", "avg_cesd_elevated")
  
  KAPPA_mat <- matrix(kappa_value, nrow = 3, ncol = 3, 
                      dimnames = list(rowname, colname))
  
  # melt to get the plotable dataset
  KAPPA_dat <- reshape2::melt(KAPPA_mat)
  return(KAPPA_dat)
  # heatmap(KAPPA_mat)
}

#---- Get the plotable datasets ----
# Overall
KAPPA_dat <- KAPPA_dat_func(CESD_data_wide$total_elevated_cesd)

# Stratified
KAPPA_dat_1 <- KAPPA_dat_func(1) %>% mutate("occasion" = 1)
KAPPA_dat_2 <- KAPPA_dat_func(2) %>% mutate("occasion" = 2)
KAPPA_dat_3 <- KAPPA_dat_func(3) %>% mutate("occasion" = 3)
KAPPA_dat_4 <- KAPPA_dat_func(4) %>% mutate("occasion" = 4)
KAPPA_dat_5 <- KAPPA_dat_func(5) %>% mutate("occasion" = 5)

KAPPA_dat_stratified <- 
  rbind(KAPPA_dat_1, KAPPA_dat_2, KAPPA_dat_3, KAPPA_dat_4, KAPPA_dat_5)

#---- heat map ----
# Overall
ggplot(data = KAPPA_dat,
       aes(x = Var1, y = Var2, fill = value)) +
  theme_minimal()+
  geom_tile(color = "white") + 
  scale_fill_gradient2(low = "blue", high = "red", 
                       limit = c(-1,1), name = "KAPPA\nStatistic")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 6)+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5)) + 
  labs(title = "KAPPA statistic heatmap")

# Stratified
ggplot(data = KAPPA_dat_stratified,
       aes(x = Var1, y = Var2, fill = value)) +
  theme_minimal() +
  facet_wrap(~ occasion, nrow = 2) +
  geom_tile(color = "white") + 
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       limit = c(-1,1), name = "KAPPA\nStatistic") +
  geom_text(
    aes(Var2, Var1, label = sprintf("%.3f", value)), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90),
    plot.title = element_text(hjust = 0.5)) + 
  labs(title = "KAPPA statistic heatmap stratified by # of elevated CESD")


#---- fix the problem for positivity issue
# CESD_data_subset <- CESD_data_wide %>%
#   select(contains("elevated")) %>%
#   filter(total_elevated_cesd == 5)
# 
# kappa_4_9 <- fmsb::Kappa.test(
#   table(CESD_data_subset$r4cesd_elevated, CESD_data_subset$r9cesd_elevated))$
#   Result$estimate
# kappa_4_avg <- fmsb::Kappa.test(
#   table(CESD_data_subset$r4cesd_elevated, CESD_data_subset$avg_cesd_elevated))$
#   Result$estimate
# 
# kappa_9_avg <- fmsb::Kappa.test(
#   table(CESD_data_subset$r9cesd_elevated, CESD_data_subset$avg_cesd_elevated))$
#   Result$estimate
# 
# a <- table(CESD_data_subset$r4cesd_elevated, CESD_data_subset$avg_cesd_elevated)
# 
# 
# #---- for occasion = 5 and 6 ----
# # Pending
# # Manually write it for now
# mat_func <- function(occasion){
#   CESD_data_subset <- CESD_data_wide %>%
#     select(contains("elevated")) %>%
#     filter(total_elevated_cesd == occasion)
#   
#   cross_tab_val <- c(sum(CESD_data_subset$r4cesd_elevated), )
#   mat_4_avg <- matrix(nrow = 2, nrow = 2)
#   
#   
#   ncol(a)
#   `0` <- c(0,0)
#   a <- cbind(a, `0`)
#   
#   CESD_data_subset <- CESD_data_wide %>%
#     select(contains("elevated")) %>%
#     filter(total_elevated_cesd == 6)
#   
#   
#   b <- matrix(c(0, 0, 0, 73), nrow = 2, ncol = 2)
#   fmsb::Kappa.test(b)
#   fmsb::Kappa.test(table(CESD_data_subset$r4cesd_elevated, CESD_data_subset$r9cesd_elevated))
#   
#   table(CESD_data_subset$r4cesd_elevated, CESD_data_subset$r9cesd_elevated)
