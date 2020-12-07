#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "lcmm")

#No scientific notation
options(scipen = 999)

set.seed(20200819)

#---- color palette ----
light_blue <- "#B4DAE5FF"

#---- Note ----
# Since the difference between win and OS, put substituted directory here
# Yingyan's directory: C:/Users/yingyan_wu
#                      C:/Users/yingyan_wu/Dropbox
# Crystal's directory: /Users/CrystalShaw
#                     ~/Dropbox/Projects

#Changing directories here will change them throughout the script
path_to_box <- "/Users/CrystalShaw"
path_to_dropbox <- "~/Dropbox/Projects"

#---- Read in analytical sample ----
CESD_data_wide <- 
  read_csv(paste0(path_to_dropbox, 
                  "/exposure_trajectories/data/", 
                  "hrs_samp_6CESD_waves4-9.csv"), 
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

#---- Cap age at baseline at 90 ----
# #Sanity check
# max(CESD_data_wide$`4age_y_int`)

CESD_data_wide %<>% filter(`4age_y_int` <= 90)

# #Sanity check
# hist(BMI_data_wide$`4age_y_int`)
# test <- BMI_data_wide %>% dplyr::select(paste0(seq(4, 9, by = 1), "BMI")) %>% 
#   is.na() %>% sum()

#---- Sample sizes ----
num_people = nrow(CESD_data_wide)
num_obs = sum(!is.na(CESD_data_wide %>% 
                       dplyr::select(paste0("r", seq(4, 9), "cesd"))))

#stratified by age at baseline
by_age_baseline <- data.frame("start" = seq(50, 85, by = 5)) %>%
  mutate("end" = start + 4, 
         "n" = 0)
by_age_baseline[nrow(by_age_baseline), "end"] <- 90

for(i in 1:nrow(by_age_baseline)){
  by_age_baseline[i, "n"] <- CESD_data_wide %>% 
    filter(`4age_y_int` %in% 
             seq(by_age_baseline[i, "start"], 
                 by_age_baseline[i, "end"], by = 1)) %>% nrow()
}

#stratified by overall age
overall_ages <- CESD_data_wide %>% 
  dplyr::select(paste0(seq(4, 9, by = 1), "age_y_int")) %>% 
  pivot_longer(everything())

by_age_overall <- data.frame("start" = seq(50, 95, by = 5)) %>%
  mutate("end" = start + 4, 
         "n" = 0)
by_age_overall[nrow(by_age_overall), "end"] <- max(overall_ages$value)

for(i in 1:nrow(by_age_overall)){
  by_age_overall[i, "n"] <- overall_ages %>% 
    filter(value %in% 
             seq(by_age_overall[i, "start"], 
                 by_age_overall[i, "end"], by = 1)) %>% nrow()
}

# #Sanity check
# sum(by_age_baseline$n)
# sum(by_age_overall$n)

#---- E1a Def: CESD at HRS wave 4 (1998) ----
#Effect of E1a on survival to HRS wave 14 (2018) 
CESD_data_wide %<>% 
  mutate("r4cesd_elevated" = ifelse(r4cesd > 4, 1, 0))

# #Sanity check
# table(CESD_data_wide$r4cesd, CESD_data_wide$r4cesd_elevated, useNA = "ifany")

#---- E1b Def: CESD at HRS wave 9 (2008) ----
#Effect of E1a on survival to HRS wave 14 (2018) 
CESD_data_wide %<>% 
  mutate("r9cesd_elevated" = ifelse(r9cesd > 4, 1, 0))

# #Sanity check
# table(CESD_data_wide$r9cesd, CESD_data_wide$r9cesd_elevated, useNA = "ifany")

#---- E2a Def: Cumulative Exposure (number occasions) ----
#Number of occasions with elevated depressive symptoms in HRS waves 4-9
elevated_cesd <- CESD_data_wide %>% 
  dplyr::select(paste0("r", seq(4, 9, by = 1), "cesd"))

elevated_cesd <- (elevated_cesd > 4)*1

CESD_data_wide %<>% mutate("total_elevated_cesd" = rowSums(elevated_cesd))

# #Sanity check
# head(elevated_cesd)
# head(CESD_data_wide$total_elevated_cesd)

#---- E2b Def: Cumulative Exposure (average CESD score) ----
CESD_data_wide %<>% 
  mutate("avg_cesd" = CESD_data_wide %>% 
           dplyr::select(paste0("r", seq(4, 9, by = 1), "cesd")) %>% 
           rowMeans(), 
         "avg_cesd_elevated" = ifelse(avg_cesd > 4, 1, 0))

# #Sanity check
# View(CESD_data_wide %<>% 
#        dplyr::select(paste0("r", seq(4, 9, by = 1), "cesd"), "avg_cesd", 
#                      "avg_cesd_elevated"))
# plot(CESD_data_wide$avg_cesd, CESD_data_wide$avg_cesd_elevated)

#---- E3 Def: Latent Classes

#---- Survival times from HRS wave 9 (2008) to HRS wave 14 (2018) ----
CESD_data_wide %<>% mutate(survtime = age_death_y - `9age_y_int`) %>% 
  mutate(survtime = ifelse(is.na(survtime), 10, survtime), 
         observed = ifelse(is.na(age_death_y), 0, 1))

# #Sanity check
# View(CESD_data_wide %>% dplyr::select("age_death_y", "survtime", "observed"))

#---- Figures ----
plot_data <- CESD_data_wide %>% 
  dplyr::select(paste0("r", seq(4, 9), "cesd")) %>% 
  pivot_longer(everything())

ggplot(data = plot_data, aes(value)) + 
  geom_bar(color = light_blue, fill = light_blue) + 
  theme_minimal() + xlab("CES-D Score")
ggsave(filename = "all_CESD.jpeg", device = "jpeg", width = 7, height = 5, 
       units = "in", 
       path = paste0(path_to_dropbox,
                     "/exposure_trajectories/manuscript/figures/"))

#---- save formatted dataset ----
write_csv(CESD_data_wide, paste0(path_to_dropbox,
                                   "/exposure_trajectories/data/",
                                   "CESD_data_wide.csv"))









