#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "broom")

#No scientific notation
options(scipen = 999)

set.seed(20200819)

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
                            drop = col_logical(), r9mstat_cat = col_factor(),
                            drinking9_cat = col_factor(),
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

#---- E1 Def: Median CESD within 5-year age bands ----
#Effect of E1 on mortality within a decade of the end of follow-up

#---- **format dataset ----
for_dataset <- c("HHIDPN", "4age_y_int", "r4cesd", "age_death_y")

E1_CESD_data_long <- CESD_data_wide %>% 
  dplyr::select(all_of(for_dataset)) %>% 
  set_colnames(c("HHIDPN", "age_CESD_y", "CESD", "age_death_y"))

for(wave in 5:9){
  for_dataset <- c("HHIDPN", paste0(wave, "age_y_int"), 
                   paste0("r", wave, "cesd"), 
                   "age_death_y")
  
  subset <- CESD_data_wide %>% dplyr::select(all_of(for_dataset)) %>% 
    set_colnames(c("HHIDPN", "age_CESD_y", "CESD", "age_death_y"))
  
  E1_CESD_data_long %<>% rbind(subset)
}

# #Sanity check
# dim(E1_CESD_data_long)
# View(E1_CESD_data_long)
# colSums(!is.na(E1_CESD_data_long))

E1_CESD_data_wide <- E1_CESD_data_long %>% 
  pivot_wider(names_from = age_CESD_y, values_from = CESD, 
              names_prefix = "CESD_")

# #Sanity check
# dim(E1_CESD_data_wide)
# colnames(E1_CESD_data_wide)
# colSums(!is.na(E1_CESD_data_wide))

#---- ****average CESD within age bands ----
for(i in seq(50, 95, by = 5)){
  if(i == 95){j = max(E1_CESD_data_long$age_CESD_y)} 
  else{j = i + 4}
  
  E1_CESD_data_wide[, paste0("CESD_", i, "-", j)] = 
    apply(E1_CESD_data_wide %>% 
            dplyr::select(paste0("CESD_", seq(i, j, by = 1))), 
          1, function(x) mean(x, na.rm = TRUE))
  
  E1_CESD_data_wide[, paste0("CESD_", i, "-", j, "_elev_dep_sx")] = 
    apply(E1_CESD_data_wide %>% 
            dplyr::select(paste0("CESD_", i, "-", j)), 
          1, function(x) case_when(x < 4 ~ 0, 
                                   x >= 4 ~ 1))
  # #Sanity check
  # test <- E1_CESD_data_wide %>% 
  #   dplyr::select(paste0("CESD_", i, "-", j), 
  #                 paste0("CESD_", seq(i, j, by = 1)))
  # View(test)
  # colSums(!is.na(test))
}

# #Sanity check
# View(E1_CESD_data_wide %>% dplyr::select(contains("CESD_")))

#Get rid of columns for individual ages
E1_CESD_data_wide %<>% 
  dplyr::select(-paste0("CESD_", 
                        seq(min(E1_CESD_data_long$age_CESD_y), 
                            max(E1_CESD_data_long$age_CESD_y), by = 1)))

# #Count people
# colSums(!is.na(E1_CESD_data_wide %>% dplyr::select(contains("CESD_"))))

#---- ****death indicator ----
#Choose the age bands
# #Sanity check
# max(E1_BMI_data_wide$age_death_y, na.rm = TRUE)
# min(E1_BMI_data_wide$age_death_y, na.rm = TRUE)

#No analyses in [50, 54] because we need this as the first "previous wave"
#Thus death ages start at [65, 69]
for(i in seq(65, 105, by = 5)){
  if(i == 105){j = max(E1_BMI_data_long$age_death_y, na.rm = TRUE)} 
  else{j = i + 4}
  
  E1_BMI_data_wide[, paste0("death_by_", i, "-", j)] = 
    apply(E1_BMI_data_wide %>% 
            dplyr::select("age_death_y"), 
          1, function(x) ifelse(x > j | is.na(x), 0, 1))
}

# #Sanity check
# View(E1_BMI_data_wide %>% dplyr::select("age_death_y", contains("death_by_")))
# colSums(!is.na(E1_BMI_data_wide))

#---- ****time-invariant covariates ----
time_invariant <- c("HHIDPN", "female", "white", "black", "hispanic", "other", 
                    "smoker")

E1_BMI_data_wide %<>% 
  left_join(., BMI_data_wide %>% dplyr::select(all_of(time_invariant)), 
            by = "HHIDPN")

# #Sanity check
# dim(E1_BMI_data_wide)
# colnames(E1_BMI_data_wide)

#---- ****save formatted dataset ----
write_csv(E1_BMI_data_wide, paste0(path_to_dropbox,
                                   "/exposure_trajectories/data/",
                                   "E1_BMI_data_wide.csv"))








