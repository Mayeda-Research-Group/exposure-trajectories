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
BMI_data_wide <- 
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
# max(BMI_data_wide$`4age_y_int`)

BMI_data_wide %<>% filter(`4age_y_int` <= 90)

# #Sanity check
# hist(BMI_data_wide$`4age_y_int`)
# test <- BMI_data_wide %>% dplyr::select(paste0(seq(4, 9, by = 1), "BMI")) %>% 
#   is.na() %>% sum()

#---- Sample sizes ----
num_people = nrow(BMI_data_wide)
num_obs = num_people*6

#stratified by age at baseline
by_age_baseline <- data.frame("start" = seq(50, 85, by = 5)) %>%
  mutate("end" = start + 4, 
         "n" = 0)
by_age_baseline[nrow(by_age_baseline), "end"] <- 90

for(i in 1:nrow(by_age_baseline)){
  by_age_baseline[i, "n"] <- BMI_data_wide %>% 
    filter(`4age_y_int` %in% 
             seq(by_age_baseline[i, "start"], 
                 by_age_baseline[i, "end"], by = 1)) %>% nrow()
}

#stratified by overall age
overall_ages <- BMI_data_wide %>% 
  dplyr::select(paste0(seq(4, 9, by = 1), "age_y_int")) %>% 
  pivot_longer(everything())

by_age_overall <- data.frame("start" = seq(50, 95, by = 5)) %>%
  mutate("end" = start + 4, 
         "n" = 0)
by_age_overall[nrow(by_age_overall), "end"] <- 100

for(i in 1:nrow(by_age_overall)){
  by_age_overall[i, "n"] <- overall_ages %>% 
    filter(value %in% 
             seq(by_age_overall[i, "start"], 
                 by_age_overall[i, "end"], by = 1)) %>% nrow()
}

# #Sanity check
# sum(by_age_baseline$n)
# sum(by_age_overall$n)

#---- E1 Def: Average BMI within 5-year age bands ----
#Effect of E1 on mortality within a decade of the end of follow-up

#---- **format dataset ----
for_dataset <- c("HHIDPN", "4age_y_int", "4BMI", "age_death_y")

E1_BMI_data_long <- BMI_data_wide %>% 
  dplyr::select(all_of(for_dataset)) %>% 
  set_colnames(c("HHIDPN", "age_BMI_y", "BMI", "age_death_y"))

for(wave in 5:9){
  for_dataset <- c("HHIDPN", paste0(wave, "age_y_int"), paste0(wave, "BMI"), 
                   "age_death_y")
  
  subset <- BMI_data_wide %>% dplyr::select(all_of(for_dataset)) %>% 
    set_colnames(c("HHIDPN", "age_BMI_y", "BMI", "age_death_y"))
  
  E1_BMI_data_long %<>% rbind(subset)
}

# #Sanity check
# dim(E1_BMI_data_long)
# View(E1_BMI_data_long)
# colSums(!is.na(E1_BMI_data_long))

E1_BMI_data_wide <- E1_BMI_data_long %>% 
  pivot_wider(names_from = age_BMI_y, values_from = BMI, 
              names_prefix = "BMI_")

# #Sanity check
# dim(E1_BMI_data_wide)
# colnames(E1_BMI_data_wide)
# colSums(!is.na(E1_BMI_data_wide))

#---- ****average BMI within age bands ----
for(i in seq(50, 95, by = 5)){
  if(i == 95){j = max(E1_BMI_data_long$age_BMI_y)} 
  else{j = i + 4}
  
  E1_BMI_data_wide[, paste0("BMI_", i, "-", j)] = 
    apply(E1_BMI_data_wide %>% 
            dplyr::select(paste0("BMI_", seq(i, j, by = 1))), 
          1, function(x) mean(x, na.rm = TRUE))
  
  E1_BMI_data_wide[, paste0("BMI_", i, "-", j, "_cat")] = 
    apply(E1_BMI_data_wide %>% 
            dplyr::select(paste0("BMI_", i, "-", j)), 
          1, function(x) case_when(x < 18.5 ~ "Underweight", 
                                   x >= 18.5 & x < 25 ~ "Normal", 
                                   x >= 25 & x < 30 ~ "Overweight", 
                                   x >= 30 ~ "Obese"))
  # #Sanity check
  # test <- E1_BMI_data_wide %>% 
  #   dplyr::select(paste0("BMI_", i, "-", j), 
  #                 paste0("BMI_", seq(i, j, by = 1)))
  # View(test)
  # colSums(!is.na(test))
}

# #Sanity check
# View(E1_BMI_data_wide %>% dplyr::select(contains("BMI_")))

#Get rid of columns for individual ages
E1_BMI_data_wide %<>% 
  dplyr::select(-paste0("BMI_", 
                        seq(min(E1_BMI_data_long$age_BMI_y), 
                            max(E1_BMI_data_long$age_BMI_y), by = 1)))

#---- ****death indicators ----
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








