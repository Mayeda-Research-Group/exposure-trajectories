if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("table1", "tidyverse", "plyr", "labelled")

#---- Note ----
# Since the difference between win and OS, put substituted directory here
# Yingyan's directory: C:/Users/yingyan_wu
#                      C:/Users/yingyan_wu/Dropbox
# Crystal's directory: /Users/CrystalShaw
#                     ~/Dropbox/Projects

#Changing directories here will change them throughout the script
path_to_box <- "C:/Users/yingyan_wu"
path_to_dropbox <- "C:/Users/yingyan_wu/Dropbox"

#---- read in analytical sample ----
CESD_data <- 
  read_csv(paste0(path_to_dropbox, 
                  "/exposure_trajectories/data/", 
                  "hrs_samp_6CESD_waves4-9.csv"), 
           col_types = cols(.default = col_double(), HHIDPN = col_character(), 
                            death2018 = col_integer(), ed_cat = col_factor(), 
                            r4mstat_cat = col_factor(), 
                            r9mstat_cat = col_factor(), 
                            drinking4_cat_impute = col_factor(),
                            drinking9_cat_impute = col_factor(),
                            female = col_factor(), hispanic = col_factor(), 
                            black = col_factor(), other = col_factor(), 
                            smoker = col_integer()))
#---- Label the data ----
# colnames(CESD_data)

CESD_data %<>% 
  mutate("raceeth" = case_when(hispanic == 0 & white == 1 ~ "Non-hispanic White",
                               hispanic == 0 & black == 1 ~ "Non-hispanic Black",
                               hispanic == 1 ~ "Hispanic",
                               other == 1 ~ "Other")) %>%
  set_variable_labels(
    r4age_y_int = "Baseline age in years",
    female = "Sex/Gender",
    raceeth = "Race/Ethnicity",
    ed_cat = "Education level",
    r4mstat_cat = "Married status",
    r4hibpe_impute = "Hypertension(Diagnosed)",
    r4diabe_impute = "Diabetes",
    r4hearte_impute = "Heart diseases",
    r4stroke_impute = "Stroke",
    r4cancre_impute = "Cancer",
    r4lunge_impute = "Lung diseases",
    r4memrye_impute = "Memory problems",
    r4conde_impute = "Count of self-report chronic diseases",
    r4BMI = "BMI",
    drinking4_cat_impute = "Alcohol intake",
    smoker = "Smoking status",
    r4shlt = "Self-reported health",
    r4cesd = "CESD score"
  ) %>%
  drop_unused_value_labels() %>%
  set_value_labels(female = c("Female" = 1, "Male" = 0),
                   r4hibpe_impute = c("Yes" = 1, "No" = 0),
                   r4diabe_impute = c("Yes" = 1, "No" = 0),
                   r4hearte_impute = c("Yes" = 1, "No" = 0),
                   r4stroke_impute = c("Yes" = 1, "No" = 0),
                   r4cancre_impute = c("Yes" = 1, "No" = 0),
                   r4lunge_impute = c("Yes" = 1, "No" = 0),
                   r4memrye_impute = c("Yes" = 1, "No" = 0),
                   drinking4_cat_impute = c("Heavy Drinking" = 3,
                                            "Moderate Drinking" = 2,
                                            "No Drinking" = 1),
                   smoker = c("Ever smoke" = 1, "No smoking" = 0),
                   r4shlt = c("Excellent" = 1, "Very Good" = 2, "Good" = 3,
                              "Fair" = 4, "Poor" = 5)) %>%
  modify_if(is.labelled, to_factor)

table1(~ r4age_y_int + female + raceeth + 
         ed_cat + r4mstat_cat + 
         r4hibpe_impute + r4diabe_impute + r4hearte_impute + r4stroke_impute + 
         r4cancre_impute + r4lunge_impute + r4memrye_impute +
         r4conde_impute + 
         r4BMI + drinking4_cat_impute + smoker + r4shlt|r4cesd, data=CESD_data)

