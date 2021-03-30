if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("table1", "tidyverse", "dplyr", "plyr", "labelled", "gt", "rvest")

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
CESD_data_wide <- 
  read_csv(paste0(path_to_dropbox, 
                  "/exposure_trajectories/data/", 
                  "CESD_data_wide.csv"), 
           col_types = cols(.default = col_double()))
           
#---- Label the data ----
colnames(CESD_data_wide)

table1_data <- CESD_data_wide %>% 
  mutate("raceeth" = case_when(hispanic == 0 & white == 1 ~ "Non-hispanic White",
                               hispanic == 0 & black == 1 ~ "Non-hispanic Black",
                               hispanic == 1 ~ "Hispanic",
                               other == 1 ~ "Other"),
         "r4mstat_cat" = case_when(
           r4married_partnered == 1 ~ "Married/Partnered",
           r4not_married_partnered == 1 ~ "Not Married/Partnered",
           r4widowed == 1 ~ "Widowed")) %>%
  labelled::set_variable_labels(
    r4age_y_int = "Baseline age in years",
    female = "Sex/Gender (%)",
    raceeth = "Race/Ethnicity (%)",
    ed_cat = "Education level (%)",
    r4mstat_cat = "Married status (%)",
    r4hibpe_impute = "Hypertension(Diagnosed) (%)",
    r4diabe_impute = "Diabetes (%)",
    r4hearte_impute = "Heart diseases (%)",
    r4stroke_impute = "Stroke (%)",
    r4cancre_impute = "Cancer (%)",
    r4lunge_impute = "Lung diseases (%)",
    r4memrye_impute = "Memory problems (%)",
    r4conde_impute = "Count of self-report chronic diseases",
    r4BMI = "BMI",
    r4drinking_cat = "Alcohol intake (%)",
    smoker = "Smoking status (%)",
    r4shlt = "Self-reported health (%)",
    r4cesd_elevated = "Baseline Elevated CES-D"
  ) %>%
  labelled::drop_unused_value_labels() %>%
  labelled::set_value_labels(female = c("Female" = 1, "Male" = 0),
                   ed_cat = c("Less than High School" = 1, "High School" = 2, 
                              "Some college" = 3, "Bachelors" = 4, 
                              "Grad studies" = 5),
                   r4hibpe_impute = c("Yes" = 1, "No" = 0),
                   r4diabe_impute = c("Yes" = 1, "No" = 0),
                   r4hearte_impute = c("Yes" = 1, "No" = 0),
                   r4stroke_impute = c("Yes" = 1, "No" = 0),
                   r4cancre_impute = c("Yes" = 1, "No" = 0),
                   r4lunge_impute = c("Yes" = 1, "No" = 0),
                   r4memrye_impute = c("Yes" = 1, "No" = 0),
                   r4drinking_cat = c("Heavy Drinking" = 2,
                                            "Moderate Drinking" = 1,
                                            "No Drinking" = 0),
                   smoker = c("Ever smoke" = 1, "No smoking" = 0),
                   r4shlt = c("Excellent" = 1, "Very Good" = 2, "Good" = 3,
                              "Fair" = 4, "Poor" = 5),
                   r4cesd_elevated = c("Elevated CES-D" = 1, 
                                       "Not Elevated CES-D" = 0)) %>%
  modify_if(is.labelled, to_factor)

#---- Use table1 pacakge ----
# For showing +/- SD in the table (save it here in case future needed)
# render_cont <- function(x){
#   with(stats.apply.rounding(stats.default(x), digits = 2), 
#        c("", "Mean (SD)" = sprintf("%s (&plusmn; %s)", MEAN, SD)))
# }

# Just showing mean(sd) for continuous variables
render_cont <- function(x){
  with(stats.apply.rounding(stats.default(x), digits = 2),
       c("", "Mean (SD)" = sprintf("%s (%s)", MEAN, SD)))
}
# just showing numbers(percentage) without "%" for categorical variables
render_cat <- function(x){
  c("", sapply(stats.default(x), 
               function(y) with(y, sprintf("%d (%.1f)", FREQ, PCT))))
}

(table_1 <- table1::table1(~ r4age_y_int + female + raceeth + 
                     ed_cat + r4mstat_cat + 
                     r4hibpe_impute + r4diabe_impute + r4hearte_impute + 
                     r4stroke_impute + r4cancre_impute + r4lunge_impute + 
                     r4memrye_impute +
                     r4conde_impute + 
                     r4BMI + r4drinking_cat + smoker + r4shlt|r4cesd_elevated, 
                   render.continuous = render_cont,
                   render.categorical = render_cat,
                   footnote = 
                     "Stratified by Elevated baseline CES-D (CES-D > 4)",
                   data = table1_data))

table_1_df <- as.data.frame(read_html(table_1) %>% html_table(fill=TRUE))

readr::write_csv(table_1_df, path = paste0(path_to_dropbox, 
                                    "/exposure_trajectories/manuscript/tables", 
                                    "/table_1_df.csv"))

# Attempt to save the html using gtsave

# table_1_gt_tbl <- as_tibble(read_html(table_1) %>% html_table(fill=TRUE), 
#                             .name_repair = c())
# 
# gt(table_1_gt_tbl)
# 
# gt(
#   data,
#   rowname_col = "rowname",
#   groupname_col = dplyr::group_vars(data),
#   rownames_to_stub = FALSE,
#   auto_align = TRUE,
#   id = NULL,
#   row_group.sep = getOption("gt.row_group.sep", " - ")
# )
# 
# gtsave(, "tab_1.html", inline_css = TRUE,
#        path = paste0(path_to_dropbox, "/manuscript/tables")
# )

#---- Use furniture package ----
# pacman::p_load("haven", "tidyverse", "magrittr", "foreign", "ggplot2", "dplyr", 
#        "survey", "tidyr", "lme4", "lmerTest", "leaps", "DescTools", "locfit", 
#        "jtools", "LMERConvenienceFunctions", "sjPlot", "sjmisc", "sjlabelled", 
#        "ggplot2", "furniture")
# library(dplyr)
# library(furniture)
# # p_load("furniture")
# # detach("package:plyr", unload = TRUE)
# # detach("package:table1", unload = TRUE)
# table_1 <- table1_data %>% 
#   furniture::table1(r4age_y_int, female, 
#                             #  raceeth, ed_cat, r4mstat_cat, 
#                             # r4hibpe_impute, r4diabe_impute, r4hearte_impute, 
#                             # r4stroke_impute, r4cancre_impute, r4lunge_impute, 
#                             # r4memrye_impute, r4conde_impute, r4BMI, 
#                             # r4drinking_cat, smoker, r4shlt,
#                             # splitby = ~r4cesd_elevated, 
#                             # na.rm = FALSE,
#                             # formate_number = TRUE,
#                             type = c("condensed"))
