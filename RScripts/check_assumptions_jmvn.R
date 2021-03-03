#This script checks the normality of residuals assumption necessary for FCS with 
#the "norm" option to be equivalent to Joint Multivariate Normal Modeling
#
#Transformations will be made if they help make the assumption more plausible

#---- package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here")

#No scientific notation
options(scipen = 999)

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
           col_types = cols(HHIDPN = col_character())) %>% 
  mutate_if(is.character, as.factor)

#---- select variables ----
vars <- c(paste0("r", seq(4, 9), "mstat_impute"), "ed_cat", 
          paste0("drinking", seq(4, 9), "_impute"), 
          paste0("r", seq(4, 9), "memrye_impute"),
          paste0("r", seq(4, 9), "stroke_impute"),
          paste0("r", seq(4, 9), "hearte_impute"),
          paste0("r", seq(4, 9), "lunge_impute"),
          paste0("r", seq(4, 9), "cancre_impute"),
          paste0("r", seq(4, 9), "hibpe_impute"),
          paste0("r", seq(4, 9), "diabe_impute"), "smoker", 
          paste0("r", seq(4, 9), "BMI"), "hispanic", "white", "black", "other", 
          "female", paste0("r", seq(4, 9), "age_y_int"), "death2018", 
          paste0("r", seq(3, 9), "cesd"), paste0("r", seq(4, 9), "shlt"))

model_subset <- CESD_data_wide %>% dplyr::select(all_of(vars))

#---- marital status ----
for(i in (4:9)){
  outcome <- paste0("r", i, "mstat_impute")
  model <- 
    lm(as.formula(paste0(outcome, "~ .")), data = model_subset)
  show(hist(model$residuals, main = outcome))
}

#---- drinking status ----
for(i in (4:9)){
  outcome <- paste0("drinking", i, "_impute")
  model <- 
    lm(as.formula(paste0(outcome, "~ .")), data = model_subset)
  show(hist(model$residuals, main = outcome))
}

#---- chronic conditions ----
#Waves 4 and 9 tend to deviate the most from normality, but doesn't seem like 
# it's worth it to do a transformation
conditions <- c("memrye", "stroke", "hearte", "lunge", "cancre", "hibpe", 
                "diabe")

for(condition in conditions){
  for(i in (4:9)){
    outcome <- paste0("r", i, condition, "_impute")
    model <- 
      lm(as.formula(paste0(outcome, "~ .")), data = model_subset)
    show(hist(model$residuals, main = outcome))
  }
}

#---- BMI ----
for(i in (4:9)){
  outcome <- paste0("r", i, "BMI")
  model <- 
    lm(as.formula(paste0(outcome, "~ .")), data = model_subset)
  show(hist(model$residuals, main = outcome))
}

#---- CESD ----
for(i in (4:9)){
  outcome <- paste0("r", i, "cesd")
  model <- 
    lm(as.formula(paste0(outcome, "~ .")), data = model_subset)
  show(hist(model$residuals, main = outcome))
}

#---- self-reported health ----
for(i in (4:9)){
  outcome <- paste0("r", i, "shlt")
  model <- 
    lm(as.formula(paste0(outcome, "~ .")), data = model_subset)
  show(hist(model$residuals, main = outcome))
}




