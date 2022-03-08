# 1. Data construction for the exposure: wave 4 - 9 CESD
# Data input: trk2018tr_r.sas7bdat, randhrs1992_2018v1.dta, cses_measures.dta
# Data output: CESD_data_wide.csv, 
#              CESD_missing_model_betas_in_sampled.csv (in commented code)
#              predict_death2018_betas.csv (in commented code)
# Author: CS & YW


#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}
p_load("here", "readr", "tidyverse", "magrittr", "plyr", "haven", "labelled",
       "lubridate")

#No scientific notation
options(scipen = 999)

#---- Note ----
# Since the difference between win and OS, put substituted directory here
# Yingyan's directory: C:/Users/Yingyan Wu/Box
#                      C:/Users/Yingyan Wu/Dropbox
# Crystal's directory: /Users/crystalshaw/Library/CloudStorage/Box-Box
#                     ~/Dropbox/Projects

#Changing directories here will change them throughout the script
path_to_box <- "C:/Users/Yingyan Wu/Box"
path_to_dropbox <- "C:/Users/Yingyan Wu/Dropbox"

#---- source scripts ----
source(here::here("RScripts", "functions", "impute_ages.R"))
source(here::here("RScripts", "functions", "measured_self_report.R"))
source(here::here("Rscripts", "functions", "impute_condition.R"))

#---- wave mapping between HRS and RAND ----
#Wave Year | HRS Core Data | RAND
# 1992 | V | 1
# 1994 | W | 2
# 1996 | E | 3
# 1998 | F | 4
# 2000 | G | 5
# 2002 | H | 6
# 2004 | J | 7
# 2006 | K | 8
# 2008 | L | 9
# 2010 | M | 10
# 2012 | N | 11
# 2014 | O | 12
# 2016 | P | 13

number_waves <- seq(1, 13, by = 1) 

#---- read in HRS tracker ----
hrs_tracker <-
  read_sas(paste0(path_to_box, "/HRS/tracker/trk2018v2a/",
                  "trk2018tr_r.sas7bdat")) %>%
  select("HHID", "PN", "PIWTYPE", "PALIVE", "QIWTYPE", "QALIVE", 
         paste0(c("F", "G", "H", "J", "K", "L"), "ALIVE")) %>%
  unite("HHIDPN", c("HHID", "PN"), sep = "", remove = TRUE) %>%
  mutate_at("HHIDPN", as.numeric)

#---- read in RAND file ----
#reading in STATA file because SAS file wouldn't load
#Variables of interest: 
## Demographics: HHIDPN, gender, race, hispanic, birth month, 
##               birth year, birth date, death month, death year, death date,
##               age data in months (biomarker waves), years of education, 
##               highest degree (masked),  
## Health: weight (kg; measured; waves 8+),
##         weight (kg; self-report),
##         height (m; measured; waves 8+),
##         height (m; self-report),
##         BMI (measured; waves 8+),
##         BMI (self-report), 
##         waist circumference,
##         BP (systolic; waves 8+), 
##         BP (diastolic; waves 8+), 
##         reports high blood pressure this wave,
##         reports diabetes this wave,
##         diabetes ever/never,
##         reports cancer this wave,
##         reports stroke this wave,
##         reports heart problems this wave，
##         report memory prob this wv(4-9), 
##         report lung disease this wv,
##         No(report arthritis since last wv (2-13),)
##         count of chronic conditions
##         CESD (depression; from wave 2-13)
##         Self-reported health
## Health Behaviors: current smoker 
##                   number of days drinking per week (waves 3+)
##                   number of drinks per day (waves 3+)
##                   frequency of vigr/modr/light physical activity (waves 7+)
## 
# Note: Dates are formatted as SAS dates (days from January 1, 1960)

rand_variables <- c("hhidpn", "ragender", "raracem", "rahispan", "rabmonth", 
                    "rabyear", "rabdate", "radmonth", "radyear", "raddate",
                    paste0("r", number_waves, "agem_e"), "raedyrs", 
                    "raedegrm", 
                    paste0("r", number_waves, "mstat"),
                    paste0("r", seq(8, 13, by = 1), "pmwght"), 
                    paste0("r", number_waves, "weight"),
                    paste0("r", seq(8, 13, by = 1), "pmhght"), 
                    paste0("r", number_waves, "height"),
                    paste0("r", number_waves, "bmi"), 
                    paste0("r", seq(8, 13, by = 1), "pmbmi"),
                    paste0("r", seq(8, 13, by = 1), "pmwaist"),
                    paste0("r", seq(8, 13, by = 1), "bpsys"), 
                    paste0("r", seq(8, 13, by = 1), "bpdia"),
                    paste0("r", number_waves, "hibpe"),
                    paste0("r", number_waves, "diabe"),
                    paste0("r", number_waves, "cancre"),
                    paste0("r", number_waves, "stroke"),
                    paste0("r", number_waves, "hearte"),
                    paste0("r", number_waves, "lunge"),
                    paste0("r", seq(4 ,9 , by = 1), "memrye"),
                    paste0("r", number_waves, "conde"),
                    paste0("r", number_waves, "smoken"), 
                    paste0("r", seq(3, 13, by = 1), "drinkd"),
                    paste0("r", seq(3, 13, by = 1), "drinkn"),
                    paste0("r", seq(2, 13, by = 1), "cesd"),
                    paste0("r", number_waves, "shlt"))

RAND <- read_dta(paste0(path_to_box, "/HRS/RAND_longitudinal/STATA/", 
                        "randhrs1992_2018v1.dta"), 
                 col_select = all_of(rand_variables)) 

colnames(RAND)[1] <- "HHIDPN" #For merging

#Remove labeled data format
val_labels(RAND) <- NULL

#---- read in Anusha Vable's CSES index ----
cSES <- read_dta(paste0(path_to_dropbox, "/exposure_trajectories/data/", 
                        "cSES measures/cses_measures.dta"), 
                 col_select = c("hhid", "pn", "cses_index")) %>% 
  unite(col = "HHIDPN", c("hhid", "pn"), sep = "") %>% 
  mutate_at("HHIDPN", as.numeric)

#---- merge datasets ----
#Use this to subset RAND data
hrs_samp <- 
  join_all(c(list(RAND, cSES, hrs_tracker)), by = "HHIDPN", type = "left") 

# Create indicators for whether CESD is missing at each wave
hrs_samp[, paste0("r", seq(2, 13), "cesd_missing")] <- NA

cesd_mat <- hrs_samp %>% select(contains("cesd"))
for (j in 1:length(seq(2, 13))) {
  cesd_mat[, j + length(seq(2, 13))] <- ifelse(is.na(cesd_mat[, j]), 1, 0)
}

hrs_samp[, colnames(cesd_mat)] <- cesd_mat

# CESD missingness proportion
cesd_missing <- cesd_mat %>%
  select(contains("missing")) %>%
  colSums(is.na(.))
cesd_missing/nrow(hrs_samp) 
# Sanity check
# table(cesd_mat$r2cesd, cesd_mat$r2cesd_missing, useNA = "ifany")
# table(cesd_mat$r13cesd, cesd_mat$r13cesd_missing, useNA = "ifany")
# table(cesd_mat$r6cesd, cesd_mat$r6cesd_missing, useNA = "ifany")

#---- death ----
hrs_samp %<>% mutate("death2018" = ifelse(is.na(raddate), 0, 1))

#format RAND dates with lubridate
hrs_samp %<>% mutate("DOD" = as.Date(hrs_samp$raddate, origin = "1960-01-01"), 
                     "Bday" = as.Date(hrs_samp$rabdate, origin = "1960-01-01"))

# #Sanity check
# View(hrs_samp[, c("Bday", "rabmonth", "rabyear")] %>% na.omit())
# View(hrs_samp[, c("DOD", "radmonth", "radyear")] %>% na.omit())

#age at death
hrs_samp %<>% 
  mutate("age_death_d" = difftime(DOD, Bday, units = "days"), 
         "age_death_y" = floor(as.numeric(age_death_d/365.25)))

# #Sanity check
# View(hrs_samp[, c("Bday", "DOD", "age_death_y")] %>% na.omit())
# View(hrs_samp[, c("age_death_y", "death2016", "death2018")] %>% 
#        filter(death2018 == 1))

#Drop RAND birth date, death date variables, and derived death in days
hrs_samp %<>% dplyr::select(-c("rabmonth", "rabyear", "rabdate", 
                               "radmonth", "radyear", "raddate", 
                               "age_death_d"))

#---- age ----
age_m <- hrs_samp %>% dplyr::select(contains("agem_e")) %>% 
  apply(., 1, impute_ages) %>% t() 

#Exact ages
hrs_samp[, paste0("r", number_waves, "age_y")] <- age_m/12
#Ages rounded down to nearest year
hrs_samp[, paste0("r", number_waves, "age_y_int")] <- floor(age_m/12)

# #Sanity check
# View(hrs_samp[, c(paste0(number_waves, "age_y"), 
#                   paste0(number_waves, "age_y_int"))])

#Drop RAND age variables
hrs_samp %<>% dplyr::select(-paste0("r", number_waves, "agem_e"))

#---- gender ----
hrs_samp %<>% 
  mutate("female" = ifelse(ragender == 2, 1, 0))

# #sanity check
# table(hrs_samp$female, hrs_samp$ragender, useNA = "ifany")

#Drop RAND gender variable 
hrs_samp %<>% dplyr::select(-ragender)

#---- race-eth ----
#Code any hispanic as 1, else 0
hrs_samp %<>% 
  mutate("hispanic" = ifelse(rahispan == 0 | is.na(rahispan), 0, 1)) %>% 
  mutate("white" = ifelse(raracem == 1 & hispanic == 0, 1, 0)) %>%
  mutate("black" = ifelse(raracem == 2 & hispanic == 0, 1, 0)) %>% 
  mutate("other" = ifelse(raracem == 3 & hispanic == 0, 1, 0)) %>% 
  mutate("unknown_race_eth" = ifelse(is.na(raracem) & hispanic == 0, 1, 0))

# #Sanity check
# table(hrs_samp$hispanic, hrs_samp$rahispan, useNA = "ifany")
# table(hrs_samp$hispanic, hrs_samp$raracem, hrs_samp$black, useNA = "ifany")
# table(hrs_samp$hispanic, hrs_samp$raracem, hrs_samp$other, useNA = "ifany")
# table(hrs_samp$hispanic, hrs_samp$raracem, hrs_samp$unknown_race_eth,
#       useNA = "ifany")
# table(hrs_samp$unknown_race_eth, useNA = "ifany")
# colSums(is.na(hrs_samp[, c("hispanic", "white", "black", "other", 
#                            "unknown_race_eth")]))
# table(hrs_samp$white, hrs_samp$raracem, useNA = "ifany")
# table(hrs_samp$rahispan, hrs_samp$hispanic, useNA = "ifany")

hrs_samp %<>% 
  #Drop the RAND rahispan variable (recoded as hispanic) and race variables
  dplyr::select(-c("rahispan", "raracem"))

#---- education ----
#Create education categories
hrs_samp %<>% 
  mutate("ed_cat" = 
           case_when(raedyrs < 12 ~ 1, 
                     raedyrs == 12 ~ 2, 
                     raedyrs > 12 & raedyrs < 16 ~ 3, 
                     raedyrs == 16 ~ 4, 
                     raedyrs > 16 ~ 5)) 

# #Sanity check
# table(is.na(hrs_samp$raedyrs))
# table(hrs_samp$raedyrs)
# table(hrs_samp$raedyrs, hrs_samp$ed_cat, useNA = "ifany")

#---- cSES index ----
# #Sanity check
# table(is.na(hrs_samp$cses_index))
# 
# # none missing!

#---- height ----
#Create a "best" height variable by taking the median of measured heights 
# (waves 8+) if available or first self-reported height
hrs_samp %<>% 
  mutate("med_height" = hrs_samp %>% 
           dplyr::select(paste0("r", seq(8, 13, by = 1), "pmhght")) %>%
           apply(1, function(x) median(x, na.rm = TRUE)), 
         "self_height" = hrs_samp %>%
           dplyr::select(paste0("r", number_waves, "height")) %>%
           apply(1, function(x) 
           {if (sum(is.na(x)) == length(x)){
             return(NA)} else {x[min(which(!is.na(x)))]}})) %>%
  mutate("height_measured" = ifelse(!is.na(med_height), 1, 0)) %>%
  mutate("height" = ifelse(height_measured == 1, med_height, self_height))

# #Sanity check
# View(hrs_samp[, c(paste0("r", seq(8, 13, by = 1), "pmhght"),
#                   paste0("r", number_waves, "height"),
#                   "med_height", "self_height", "height_measured", "height")] %>%
#        filter(is.na(height)))

#---- weight ----
hrs_samp %<>% 
  cbind(measured_self_report(data = hrs_samp, 
                             measured_cols = 
                               paste0("r", seq(8, 13, by = 1), "pmwght"), 
                             self_cols = 
                               paste0("r", number_waves, "weight"), 
                             derived_variable = "weight",
                             measured_waves_start = 8, all_waves_end = 13))

# #Checking weird weight values (from YW's analysis)
# weird_values <- c(12738020, 16381010, 31605040, 46769010, 47242010, 73050020,
#                   75888040, 76793010, 79557010, 112747020, 146359010,
#                   207810010, 208021020, 208024010, 210114010, 35201010, 
#                   25391010)

# View(hrs_samp %>% filter(HHIDPN %in% weird_values) %>%
#        dplyr::select("HHIDPN", c(paste0(seq(4, 9, by = 1), "weight"))))

#Set the weird measures to NA
hrs_samp[which(hrs_samp$HHIDPN %in% 
                 c(12738020, 16381010, 46769010, 73050020, 
                   146359010, 208024010)), "4weight"] <- NA
hrs_samp[which(hrs_samp$HHIDPN == 75888040), "7weight"] <- NA
hrs_samp[which(hrs_samp$HHIDPN == 47242010), "8weight"] <- NA
hrs_samp[which(hrs_samp$HHIDPN %in% 
                 c(31605040, 76793010, 79557010, 112747020, 
                   207810010, 208021020, 210114010, 35201010, 
                   25391010)), "9weight"] <- NA

#Fix measured indicators
hrs_samp[which(hrs_samp$HHIDPN %in% 
                 c(12738020, 16381010, 46769010, 73050020, 
                   146359010, 208024010)), "4weight_measured"] <- 0
hrs_samp[which(hrs_samp$HHIDPN == 75888040), "7weight_measured"] <- 0
hrs_samp[which(hrs_samp$HHIDPN == 47242010), "8weight_measured"] <- 0
hrs_samp[which(hrs_samp$HHIDPN %in% 
                 c(31605040, 76793010, 79557010, 112747020, 
                   207810010, 208021020, 210114010, 35201010, 
                   25391010)), "9weight_measured"] <- 0

#Drop RAND's weight variables
hrs_samp %<>% dplyr::select(-c(paste0("r", seq(8, 13, by = 1), "pmwght"), 
                               paste0("r", number_waves, "weight")))

#---- derived BMI ----
# #variable check
# weight <- hrs_samp %>% 
#   dplyr::select(paste0(seq(4, 9), "weight")) 
# missing_weight <- rowSums(is.na(weight))
# table(missing_weight, useNA = "ifany")

bmi_mat <- hrs_samp %>% 
  dplyr::select(paste0(seq(4, 9), "weight"))
for(i in 1:ncol(bmi_mat)){
  bmi_mat[, i] <- bmi_mat[, i]/(hrs_samp[, "height"])^2
}
bmi_mat %<>% set_colnames(paste0("r", seq(4, 9), "BMI"))
missing_bmi <- rowSums(is.na(bmi_mat))

# #Sanity check
# head(bmi_mat)
# head(hrs_samp %>% dplyr::select(c("height", paste0(seq(4, 9), "weight"))))

hrs_samp %<>% cbind(bmi_mat) %>% cbind(missing_bmi) 

#Drop RAND's BMI variables
hrs_samp %<>% dplyr::select(-c(paste0("r", number_waves, "bmi"), 
                               paste0("r", seq(8, 13, by = 1), "pmbmi")))

#---- smoking ----
hrs_samp %<>% 
  mutate("smoker" = hrs_samp %>% 
           dplyr::select(paste0("r", number_waves, "smoken")) %>%
           apply(1, function(x) 
           {if (sum(is.na(x)) == length(x)){
             return(NA)} else {x[min(which(!is.na(x)))]}}))
# #Sanity check
# View(hrs_samp[, c(paste0("r", number_waves, "smoken"), "smoker")])
# table(hrs_samp$smoker, useNA = "ifany")

#Drop RAND's smoking variables
hrs_samp %<>% dplyr::select(-paste0("r", number_waves, "smoken"))

#---- chronic conditions ----
# #Check missingness in wave-updated ever/never chronic conditions
# conditions <- c("diabe", "hibpe", "cancre", "lunge", "hearte", "stroke", 
#                 "memrye")
# 
# for(condition in conditions){
#   print(condition)
#   subset <- hrs_samp %>% dplyr::select(paste0("r", seq(4, 9), condition))
#   counts <- rowSums(is.na(subset))
#   print(table(counts, useNA = "ifany"))
# }
# view(hrs_samp[hrs_samp$r5hibpe == 6, c("HHIDPN", colnames(cond_mat))])
# #Check variability in chronic conditions
# for(condition in conditions){
#   print(condition)
#   subset <- hrs_samp %>% dplyr::select(paste0("r", c(4, 9), condition))
#   print(table(subset[, 1], subset[, 2]))
# }

# Imputing chronic conditions
#---- ** diabetes ----
hrs_samp <- impute_chronic_condition("diabe", paste0("r", seq(1, 9), "diabe"),
                                     seq(1, 9), hrs_samp)

#---- **high bp ----
hrs_samp <- impute_chronic_condition("hibpe", paste0("r", seq(1, 9), "hibpe"),
                                     seq(1, 9), hrs_samp)

#---- **cancer ----
hrs_samp <- impute_chronic_condition("cancre", paste0("r", seq(1, 9), "cancre"),
                                     seq(1, 9), hrs_samp)

#---- **lung ----
hrs_samp <- impute_chronic_condition("lunge", paste0("r", seq(1, 9), "lunge"),
                                     seq(1, 9), hrs_samp)

#---- **heart ----
hrs_samp <- impute_chronic_condition("hearte", paste0("r", seq(1, 9), "hearte"),
                                     seq(1, 9), hrs_samp)

#---- **stroke ----
hrs_samp <- impute_chronic_condition("stroke", paste0("r", seq(1, 9), "stroke"),
                                     seq(1, 9), hrs_samp)

#---- **memory ----
#For memory problems, data starts from wave 4.
hrs_samp <- impute_chronic_condition("memrye", paste0("r", seq(4, 9), "memrye"),
                                     seq(4, 9), hrs_samp)

#Sanity check
# 
# for(condition in conditions){
#   print(condition)
#   subset <- hrs_samp %>% dplyr::select(paste0("r", seq(4, 9), condition,
#                                               "_impute"))
#   counts <- rowSums(is.na(subset))
#   print(table(counts, useNA = "ifany"))
# }

# #---- data check ----
# #How many people are missing all data for a chronic condition
# missing_data_check <- c("diab",  "hibp", "cancr", "lung", "heart", "strok",
#                         "arthrs", "memry",
#                         "diabetes_rx_swallowed", "diabetes_rx_insulin",
#                         "bp_rx",  "lung_rx", "heart_rx", "stroke_rx",
#                         "arthritis_rx")
# 
# for(var in missing_data_check){
#   print(var)
#   if(grepl("rx", var)){
#     test <- hrs_samp %>%
#       dplyr::select(paste0(var, seq(4, 9)))
#     test_table <- table(rowSums(is.na(test)))
#     print(test_table)
#     print(sum(test_table))
#   } else{
#     test <- hrs_samp %>%
#       dplyr::select(paste0("r", seq(4, 9), var))
#     test_table <- table(rowSums(is.na(test)))
#     print(test_table)
#     print(sum(test_table))
#   }
# }
# 
# #Sanity check
# test <- hrs_samp %>%
#   dplyr::select(paste0("r", seq(4, 9), "diab"))
# test %<>% cbind(., rowSums(is.na(.)))
# table(test[, 7])
# sum(table(test[, 7]))

#---- sum of conditions ----
# wave-specific r(wave)conde
cond_mat <- hrs_samp %>%
  dplyr::select(contains("_impute"))

waves <- seq(1, 9)
for(j in 1:length(waves)){
  wave <- waves[j] 
  cond_mat[, paste0("r", wave , "conde", "_impute")] <- 
    rowSums(cond_mat %>% dplyr::select(contains(paste0("r", wave ))), 
            na.rm = TRUE)
}

hrs_samp[, colnames(cond_mat %>% select(contains("conde_impute")))] <- 
  cond_mat %>% select(contains("conde_impute"))

#---- marital status ----
# #Variable check
# table(hrs_samp$r4mstat, useNA = "ifany")
# table(hrs_samp$r9mstat, useNA = "ifany")

# Impute mstat with closest non-missing value
hrs_samp <- impute_status("mstat", paste0("r", seq(1, 13), "mstat"),
                          seq(1, 13), seq(4, 9), hrs_samp)

# Create marital status categories
mstat_mat <- hrs_samp %>% select(contains("mstat_impute"))

#Consolidate categories
mstat_mat[mstat_mat == 1 | mstat_mat == 2 | mstat_mat == 3] <- 1
mstat_mat[mstat_mat == 4 | mstat_mat == 5 | mstat_mat == 6 | 
            mstat_mat == 8] <- 2
mstat_mat[mstat_mat == 7] <- 3

hrs_samp[, colnames(mstat_mat)] <- mstat_mat

for(wave in seq(4, 9)){
  level = 1
  for(cat in c("married_partnered", "not_married_partnered", "widowed")){
    hrs_samp[, paste0("r", wave, cat)] <- 
      (hrs_samp[, paste0("r", wave, "mstat_impute")] == level)*1
    level = level + 1
  }
}

# #Sanity check
# View(hrs_samp %>% 
#        dplyr::select("r4mstat_impute", "r4married_partnered", 
#                      "r4not_married_partnered", "r4widowed"))
# View(hrs_samp %>% 
#        dplyr::select("r9mstat_impute", "r9married_partnered", 
#                      "r9not_married_partnered", "r9widowed"))

#---- drinking ----
# Imputing drinking per day and drinking per week
hrs_samp <- impute_status("drinkd", paste0("r", seq(3, 13), "drinkd"),
                          seq(3, 13), seq(4, 9),  hrs_samp)
hrs_samp <- impute_status("drinkn", paste0("r", seq(3, 13), "drinkn"),
                          seq(3, 13), seq(4, 9),  hrs_samp)
#Sanity check
# table(hrs_samp$r4drinkn, hrs_samp$r4drinkn_impute, useNA = "ifany")

drinks_per_week_mat <- (hrs_samp %>% dplyr::select(contains("drinkd_impute")))*
  (hrs_samp %>% dplyr::select(contains("drinkn_impute")))
ndrinks_mat <- hrs_samp %>% dplyr::select(contains("drinkn_impute"))

drinking_cat_mat <- 
  matrix(nrow = nrow(drinks_per_week_mat), ncol = ncol(drinks_per_week_mat)) %>% 
  set_colnames(paste0("r", seq(4, 9), "drinking_cat"))

for(i in 1:ncol(drinking_cat_mat)){
  for(j in 1:nrow(drinking_cat_mat)){
    drinking_cat_mat[j, paste0("r", (i + 3), "drinking_cat")] = 
      case_when(drinks_per_week_mat[j, i] == 0 ~ 0, 
                (drinks_per_week_mat[j, i] >= 7 | ndrinks_mat[j, i] >= 3) & 
                  hrs_samp[j, "female"] == 1 ~ 2, 
                (drinks_per_week_mat[j, i] >= 14 | ndrinks_mat[j, i] >= 4) & 
                  hrs_samp[j, "female"] == 0 ~ 2, 
                (drinks_per_week_mat[j, i] >= 1 & 
                   drinks_per_week_mat[j, i] < 7) & 
                  hrs_samp[j, "female"] == 1 ~ 1, 
                (drinks_per_week_mat[j, i] >= 1 & 
                   drinks_per_week_mat[j, i] < 14) & 
                  hrs_samp[j, "female"] == 0 ~ 1)
  }
}

hrs_samp %<>% cbind(drinking_cat_mat)

# #Sanity Check
# View(hrs_samp %>% dplyr::select("r9drinkn", "female",
#                                 "drinking9_cat") %>%
#        filter(drinking9_cat == "Heavy Drinking"))
# View(hrs_samp %>% dplyr::select("r9drinkn", "drinks_per_week9", "female",
#                                 "drinking9_cat") %>%
#        filter(drinking9_cat == "Heavy Drinking"))

# #---- Preliminary models ----
# # DO NOT DELETE DURING CODE CLEAN-UP
# # logit(P(missing CESD this wave)) =
# # \beta_0 + \beta_1*age at current wave + \beta_2*value of previous CESD +
# # \beta_3* chronic condition count (at last wave)
# 
# #num missing cesd among those sampled
# wave4_sampled <- hrs_samp %>% filter(!is.na(FALIVE))
# mean(wave4_sampled$r4cesd_missing)
# summary(r4cesdmissing_mod <- glm(r4cesd_missing ~
#                                    r4age_y_int + r3cesd + r3conde_impute +
#                                    death2018 - 1,
#                                  family = binomial(link = "logit"),
#                                  data = wave4_sampled))
# r4results <- tidy(r4cesdmissing_mod, exponentiate = TRUE, conf.int = TRUE)
# 
# wave5_sampled <- hrs_samp %>% filter(!is.na(GALIVE))
# mean(wave5_sampled$r5cesd_missing)
# summary(r5cesdmissing_mod <- glm(r5cesd_missing ~
#                                    r5age_y_int + r4cesd + r4conde_impute +
#                                    death2018 - 1,
#                                  family = binomial(link = "logit"),
#                                  data = wave5_sampled))
# r5results <- tidy(r5cesdmissing_mod, exponentiate = TRUE, conf.int = TRUE)
# 
# wave6_sampled <- hrs_samp %>% filter(!is.na(HALIVE))
# mean(wave6_sampled$r6cesd_missing)
# summary(r6cesdmissing_mod <- glm(r6cesd_missing ~
#                                    r6age_y_int + r5cesd + r5conde_impute +
#                                    death2018 - 1,
#                                  family = binomial(link = "logit"),
#                                  data = wave6_sampled))
# r6results <- tidy(r6cesdmissing_mod, exponentiate = TRUE, conf.int = TRUE)
# 
# wave7_sampled <- hrs_samp %>% filter(!is.na(JALIVE))
# mean(wave7_sampled$r7cesd_missing)
# summary(r7cesdmissing_mod <- glm(r7cesd_missing ~
#                                    r7age_y_int + r6cesd + r6conde_impute +
#                                    death2018 - 1,
#                                  family = binomial(link = "logit"),
#                                  data = wave7_sampled))
# r7results <- tidy(r7cesdmissing_mod, exponentiate = TRUE, conf.int = TRUE)
# 
# wave8_sampled <- hrs_samp %>% filter(!is.na(KALIVE))
# mean(wave8_sampled$r8cesd_missing)
# summary(r8cesdmissing_mod <- glm(r8cesd_missing ~
#                                    r8age_y_int + r7cesd + r7conde_impute +
#                                    death2018 - 1,
#                                  family = binomial(link = "logit"),
#                                  data = wave8_sampled))
# r8results <- tidy(r8cesdmissing_mod, exponentiate = TRUE, conf.int = TRUE)
# 
# wave9_sampled <- hrs_samp %>% filter(!is.na(LALIVE))
# mean(wave9_sampled$r9cesd_missing)
# summary(r9cesdmissing_mod <- glm(r9cesd_missing ~
#                                    r9age_y_int + r8cesd + r8conde_impute +
#                                    death2018 - 1,
#                                  family = binomial(link = "logit"),
#                                  data = wave9_sampled))
# r9results <- tidy(r9cesdmissing_mod, exponentiate = TRUE, conf.int = TRUE)
# 
# results_tbl <- tibble(
#   variables = c("age at current wave", "previous CESD value",
#                 "previous chronic condition count", "death2018"),
#   r4beta = round(r4results$estimate, 4),
#   r5beta = round(r5results$estimate, 4),
#   r6beta = round(r6results$estimate, 4),
#   r7beta = round(r7results$estimate, 4),
#   r8beta = round(r8results$estimate, 4),
#   r9beta = round(r9results$estimate, 4)
# )
# 
# results_tbl %>%
#   kbl(caption =
#         "Exponentiated betas of the CESD missing model (wave 4 - 9)") %>%
#   kable_classic(full_width = F, html_font = "Arial")
# 
# write_csv(results_tbl, paste0(path_to_dropbox,
#                               "/exposure_trajectories/data/",
#                               "CESD_missing_model_betas_in_sampled.csv"))
# 
# #---- **predict outcome?? ----
# # DO NOT DELETE DURING CODE CLEAN-UP
# summary(r4outcome_mod <- glm(death2018 ~
#                                r4age_y_int + r3cesd + r3conde_impute,
#                              family = binomial(link = "logit"),
#                              data = hrs_samp))
# r4outcome_results <- tidy(r4outcome_mod, exponentiate = TRUE, conf.int = TRUE)
# 
# summary(r5outcome_mod <- glm(death2018 ~
#                                r5age_y_int + r4cesd + r4conde_impute,
#                              family = binomial(link = "logit"),
#                              data = hrs_samp))
# r5outcome_results <- tidy(r5outcome_mod, exponentiate = TRUE, conf.int = TRUE)
# 
# summary(r6outcome_mod <- glm(death2018 ~
#                                r6age_y_int + r5cesd + r5conde_impute,
#                              family = binomial(link = "logit"),
#                              data = hrs_samp))
# r6outcome_results <- tidy(r6outcome_mod, exponentiate = TRUE, conf.int = TRUE)
# 
# summary(r7outcome_mod <- glm(death2018 ~
#                                    r7age_y_int + r6cesd + r6conde_impute,
#                                  family = binomial(link = "logit"),
#                                  data = hrs_samp))
# r7outcome_results <- tidy(r7outcome_mod, exponentiate = TRUE, conf.int = TRUE)
# 
# summary(r8outcome_mod <- glm(death2018 ~
#                                    r8age_y_int + r7cesd + r7conde_impute,
#                                  family = binomial(link = "logit"),
#                                  data = hrs_samp))
# r8outcome_results <- tidy(r8outcome_mod, exponentiate = TRUE, conf.int = TRUE)
# 
# summary(r9outcome_mod <- glm(death2018 ~
#                                    r9age_y_int + r8cesd + r8conde_impute,
#                                  family = binomial(link = "logit"),
#                                  data = hrs_samp))
# r9outcome_results <- tidy(r9outcome_mod, exponentiate = TRUE, conf.int = TRUE)
# 
# results_tbl <- tibble(
#   variables = c("Intercept", "age at current wave", "previous CESD value",
#                 "Previous chronic condition count"),
#   r4beta = round(r4outcome_results$estimate, 4),
#   r5beta = round(r5outcome_results$estimate, 4),
#   r6beta = round(r6outcome_results$estimate, 4),
#   r7beta = round(r7outcome_results$estimate, 4),
#   r8beta = round(r8outcome_results$estimate, 4),
#   r9beta = round(r9outcome_results$estimate, 4)
# )
# 
# results_tbl %>%
#   kbl(caption = "Exponentiated betas association with death2018 missing model (wave 4 - 9)") %>%
#   kable_classic(full_width = F, html_font = "Arial")
# 
# write_csv(results_tbl, paste0(path_to_dropbox,
#                               "/exposure_trajectories/data/",
#                               "predict_death2018_betas.csv"))

#---- Dropping people ----
# 1. full HRSsample (n = 42233)

# 2. remove people who were not sampled in 2018
sum(is.na(hrs_samp$QALIVE))
hrs_samp %<>% filter(!is.na(QALIVE))

# 3. drop those with missing CESD observations in wave 4 - 9
drop <- hrs_samp %>% dplyr::select(paste0("r", seq(4, 9, by = 1), "cesd")) %>% 
  mutate("drop" = apply(., 1, function(x) sum(is.na(x)) > 0)*1)
table(drop$drop, useNA = "ifany")
hrs_samp %<>% mutate("drop" = drop$drop) %>% filter(drop == 0)

# 4.1 Check those missing age data-- these people have no birth date data so I am 
# dropping them; looks like no one is in this group
still_missing <- 
  which(is.na(rowSums(hrs_samp %>% dplyr::select(contains("age_y")))))
sum(is.na(hrs_samp[still_missing, "Bday"])) == length(still_missing)
if(length(still_missing) > 0){
  hrs_samp <- hrs_samp[-c(still_missing), ]
}

# 4.2 Drop those not age-eligible at HRS wave 4 and those who are 91+
sum(!hrs_samp$r4age_y_int %in% c(seq(50, 90)))
hrs_samp %<>% filter(r4age_y_int %in% c(seq(50, 90)))

# 5. Dropping those missing race/ethnicity data
sum(!hrs_samp$unknown_race_eth == 0)
hrs_samp %<>% filter(unknown_race_eth == 0)

# 6. Drop people missing height data
#Drop RAND's height variables + extra derived variables
# sum(is.na(hrs_samp$height))
hrs_samp %<>% filter(!is.na(height)) %>% 
  dplyr::select(-c(paste0("r", seq(8, 13, by = 1), "pmhght"), 
                   paste0("r", number_waves, "height"), 
                   "med_height", "self_height"))

# 7. Drop those missing BMI （weight)
sum(hrs_samp$missing_bmi != 0)
hrs_samp %<>% filter(missing_bmi == 0)

# 8. Drop those who miss drinking status at any wave after imputation
subset <- hrs_samp %>% dplyr::select(contains("drinking"))
table(rowSums(is.na(subset)))
hrs_samp %<>% filter(rowSums(is.na(subset)) == 0)

# 9. Drop those without self-reported health
drop <- rowSums(is.na(hrs_samp %>%
                        dplyr::select(paste0("r", seq(4, 9), "shlt"))))
# table(drop)
hrs_samp[, "drop"] <- drop
# #Sanity check
# table(hrs_samp$drop, useNA = "ifany")
hrs_samp %<>% filter(drop == 0)

#9. Drop people missing wave-updated ever/never chronic conditions
conditions <- c("diabe", "hibpe", "cancre", "lunge", "hearte", "stroke",
                "memrye")
subset <- hrs_samp %>% 
  dplyr::select(do.call(paste, 
                        c(expand.grid("r", seq(4, 9), conditions, "_impute"), 
                          sep = ""))) 
hrs_samp %<>% mutate("drop" = rowSums(is.na(subset))) %>% filter(drop == 0)

# Checking other variables missingness

# No missing education data
# table(hrs_samp$ed_cat, useNA = "ifany")

# #No missing marital status for waves 5-8
# subset <- hrs_samp %>% dplyr::select(paste0("r", seq(5, 8), "mstat_cat"))
# drop <- rowSums(is.na(subset))

#---- select variables ----
vars <- c("HHIDPN", paste0("r", seq(4, 9), "married_partnered"),
          paste0("r", seq(4, 9), "not_married_partnered"),
          paste0("r", seq(4, 9), "widowed"), "ed_cat", 
          paste0("r", seq(4, 9), "drinking_cat"), 
          paste0("r", seq(4, 9), "memrye_impute"),
          paste0("r", seq(4, 9), "stroke_impute"),
          paste0("r", seq(4, 9), "hearte_impute"),
          paste0("r", seq(4, 9), "lunge_impute"),
          paste0("r", seq(4, 9), "cancre_impute"),
          paste0("r", seq(4, 9), "hibpe_impute"),
          paste0("r", seq(4, 9), "diabe_impute"),
          paste0("r", seq(3, 9), "conde_impute"), "smoker", 
          paste0("r", seq(4, 9), "BMI"), "hispanic", "white", "black", "other", 
          "female", paste0("r", seq(4, 9), "age_y_int"), "death2018", 
          paste0("r", seq(3, 9), "cesd"), paste0("r", seq(3, 9), "shlt"), 
          "age_death_y")

hrs_samp %<>% dplyr::select(all_of(vars))

#---- check missingness ----
#Should have no missingness except in 
# r3cesd (not part of optimal wave range)
# age_death_y (for those who are still living)
# self reported health  (rwshlt, since no longer using this in the models)
colSums(is.na(hrs_samp))

#---- Exposures ----
#---- **E1 Def: CESD at HRS wave 4 (1998) ----
#Effect of E1a on survival to HRS wave 14 (2018) 
hrs_samp %<>% 
  mutate("r4cesd_elevated" = ifelse(r4cesd >= 4, 1, 0), 
         "r4cesd_elevated_sens" = ifelse(r4cesd >= 1, 1, 0)) 

# #Sanity check
# table(hrs_samp$r4cesd, hrs_samp$r4cesd_elevated, useNA = "ifany")

#---- **E2 Def: CESD at HRS wave 9 (2008) ----
#Effect of E1a on survival to HRS wave 14 (2018) 
hrs_samp %<>% 
  mutate("r9cesd_elevated" = ifelse(r9cesd >= 4, 1, 0), 
         "r9cesd_elevated_sens" = ifelse(r9cesd >= 1, 1, 0))

# #Sanity check
# table(hrs_samp$r9cesd, hrs_samp$r9cesd_elevated, useNA = "ifany")

#---- **E3 Def: Cumulative Exposure (proportion of occasions) ----
#Number of occasions with elevated depressive symptoms in HRS waves 4-9
hrs_samp %<>% 
  mutate("prop_elevated_cesd" = hrs_samp %>% 
           dplyr::select(paste0("r", seq(4, 9), "cesd")) %>% 
           mutate_all(function(x) ifelse(x >= 4, 1, 0)) %>% rowMeans(),
         "prop_elevated_cesd_sens" = hrs_samp %>% 
           dplyr::select(paste0("r", seq(4, 9), "cesd")) %>% 
           mutate_all(function(x) ifelse(x >= 1, 1, 0)) %>% rowMeans())

# #Sanity check
# head(elevated_cesd)
# head(hrs_samp$prop_elevated_cesd)
# table(hrs_samp$prop_elevated_cesd, useNA = "ifany")

#---- **E4 Def: Cumulative Exposure (average CESD score) ----
hrs_samp %<>% 
  mutate("avg_cesd" = hrs_samp %>% 
           dplyr::select(paste0("r", seq(4, 9, by = 1), "cesd")) %>% 
           rowMeans(), 
         "avg_cesd_elevated" = ifelse(avg_cesd >= 4, 1, 0), 
         "avg_cesd_elevated_sens" = ifelse(avg_cesd >= 1, 1, 0))

# #Sanity check
# View(hrs_samp %<>% 
#        dplyr::select(paste0("r", seq(4, 9, by = 1), "cesd"), "avg_cesd", 
#                      "avg_cesd_elevated"))
# plot(hrs_samp$avg_cesd, hrs_samp$avg_cesd_elevated)

#---- Outcome ----
#---- ** survival times from HRS wave 9 (2008) to HRS wave 14 (2018) ----
hrs_samp %<>% mutate(survtime = age_death_y - r9age_y_int) %>% 
  mutate(survtime = ifelse(is.na(survtime), 10, survtime), 
         observed = ifelse(is.na(age_death_y), 0, 1))

# #Sanity check
# View(hrs_samp %>% dplyr::select("age_death_y", "survtime", "observed"))

#---- extra variables ----
#intercept for fast mice methods
hrs_samp %<>% mutate("intercept" = 1)

#variables for missingness mechanisms
for(wave in seq(4, 9)){
  hrs_samp %<>% 
    dplyr::mutate(!!paste0("r", wave, "cesd_death2018") := 
                    !!sym(paste0("r", wave, "cesd"))*death2018, 
                  !!paste0("r", wave - 1, "cesd_conde_impute") := 
                    !!sym(paste0("r", wave - 1, "cesd"))*
                    !!sym(paste0("r", wave - 1, "conde_impute")))
}

#---- save dataset ----
write_csv(hrs_samp, paste0(path_to_dropbox,
                           "/exposure_trajectories/data/",
                           "CESD_data_wide.csv"))

#---- OLD ----
#source(here::here("RScripts", "non_missing.R"))
#source(here::here("RScripts", "read_da_dct.R"))
#source(here::here("Rscripts", "chronic_condition.R"))

#---- Check the variability of drinking behavior ----
# #Run this script after dropping people!
# # Check variability for drinking behavior in the analytic sample
# 
# temp <- hrs_samp %>% 
#   select(HHIDPN, paste0("r", seq(4, 9), "drinkd"), 
#          paste0("r", seq(4, 9), "drinkn"),
#          paste0("r", seq(4, 9), "drinkd_impute"), 
#          paste0("r", seq(4, 9), "drinkn_impute"),
#          paste0("r", seq(4, 9), "drinking_cat"),
#          female) 
# 
# for (wave in seq(4, 9)){
#   temp %<>%
#     dplyr::rename(!!paste0("r", wave, "drink_cat_impute") :=
#                     paste0("r", wave, "drinking_cat"))
#   temp[, paste0("r", wave, "drinkw")] <- 
#     temp[, paste0("r", wave, "drinkd")] *
#     temp[, paste0("r", wave, "drinkn")]
#   temp[, paste0("r", wave, "drinkw_impute")] <- 
#     temp[, paste0("r", wave, "drinkd_impute")] *
#     temp[, paste0("r", wave, "drinkn_impute")]
#   temp[, paste0("r", wave, "drink_cat")] <- 
#     case_when(temp[, "female"] ==  1 & 
#                 (temp[, paste0("r", wave, "drinkw")] >= 7 | 
#                    temp[, paste0("r", wave, "drinkn")] >= 3) ~ 2,
#               temp[, "female"] ==  1 & 
#                 (temp[, paste0("r", wave, "drinkw")] <= 7 & 
#                    temp[, paste0("r", wave, "drinkw")] >= 1) ~ 1,
#               temp[, paste0("r", wave, "drinkw")] == 0 ~ 0,
#               temp[, "female"] ==  0 & 
#                 (temp[, paste0("r", wave, "drinkw")] >= 14 | 
#                    temp[, paste0("r", wave, "drinkn")] >= 4) ~ 2,
#               temp[, "female"] ==  0 & 
#                 (temp[, paste0("r", wave, "drinkw")] <= 14 & 
#                    temp[, paste0("r", wave, "drinkw")] >= 1) ~ 1)
# }
# # Somehow the !!paste0("r", wave, "drinkw") := !!sym(paste0("r", wave, "drinkd"))
# # * !!sym(paste0("r", wave, "drinkn")) is not working
# 
# p_load("scales")
# n <- length(levels(as.factor(hrs_samp$HHIDPN)))# number of colors
# cols <- hue_pal(h = c(0, 360) + 15, 
#                 c = 100, l = 65, 
#                 h.start = 0, direction = 1)(n)[order(sample(1:n, n))] 
# 
# frac <- 0.1
# 
# temp %>%
#   select(HHIDPN, paste0("r", seq(4, 9), "drinkw")) %>%
#   sample_frac(frac, replace = FALSE) %>%
#   pivot_longer(!HHIDPN,
#                names_to = "wave",
#                names_pattern = "r(.*)drinkw",
#                values_to = "drink") %>%
#   mutate(HHIDPN = as.factor(HHIDPN)) %>%
#   ggplot(aes(x = wave, y = drink, group = HHIDPN, color = HHIDPN)) +
#   geom_point() + geom_line() +
#   scale_y_continuous(limits=c(0, 20), breaks = seq(0, 14, 7)) +
#   theme(legend.position = "none") +
#   scale_color_manual(values = cols) +
#   labs(title = "# of Drinks per week across waves", 
#        subtitle = paste0("A ", frac*100, 
#                          "% random sample of the analytic sample"))
# 
# temp %>%
#   select(HHIDPN, paste0("r", seq(4, 9), "drinkw_impute")) %>%
#   sample_frac(frac, replace = FALSE) %>%
#   pivot_longer(!HHIDPN,
#                names_to = "wave",
#                names_pattern = "r(.*)drinkw_impute",
#                values_to = "drink_impute") %>%
#   mutate(HHIDPN = as.factor(HHIDPN)) %>%
#   ggplot(aes(x = wave, y = drink_impute, group = HHIDPN, color = HHIDPN)) +
#   geom_point() + geom_line() +
#   scale_y_continuous(limits=c(0, 20), breaks = seq(0, 14, 7)) +
#   theme(legend.position = "none") +
#   scale_color_manual(values = cols) +
#   labs(title = "# of Drinks per week across waves (imputed)", 
#        subtitle = paste0("A ", frac*100, 
#                          "% random sample of the analytic sample"))
# 
# with(temp %>%
#        mutate(drink_change = case_when(
#          abs(r5drink_cat - r4drink_cat) == 2|
#            abs(r6drink_cat - r5drink_cat) == 2|
#            abs(r7drink_cat - r6drink_cat) == 2|
#            abs(r8drink_cat - r7drink_cat) == 2|
#            abs(r9drink_cat - r8drink_cat) == 2 ~ 2,
#          abs(r5drink_cat - r4drink_cat) == 1|
#            abs(r6drink_cat - r5drink_cat) == 1|
#            abs(r7drink_cat - r6drink_cat) == 1|
#            abs(r8drink_cat - r7drink_cat) == 1|
#            abs(r9drink_cat - r8drink_cat) == 1 ~ 1,
#          TRUE ~ 0),
#          drink_change_impute = case_when(
#            abs(r5drink_cat_impute - r4drink_cat_impute) == 2|
#              abs(r6drink_cat_impute - r5drink_cat_impute) == 2|
#              abs(r7drink_cat_impute - r6drink_cat_impute) == 2|
#              abs(r8drink_cat_impute - r7drink_cat_impute) == 2|
#              abs(r9drink_cat_impute - r8drink_cat_impute) == 2 ~ 2,
#            abs(r5drink_cat_impute - r4drink_cat_impute) == 1|
#              abs(r6drink_cat_impute - r5drink_cat_impute) == 1|
#              abs(r7drink_cat_impute - r6drink_cat_impute) == 1|
#              abs(r8drink_cat_impute - r7drink_cat_impute) == 1|
#              abs(r9drink_cat_impute - r8drink_cat_impute) == 1 ~ 1,
#            TRUE ~ 0),), table(drink_change, 
#                               drink_change_impute, useNA = "ifany"))

# #---- physical activity ----
# PA_mat <- hrs_samp %>% 
#   dplyr::select(contains("ltactx"), contains("mdactx"), contains("vgactx"))
# 
# waves_PA <-  seq(7, 13, by = 1)
# 
# for(i in 1:nrow(PA_mat)){
#   for(j in 1:length(waves_PA)){
#     wave <- waves_PA[j]
#     PA_mat[i, paste0("r", wave, "MET")] = 
#       case_when(PA_mat[i, paste0("r", wave, "vgactx")] == 1 ~ 17,
#                 PA_mat[i, paste0("r", wave, "vgactx")] == 2 ~ 13,
#                 PA_mat[i, paste0("r", wave, "vgactx")] %in% c(3,4,5,NA) &
#                   PA_mat[i, paste0("r", wave, "mdactx")] == 1 ~ 10，
#                 PA_mat[i, paste0("r", wave, "vgactx")] == 3 &
#                   PA_mat[i, paste0("r", wave, "mdactx")] %in% c(2,3,4,5,NA) ~ 9,
#                 PA_mat[i, paste0("r", wave, "vgactx")] %in% c(4,5,NA) &
#                   PA_mat[i, paste0("r", wave, "mdactx")] == 2 ~ 7.5，
#                 PA_mat[i, paste0("r", wave, "vgactx")] %in% c(4,5,NA) &
#                   PA_mat[i, paste0("r", wave, "mdactx")] == 3 ~ 5，
#                 PA_mat[i, paste0("r", wave, "vgactx")] == 4 &
#                   PA_mat[i, paste0("r", wave, "mdactx")] %in% c(4,5,NA) ~ 4.25,
#                 PA_mat[i, paste0("r", wave, "vgactx")] %in% c(5, NA) &
#                   PA_mat[i, paste0("r", wave, "mdactx")] %in% c(4,5,NA) &
#                   PA_mat[i, paste0("r", wave, "ltactx")] == 1 ~ 4,
#                 PA_mat[i, paste0("r", wave, "vgactx")] %in% c(5, NA) &
#                   PA_mat[i, paste0("r", wave, "mdactx")] %in% c(4,5,NA) &
#                   PA_mat[i, paste0("r", wave, "ltactx")] == 2 ~ 3,
#                 PA_mat[i, paste0("r", wave, "vgactx")] %in% c(5, NA) &
#                   PA_mat[i, paste0("r", wave, "mdactx")] == 4 &
#                   PA_mat[i, paste0("r", wave, "ltactx")] %in% c(3, 4,5,NA) ~ 2.5,
#                 PA_mat[i, paste0("r", wave, "vgactx")] %in% c(5, NA) &
#                   PA_mat[i, paste0("r", wave, "mdactx")] %in% c(5,NA) &
#                   PA_mat[i, paste0("r", wave, "ltactx")] == 3 ~ 2,
#                 PA_mat[i, paste0("r", wave, "vgactx")] %in% c(5, NA) &
#                   PA_mat[i, paste0("r", wave, "mdactx")] %in% c(5,NA) &
#                   PA_mat[i, paste0("r", wave, "ltactx")] == 4 ~ 1,
#                 ifelse(is.na(PA_mat[i, paste0("r", wave, "vgactx")]),1,0) == 1 &
#                   ifelse(is.na(PA_mat[i, paste0("r", wave, "mdactx")]),1,0) == 1 &
#                   ifelse(is.na(PA_mat[i, paste0("r", wave, "ltactx")]),1,0) == 1
#                 ~ NA_real_,
#                 PA_mat[i, paste0("r", wave, "vgactx")] %in% c(5, NA) &
#                   PA_mat[i, paste0("r", wave, "mdactx")] %in% c(5,NA) &
#                   PA_mat[i, paste0("r", wave, "ltactx")] %in% c(5,NA) ~ 0
#       )
#   }
# }
# 
# # Sanity check
# # View(PA_mat %>% dplyr::select("r8ltactx", "r8mdactx", "r8vgactx", "r8MET") %>%
# #        filter(is.na(r8MET)))
# # View(PA_mat %>% dplyr::select("r13ltactx", "r13mdactx", "r13vgactx", "r13MET") %>%
# #        filter(is.na(r13MET)))
# 
# # Combine to hrs_samp
# hrs_samp %<>% cbind(PA_mat[, paste0("r", waves_PA, "MET")])

# #---- self-reported health ----
# #Variable check
# self_reported_health <- hrs_samp %>%
#   dplyr::select(paste0("r", seq(4, 9), "shlt"))
# colSums(is.na(self_reported_health))
#---- ** drinking ----
# Old code chunk
# drinks_per_week_mat <- (hrs_samp %>% dplyr::select(contains("drinkd")))*
#   (hrs_samp %>% dplyr::select(contains("drinkn")))
# ndrinks_mat <- hrs_samp %>% dplyr::select(contains("drinkn"))
# 
# 
# drinking_cat_mat <- 
#   matrix(nrow = nrow(drinks_per_week_mat), ncol = ncol(drinks_per_week_mat)) %>% 
#   set_colnames(paste0("drinking", seq(3, 13), "_cat"))
# 
# for(i in 1:ncol(drinking_cat_mat)){
#   for(j in 1:nrow(drinking_cat_mat)){
#     drinking_cat_mat[j, paste0("drinking", (i + 2), "_cat")] = 
#       case_when(drinks_per_week_mat[j, i] == 0 ~ "No Drinking", 
#                 (drinks_per_week_mat[j, i] >= 7 | ndrinks_mat[j, i] >= 3) & 
#                   hrs_samp[j, "female"] == 1 ~ "Heavy Drinking", 
#                 (drinks_per_week_mat[j, i] >= 14 | ndrinks_mat[j, i] >= 4) & 
#                   hrs_samp[j, "female"] == 0 ~ "Heavy Drinking", 
#                 (drinks_per_week_mat[j, i] >= 1 & 
#                    drinks_per_week_mat[j, i] < 7) & 
#                   hrs_samp[j, "female"] == 1 ~ "Moderate Drinking", 
#                 (drinks_per_week_mat[j, i] >= 1 & 
#                    drinks_per_week_mat[j, i] < 14) & 
#                   hrs_samp[j, "female"] == 0 ~ "Moderate Drinking")
#   }
# }

# #Sanity Check
# View(drinking_cat_mat)

# #Impute drinking4_cat and drinking9_cat with closest non-missing values
# hrs_samp %<>% 
#   mutate("drinking4_cat_impute" = 
#            ifelse(is.na(drinking4_cat), drinking3_cat, drinking4_cat)) %>% 
#   mutate("drinking4_cat_impute" = 
#            ifelse(is.na(drinking4_cat_impute), 
#                   hrs_samp %>% 
#                     dplyr::select(paste0("drinking", seq(5, 13), "_cat")) %>% 
#                     apply(., 1, function(x) x[min(which(!is.na(x)))]), 
#                   drinking4_cat_impute)) %>% 
#   mutate("drinking9_cat_impute" = 
#            ifelse(is.na(drinking9_cat), 
#                   hrs_samp %>% 
#                     dplyr::select(paste0("drinking", seq(3, 8), "_cat")) %>% 
#                     apply(., 1, function(x) x[max(which(!is.na(x)))]), 
#                   drinking9_cat)) %>% 
#   mutate("drinking9_cat_impute" = 
#            ifelse(is.na(drinking9_cat_impute), 
#                   hrs_samp %>% 
#                     dplyr::select(paste0("drinking", seq(10, 13), "_cat")) %>% 
#                     apply(., 1, function(x) x[min(which(!is.na(x)))]), 
#                   drinking9_cat_impute))

# #Sanity check
# View(hrs_samp %>% dplyr::select(contains("drinking", seq(3, 13))))
# View(hrs_samp %>% dplyr::select(contains("drinking", seq(3, 13))) %>% 
#        filter(is.na(drinking4_cat) | is.na(drinking9_cat)))
# table(hrs_samp$drinking4_cat_impute, useNA = "ifany")
# table(hrs_samp$drinking9_cat_impute, useNA = "ifany")
# table(hrs_samp$drinking4_cat, hrs_samp$drinking9_cat)

# #Drop the one person who has no information on drinking behavior
# hrs_samp %<>% filter(!is.na(drinking4_cat_impute))
#---- ** marital status ----
# Old code chunk for mstat
#Impute r4mstat with closest non-missing value
# hrs_samp %<>%
#   mutate("r4mstat_impute" =
#            ifelse(is.na(r4mstat),
#                   hrs_samp %>%
#                     dplyr::select(paste0("r", seq(1, 3), "mstat")) %>%
#                     apply(., 1, function(x) x[max(which(!is.na(x)))]),
#                   r4mstat)) %>%
#   mutate("r4mstat_impute" =
#            ifelse(is.na(r4mstat_impute),
#                   hrs_samp %>%
#                     dplyr::select(paste0("r", seq(5, 13), "mstat")) %>%
#                     apply(., 1, function(x) x[min(which(!is.na(x)))]),
#                   r4mstat_impute)) %>%


# #Sanity check
# sum(hrs_samp$r4mstat != hrs_samp$r4mstat_impute, na.rm = TRUE)
# View(hrs_samp %>%
#        dplyr::select("HHIDPN", paste0("r", number_waves, "mstat"),
#                      "r4mstat_impute"))
# View(hrs_samp %>%
#        dplyr::select("HHIDPN", paste0("r", number_waves, "mstat"),
#                      "r4mstat_impute") %>% filter(is.na(r4mstat)))

# #---- **arthritis ----
# hrs_samp <-
#   impute_chronic_condition("arthre", paste0("r", seq(4, 9), "arthre"),
#                               seq(4, 9), hrs_samp)

# # #---- Old code chunk ----
# # #---- ** diabetes ----
# hrs_samp <- chronic_condition("diabetes", paste0("r", seq(1, 13), "diab"),
#                               c(paste0("diabetes_rx_insulin", seq(4, 9)),
#                                 paste0("diabetes_rx_swallowed", seq(4, 9))),
#                               hrs_samp)
# 
# # #Sanity check
# # for(var in c(paste0("r", seq(1, 13), "diab"),
# #              paste0("diabetes_rx_insulin", seq(4, 9)))){
# #   print(var)
# #   print(table(hrs_samp[, var], useNA = "ifany"))
# # }
# 
# #---- **high bp ----
# hrs_samp <- chronic_condition("hibp", paste0("r", seq(1, 13), "hibp"),
#                               paste0("bp_rx", seq(4, 9)), hrs_samp)
# 
# #---- **cancer ----
# hrs_samp <- chronic_condition("cancer", paste0("r", seq(1, 13), "cancr"),
#                               NA, hrs_samp)
# 
# #---- **lung ----
# hrs_samp <- chronic_condition("lung", paste0("r", seq(1, 13), "lung"),
#                               paste0("lung_rx", seq(4, 9)), hrs_samp)
# 
# #---- **heart ----
# hrs_samp <- chronic_condition("heart", paste0("r", seq(1, 13), "heart"),
#                               paste0("heart_rx", seq(4, 9)), hrs_samp)
# 
# #---- **stroke ----
# hrs_samp <- chronic_condition("stroke", paste0("r", seq(1, 13), "strok"),
#                               paste0("stroke_rx", seq(4, 9)), hrs_samp)
# 
# #---- **arthritis ----
# hrs_samp <- chronic_condition("arthritis", paste0("r", seq(2, 13), "arthrs"),
#                               paste0("arthritis_rx", seq(4, 9)), hrs_samp)
# ## arthritis is reported from wave 2 to wave 13.
# 
# #---- **memory ----
# hrs_samp <- chronic_condition("mem", paste0("r", seq(4, 9), "memry"),
#                               NA, hrs_samp)
# 
# #---- **sum of conditions ----
# #We're going to create our own version of r[wave]conde from RAND, but ours will
# #   not be wave-specific
# cond_mat <- hrs_samp %>%
#   dplyr::select(contains("ever"))
# 
# #since ever_condition variables are used to derive new conde variable,
# #it should not vary across waves
# hrs_samp[, "conde"] <- rowSums(cond_mat, na.rm = TRUE)

# #---- read in HRS Core files ---- 
# #needed for medication
# 
# # '98-'00 medication in section B, 02-later medication in section C
# # memory problem medication starts from wave 2008; we won't do memory meds
# # since it starts on the last wave relevant to our analysis
# 
# years <- c("98", "00", 
#            "02", "04", "06", "08") #medication
# dataframes_list <- vector(mode = "list", length = (length(years)))
# 
# for(i in 1:length(years)){
#   year <- years[i]
#   if(i == 1 | i == 2){
#     dataframes_list[[i]] <- 
#       read_da_dct(paste0(path_to_box, "/Box/HRS/core_files/h", year,
#                          "core/h", year, "da/H", year, "B_R.da"),
#                   paste0(path_to_box, "/Box/HRS/core_files/h", year,
#                          "core/h", year, "sta/H", year, "B_R.dct"),
#                   skip_lines = 1,
#                   HHIDPN = TRUE) %>% mutate_at("HHIDPN", as.numeric) %>% 
#       #select variables of interest
#       dplyr::select("HHIDPN", 
#                     #1998 vars         2000 vars
#                     contains("F1117"), contains("G1248"),
#                     contains("F1118"), contains("G1249"),
#                     contains("F1110"), contains("G1239"),
#                     contains("F1151"), contains("G1284"),
#                     contains("F1157"), contains("G1290"),
#                     contains("F1184"), contains("G1317")) %>%
#       #contains("F1198"), contains("G1331")) 
#       set_colnames(c("HHIDPN", 
#                      paste0("diabetes_rx_swallowed", (i + 3)),
#                      paste0("diabetes_rx_insulin", (i + 3)), 
#                      paste0("bp_rx", (i + 3)), 
#                      paste0("lung_rx", (i + 3)), 
#                      paste0("heart_rx", (i + 3)), 
#                      paste0("stroke_rx", (i + 3)) 
#                      #paste0("arthritis_rx", (i + 3))
#       )) 
#   } else{
#     dataframes_list[[i]] <-
#       read_da_dct(paste0(path_to_box, "/Box/HRS/core_files/h", year,
#                          "core/h", year, "da/H", year, "C_R.da"),
#                   paste0(path_to_box, "/Box/HRS/core_files/h", year,
#                          "core/h", year, "sta/H", year, "C_R.dct"), 
#                   skip_lines = 2, 
#                   HHIDPN = TRUE) %>% mutate_at("HHIDPN", as.numeric) %>% 
#       #select variables of interest
#       dplyr::select("HHIDPN", 
#                     contains("C011"), contains("C012"), contains("C006"), 
#                     contains("C032"), contains("C037"), contains("C060")) %>% 
#       #contains("C074")) %>%
#       set_colnames(c("HHIDPN", 
#                      paste0("diabetes_rx_swallowed", (i + 3)),
#                      paste0("diabetes_rx_insulin", (i + 3)), 
#                      paste0("bp_rx", (i + 3)), 
#                      paste0("lung_rx", (i + 3)), 
#                      paste0("heart_rx", (i + 3)), 
#                      paste0("stroke_rx", (i + 3))
#                      #paste0("arthritis_rx", (i + 3))
#       ))
#   }
# }