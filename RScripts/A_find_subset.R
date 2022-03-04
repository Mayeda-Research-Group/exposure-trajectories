# A. Find the optimal subset (which will have most measures as our truth dataset)
# Data input: trk2018tr_r.sas7bdat, randhrs1992_2018v1.dta, cses_measures.dta
# Data output: CESD_complete_subsets.csv
# Author: CS & YW

#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}
p_load("here", "readr", "tidyverse", "magrittr", "plyr", "haven", "labelled")

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
                    paste0("r", number_waves, "hibp"),
                    paste0("r", number_waves, "hibpe"),
                    paste0("r", number_waves, "diab"),
                    paste0("r", number_waves, "diabe"),
                    paste0("r", number_waves, "cancr"),
                    paste0("r", number_waves, "cancre"),
                    paste0("r", number_waves, "strok"), 
                    paste0("r", number_waves, "stroke"),
                    paste0("r", number_waves, "heart"),
                    paste0("r", number_waves, "hearte"),
                    paste0("r", number_waves, "lung"),
                    paste0("r", number_waves, "lunge"),
                    paste0("r", seq(4 ,9 , by = 1), "memry"),
                    paste0("r", seq(4 ,9 , by = 1), "memrye"),
                    #We don't think this is an important confounder of 
                    # CESD --> Mortality
                    # paste0("r", seq(2, 13, by = 1), "arthrs"),
                    # paste0("r", seq(2, 13, by = 1), "arthre"),
                    paste0("r", number_waves, "conde"),
                    paste0("r", number_waves, "smoken"), 
                    paste0("r", seq(3, 13, by = 1), "drinkd"),
                    paste0("r", seq(3, 13, by = 1), "drinkn"),
                    #We don't think this is an important confounder of 
                    # CESD --> Mortality
                    # paste0("r", seq(7, 13, by = 1), "vgactx"),
                    # paste0("r", seq(7, 13, by = 1), "mdactx"), 
                    # paste0("r", seq(7, 13, by = 1), "ltactx"),
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

#---- Data cleaning ----
hrs_samp %<>% mutate("death2018" = ifelse(is.na(raddate), 0, 1))
#---- **age ----
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

#---- looking for optimal subset ----
#Drop those who are not age-eligible for HRS at the start of follow-up
subsets_data <- data.frame(matrix(nrow = 36, ncol = 8)) %>%
  set_colnames(c("CESD_start_wave", "CESD_end_wave", "num_measures",
                 "sample_size", "min_age", "max_age", "death_2018",
                 "prop_dead"))

index = 0
for(i in 2:9){
  for(j in (i + 4):13){
    index = index + 1
    subsets_data[index, c("CESD_start_wave", "CESD_end_wave")] = c(i,j)
    subsets_data[index, "num_measures"] = j - i + 1
    
    data_subset <- hrs_samp %>%
      dplyr::select(paste0("r", seq(i, j, by = 1), "cesd"), "death2018",
                    paste0("r", i, "age_y_int")) %>%
      na.omit()
    data_subset[, "too_young"] =
      ifelse(data_subset[, tail(colnames(data_subset), n = 1)] < 50, 1, 0)
    
    data_subset %<>% filter(too_young == 0)
    
    subsets_data[index, "sample_size"] = nrow(data_subset)
    subsets_data[index, "min_age"] = 
      min(data_subset[, paste0("r", i, "age_y_int")])
    subsets_data[index, "max_age"] = 
      max(data_subset[, paste0("r", i, "age_y_int")])
    subsets_data[index, "death_2018"] = sum(data_subset$death2018)
    subsets_data[index, "prop_dead"] = mean(data_subset$death2018)
  }
}

#Best subset was waves 4-9
write_csv(subsets_data, here::here("Prelim Analyses", "exp_CESD_out_mortality",
                                   "CESD_complete_subsets.csv"))

