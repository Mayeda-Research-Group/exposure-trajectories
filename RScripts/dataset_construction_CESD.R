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
# Yingyan's directory: C:/Users/yingyan_wu
#                      C:/Users/yingyan_wu/Dropbox
# Crystal's directory: /Users/CrystalShaw
#                     ~/Dropbox/Projects

#Changing directories here will change them throughout the script
path_to_box <- "/Users/CrystalShaw"
path_to_dropbox <- "~/Dropbox/Projects"

#---- source scripts ----
source(here::here("RScripts", "non_missing.R"))
source(here::here("RScripts", "impute_ages.R"))
source(here::here("RScripts", "measured_self_report.R"))
source(here::here("RScripts", "read_da_dct.R"))
source(here::here("Rscripts", "chronic_condition.R"))
source(here::here("Rscripts", "impute_condition.R"))

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
  read_sas(paste0(path_to_box, "/Box/HRS/tracker/trk2018_3/", 
                  "trk2018tr_r.sas7bdat")) %>% 
  select("HHID", "PN", "PIWTYPE", "PALIVE", "QIWTYPE", "QALIVE") %>% 
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
                    # paste0("r", seq(2, 13, by = 1), "arthrs"),
                    # paste0("r", seq(2, 13, by = 1), "arthre"),
                    paste0("r", number_waves, "conde"),
                    paste0("r", number_waves, "smoken"), 
                    paste0("r", seq(3, 13, by = 1), "drinkd"),
                    paste0("r", seq(3, 13, by = 1), "drinkn"),
                    paste0("r", seq(7, 13, by = 1), "vgactx"),
                    paste0("r", seq(7, 13, by = 1), "mdactx"), 
                    paste0("r", seq(7, 13, by = 1), "ltactx"),
                    paste0("r", seq(2, 13, by = 1), "cesd"),
                    paste0("r", number_waves, "shlt"))

RAND <- read_dta(paste0(path_to_box, "/Box/HRS/RAND_longitudinal/STATA/", 
                        "randhrs1992_2016v2.dta"), 
                 col_select = all_of(rand_variables)) 

colnames(RAND)[1] <- "HHIDPN" #For merging

#Remove labeled data format
val_labels(RAND) <- NULL

#---- read in HRS Core files ---- 
#needed for medication

# '98-'00 medication in section B, 02-later medication in section C
# memory problem medication starts from wave 2008; we won't do memory meds
# since it starts on the last wave relevant to our analysis

years <- c("98", "00", 
           "02", "04", "06", "08") #medication
dataframes_list <- vector(mode = "list", length = (length(years)))

for(i in 1:length(years)){
  year <- years[i]
  if(i == 1 | i == 2){
    dataframes_list[[i]] <- 
      read_da_dct(paste0(path_to_box, "/Box/HRS/core_files/h", year,
                         "core/h", year, "da/H", year, "B_R.da"),
                  paste0(path_to_box, "/Box/HRS/core_files/h", year,
                         "core/h", year, "sta/H", year, "B_R.dct"),
                  skip_lines = 1,
                  HHIDPN = TRUE) %>% mutate_at("HHIDPN", as.numeric) %>% 
      #select variables of interest
      dplyr::select("HHIDPN", 
                    #1998 vars         2000 vars
                    contains("F1117"), contains("G1248"),
                    contains("F1118"), contains("G1249"),
                    contains("F1110"), contains("G1239"),
                    contains("F1151"), contains("G1284"),
                    contains("F1157"), contains("G1290"),
                    contains("F1184"), contains("G1317"),
                    contains("F1198"), contains("G1331")) %>%
      set_colnames(c("HHIDPN", 
                     paste0("diabetes_rx_swallowed", (i + 3)),
                     paste0("diabetes_rx_insulin", (i + 3)), 
                     paste0("bp_rx", (i + 3)), 
                     paste0("lung_rx", (i + 3)), 
                     paste0("heart_rx", (i + 3)), 
                     paste0("stroke_rx", (i + 3)), 
                     paste0("arthritis_rx", (i + 3))
                     )) 
  } else{
    dataframes_list[[i]] <-
      read_da_dct(paste0(path_to_box, "/Box/HRS/core_files/h", year,
                         "core/h", year, "da/H", year, "C_R.da"),
                  paste0(path_to_box, "/Box/HRS/core_files/h", year,
                         "core/h", year, "sta/H", year, "C_R.dct"), 
                  skip_lines = 2, 
                  HHIDPN = TRUE) %>% mutate_at("HHIDPN", as.numeric) %>% 
      #select variables of interest
      dplyr::select("HHIDPN", 
                    contains("C011"), contains("C012"), contains("C006"), 
                    contains("C032"), contains("C037"), contains("C060"), 
                    contains("C074")) %>%
      set_colnames(c("HHIDPN", 
                     paste0("diabetes_rx_swallowed", (i + 3)),
                     paste0("diabetes_rx_insulin", (i + 3)), 
                     paste0("bp_rx", (i + 3)), 
                     paste0("lung_rx", (i + 3)), 
                     paste0("heart_rx", (i + 3)), 
                     paste0("stroke_rx", (i + 3)), 
                     paste0("arthritis_rx", (i + 3))
                     ))
  }
}

#---- read in Anusha Vable's CSES index ----
cSES <- read_dta(paste0(path_to_dropbox, "/exposure_trajectories/data/", 
                        "cSES measures/cses_measures.dta"), 
                 col_select = c("hhid", "pn", "cses_index")) %>% 
  unite(col = "HHIDPN", c("hhid", "pn"), sep = "") %>% 
  mutate_at("HHIDPN", as.numeric)

#---- merge datasets ----
#Use this to subset RAND data
hrs_samp <- join_all(c(list(hrs_tracker, RAND, cSES), dataframes_list), 
                     by = "HHIDPN", type = "left") 
             
#---- remove people who were not sampled in 2018 ----
hrs_samp %<>% filter(!is.na(QALIVE))

#---- looking for optimal subset ----
# #Drop those who are not age-eligible for HRS at the start of follow-up
# subsets_data <- data.frame(matrix(nrow = 45, ncol = 8)) %>%
#   set_colnames(c("CESD_start_wave", "CESD_end_wave", "num_measures",
#                  "sample_size", "min_age", "max_age", "death_2018",
#                  "prop_dead"))
# 
# index = 0
# for(i in 2:9){
#   for(j in (i + 4):13){
#     index = index + 1
#     subsets_data[index, c("CESD_start_wave", "CESD_end_wave")] = c(i,j)
#     subsets_data[index, "num_measures"] = j - i + 1
#     
#     data_subset <- hrs_samp %>%
#       dplyr::select(paste0("r", seq(i, j, by = 1), "cesd"), "death2018",
#                     paste0(i, "age_y_int")) %>%
#       na.omit()
#     data_subset[, "too_young"] =
#       ifelse(data_subset[, tail(colnames(data_subset), n = 1)] < 50, 1, 0)
#     
#     data_subset %<>% filter(too_young == 0)
#     
#     subsets_data[index, "sample_size"] = nrow(data_subset)
#     subsets_data[index, "min_age"] = min(data_subset[, paste0(i, "age_y_int")])
#     subsets_data[index, "max_age"] = max(data_subset[, paste0(i, "age_y_int")])
#     subsets_data[index, "death_2018"] = sum(data_subset$death2018)
#     subsets_data[index, "prop_dead"] = mean(data_subset$death2018)
#   }
# }
# 
# #Best subset was waves 4-9
# write_csv(subsets_data, here::here("Prelim Analyses", "exp_CESD_out_mortality",
#                              "CESD_complete_subsets.csv"))

drop <- hrs_samp %>% dplyr::select(paste0("r", seq(4, 9, by = 1), "cesd")) %>% 
  mutate("drop" = apply(., 1, function(x) sum(is.na(x)) > 0)*1)

hrs_samp %<>% mutate("drop" = drop$drop) %>%
  #drop those with missing CESD observations in these waves
  filter(drop == 0) 

#---- death ----
#death indicators: RAND dates of death get us to 2016 use QALIVE for 2018
hrs_samp %<>% mutate("death2016" = ifelse(is.na(raddate), 0, 1), 
                     "death2018" = ifelse(QALIVE %in% c(5, 6), 1, 0))

# #Sanity check
# #1 = alive; 2 = presumed alive; 5 = known deceased this wave; 
# #6 = known deceased prior wave 
# table(hrs_samp$QALIVE, useNA = "ifany")
# table(hrs_samp$death2016, useNA = "ifany")
# table(hrs_samp$death2018, useNA = "ifany")

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

#Check those missing age data-- these people have no birth date data so I am 
# dropping them; looks like no one is in this group
still_missing <- 
  which(is.na(rowSums(hrs_samp %>% dplyr::select(contains("age_y")))))
sum(is.na(hrs_samp[still_missing, "Bday"])) == length(still_missing)
if(length(still_missing) > 0){
  hrs_samp <- hrs_samp[-c(still_missing), ]
}

#Impute date of death for those who are dead in 2018
hrs_samp %<>% 
  mutate("age_death_y" = ifelse((is.na(age_death_y) & death2018 == 1), 
                                r13age_y_int + 2, age_death_y))

#Drop those not age-eligible at HRS wave 4 and those who are 91+
hrs_samp %<>% filter(r4age_y_int %in% c(seq(50, 90)))

# #Sanity check
# View(hrs_samp[, c("age_death_y", "death2016", "death2018", "13age_y_int")])
# table(hrs_samp$age_death_y, useNA = "ifany")
# sum(table(hrs_samp$age_death_y))
# table(hrs_samp$QALIVE %in% c(5, 6))
# #I think there's a one person discrepancy between RAND's death data and 
# # HRS's QALIVE variable-- check this after all the other data cleaning steps
# table(hrs_samp$death2018, hrs_samp$age_death_y)

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

#1 person missing race/ethnicity data
hrs_samp %<>% filter(unknown_race_eth == 0) %>% 
  #Drop the RAND rahispan variable (recoded as hispanic) and race variables
  dplyr::select(-c("rahispan", "raracem"))

#---- education ----
#Create education categories
hrs_samp %<>% 
  mutate("ed_cat" = case_when(raedyrs < 12 ~ "Less than HS", 
                              raedyrs == 12 ~ "HS", 
                              raedyrs > 12 & raedyrs < 16 ~ "Some College", 
                              raedyrs == 16 ~ "Bachelors", 
                              raedyrs > 16 ~ "Grad Studies"))

# #Sanity check
# table(is.na(hrs_samp$raedyrs))
# table(hrs_samp$raedyrs)
# table(hrs_samp$raedyrs, hrs_samp$ed_cat, useNA = "ifany")

#---- cSES index ----
# #Sanity check
# table(is.na(hrs_samp$cses_index))

#none missing!

#---- height ----
#Create a "best" height variable by taking the median of measured heights 
# (waves 8+) if available or first self-reported height
hrs_samp %<>% 
  mutate("med_height" = hrs_samp %>% 
           dplyr::select(paste0("r", seq(8, 13, by = 1), "pmhght")) %>%
           apply(1, function(x) median(x, na.rm = TRUE)), 
         "self_height" = hrs_samp %>% 
           dplyr::select(paste0("r", number_waves, "height")) %>%
           apply(1, function(x) x[min(which(!is.na(x)))])) %>% 
  mutate("height_measured" = ifelse(!is.na(med_height), 1, 0)) %>% 
  mutate("height" = ifelse(height_measured == 1, med_height, self_height))

# #Sanity check
# View(hrs_samp[, c(paste0("r", seq(8, 13, by = 1), "pmhght"),
#                   paste0("r", number_waves, "height"),
#                   "med_height", "self_height", "height_measured", "height")] %>%
#        filter(is.na(height)))

#Drop people missing height data
#Drop RAND's height variables + extra derived variables
hrs_samp %<>% filter(!is.na(height)) %>% 
  dplyr::select(-c(paste0("r", seq(8, 13, by = 1), "pmhght"), 
                               paste0("r", number_waves, "height"), 
                               "med_height", "self_height"))

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
#variable check
weight <- hrs_samp %>% 
  dplyr::select(paste0(seq(4, 9), "weight")) 
missing_weight <- rowSums(is.na(weight))
table(missing_weight, useNA = "ifany")

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
hrs_samp %<>% filter(missing_bmi == 0)

#Drop RAND's BMI variables
hrs_samp %<>% dplyr::select(-c(paste0("r", number_waves, "bmi"), 
                               paste0("r", seq(8, 13, by = 1), "pmbmi")))

#---- smoking ----
hrs_samp %<>% 
  mutate("smoker" = hrs_samp %>% 
           dplyr::select(paste0("r", number_waves, "smoken")) %>%
           apply(1, function(x) x[min(which(!is.na(x)))]))

# #Sanity check
# View(hrs_samp[, c(paste0("r", number_waves, "smoken"), "smoker")])
# table(hrs_samp$smoker, useNA = "ifany")

#Drop RAND's smoking variables
hrs_samp %<>% dplyr::select(-paste0("r", number_waves, "smoken"))

#---- chronic conditions ----
# Imputing chronic conditions

#---- ** diabetes ----
hrs_samp <- impute_chronic_condition("diabe", paste0("r", seq(4,9), "diabe"),
                              seq(4,9), hrs_samp)

#---- **high bp ----
hrs_samp <- impute_chronic_condition("hibpe", paste0("r", seq(4, 9), "hibpe"),
                             seq(4, 9), hrs_samp)

#---- **cancer ----
hrs_samp <- impute_chronic_condition("cancre", paste0("r", seq(4, 9), "cancre"),
                                     seq(4, 9), hrs_samp)

#---- **lung ----
hrs_samp <- impute_chronic_condition("lunge", paste0("r", seq(4, 9), "lunge"),
                                     seq(4, 9), hrs_samp)

#---- **heart ----
hrs_samp <- impute_chronic_condition("hearte", paste0("r", seq(4, 9), "hearte"),
                                     seq(4, 9), hrs_samp)

#---- **stroke ----
hrs_samp <- impute_chronic_condition("stroke", paste0("r", seq(4, 9), "stroke"),
                                     seq(4, 9), hrs_samp)

# #---- **arthritis ----
# hrs_samp <- 
#   impute_chronic_condition("arthre", paste0("r", seq(4, 9), "arthre"),
#                               seq(4, 9), hrs_samp)

#---- **memory ----
hrs_samp <- impute_chronic_condition("memrye", paste0("r", seq(4, 9), "memrye"),
                                     seq(4, 9), hrs_samp)

#sanity check
# table(hrs_samp$r5memrye_impute, hrs_samp$r5memrye, useNA = "ifany")
# table(hrs_samp$r6memrye_impute, hrs_samp$r6memrye, useNA = "ifany")
# table(hrs_samp$r7memrye_impute, hrs_samp$r7memrye, useNA = "ifany")
# table(hrs_samp$r8memrye_impute, hrs_samp$r8memrye, useNA = "ifany")
# table(hrs_samp$r9memrye_impute, hrs_samp$r9memrye, useNA = "ifany")
# View(hrs_samp %>% dplyr::select(contains(paste0("r", seq(4, 9), "memrye")),
#                              -contains("rx"))%>%
#        filter(!rowSums(is.na(.)) == 0))

#---- sum of conditions ----
# wave-specific r(wave)conde
cond_mat <- hrs_samp %>%
  dplyr::select(contains("_impute"), -contains("drinking"))

waves <- seq(4,9)
  for(j in 1:length(waves)){
    wave <- waves[j] 
     cond_mat[, paste0("r", wave , "conde", "_impute")] <- 
       rowSums(cond_mat %>% dplyr::select(contains(paste0("r", wave ))), 
               na.rm = TRUE)
  }

hrs_samp[, colnames(cond_mat %>% select(contains("conde_impute")))] <- 
  cond_mat %>% select(contains("conde_impute"))

# view(hrs_samp %>% select(contains("conde_impute")))

# # #---- Old code chunk
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
# #---- sum of conditions ----
# #We're going to create our own version of r[wave]conde from RAND, but ours will
# #   not be wave-specific
# cond_mat <- hrs_samp %>%
#   dplyr::select(contains("ever"))
# 
# #since ever_condition variables are used to derive new conde variable,
# #it should not vary across waves
# hrs_samp[, "conde"] <- rowSums(cond_mat, na.rm = TRUE)



#---- marital status ----
# #Variable check
# table(hrs_samp$r4mstat, useNA = "ifany")
# table(hrs_samp$r9mstat, useNA = "ifany")

#Impute r4mstat with closest non-missing value
hrs_samp %<>% 
  mutate("r4mstat_impute" = 
           ifelse(is.na(r4mstat), 
                  hrs_samp %>% 
                    dplyr::select(paste0("r", seq(1, 3), "mstat")) %>% 
                    apply(., 1, function(x) x[max(which(!is.na(x)))]), 
                  r4mstat)) %>% 
  mutate("r4mstat_impute" = 
           ifelse(is.na(r4mstat_impute), 
                  hrs_samp %>% 
                    dplyr::select(paste0("r", seq(5, 13), "mstat")) %>% 
                    apply(., 1, function(x) x[min(which(!is.na(x)))]), 
                  r4mstat_impute))

# #Sanity check
# sum(hrs_samp$r4mstat != hrs_samp$r4mstat_impute, na.rm = TRUE)
# View(hrs_samp %>%
#        dplyr::select("HHIDPN", paste0("r", number_waves, "mstat"),
#                      "r4mstat_impute"))
# View(hrs_samp %>%
#        dplyr::select("HHIDPN", paste0("r", number_waves, "mstat"),
#                      "r4mstat_impute") %>% filter(is.na(r4mstat)))
  
#Create marital status categories
hrs_samp %<>% 
  mutate("r9mstat_cat" = 
           case_when(r9mstat %in% c(1, 2, 3) ~ "Married/Partnered", 
                     r9mstat %in% c(4, 5, 6, 8) ~ "Not Married/Partnered", 
                     r9mstat == 7 ~ "Widowed"), 
         "r4mstat_cat" = 
           case_when(r4mstat_impute %in% c(1, 2, 3) ~ "Married/Partnered", 
                     r4mstat_impute %in% c(4, 5, 6, 8) ~ "Not Married/Partnered", 
                     r4mstat_impute == 7 ~ "Widowed"))

# #Sanity check
# table(hrs_samp$r4mstat_impute, hrs_samp$r4mstat_cat, useNA = "ifany")
# table(hrs_samp$r9mstat, hrs_samp$r9mstat_cat, useNA = "ifany")

#---- drinking ----
drinks_per_week_mat <- (hrs_samp %>% dplyr::select(contains("drinkd")))*
  (hrs_samp %>% dplyr::select(contains("drinkn")))
ndrinks_mat <- hrs_samp %>% dplyr::select(contains("drinkn"))

drinking_cat_mat <- 
  matrix(nrow = nrow(drinks_per_week_mat), ncol = ncol(drinks_per_week_mat)) %>% 
  set_colnames(paste0("drinking", seq(3, 13), "_cat"))

for(i in 1:ncol(drinking_cat_mat)){
  for(j in 1:nrow(drinking_cat_mat)){
    drinking_cat_mat[j, paste0("drinking", (i + 2), "_cat")] = 
      case_when(drinks_per_week_mat[j, i] == 0 ~ "No Drinking", 
                (drinks_per_week_mat[j, i] >= 7 | ndrinks_mat[j, i] >= 3) & 
                  hrs_samp[j, "female"] == 1 ~ "Heavy Drinking", 
                (drinks_per_week_mat[j, i] >= 14 | ndrinks_mat[j, i] >= 4) & 
                  hrs_samp[j, "female"] == 0 ~ "Heavy Drinking", 
                (drinks_per_week_mat[j, i] >= 1 & 
                   drinks_per_week_mat[j, i] < 7) & 
                  hrs_samp[j, "female"] == 1 ~ "Moderate Drinking", 
                (drinks_per_week_mat[j, i] >= 1 & 
                   drinks_per_week_mat[j, i] < 14) & 
                  hrs_samp[j, "female"] == 0 ~ "Moderate Drinking")
  }
}

# #Sanity Check
# View(drinking_cat_mat)

hrs_samp %<>% cbind(drinking_cat_mat)

# #Sanity Check
# View(hrs_samp %>% dplyr::select("r9drinkn", "drinks_per_week9", "female",
#                                 "drinking9_cat") %>%
#        filter(drinking9_cat == "Heavy Drinking"))

# #Variable check
# table(hrs_samp$drinking4_cat, useNA = "ifany")
# table(hrs_samp$drinking9_cat, useNA = "ifany")

#Impute drinking4_cat and drinking9_cat with closest non-missing values
hrs_samp %<>% 
  mutate("drinking4_cat_impute" = 
           ifelse(is.na(drinking4_cat), drinking3_cat, drinking4_cat)) %>% 
  mutate("drinking4_cat_impute" = 
           ifelse(is.na(drinking4_cat_impute), 
                  hrs_samp %>% 
                    dplyr::select(paste0("drinking", seq(5, 13), "_cat")) %>% 
                    apply(., 1, function(x) x[min(which(!is.na(x)))]), 
                  drinking4_cat_impute)) %>% 
  mutate("drinking9_cat_impute" = 
           ifelse(is.na(drinking9_cat), 
                  hrs_samp %>% 
                    dplyr::select(paste0("drinking", seq(3, 8), "_cat")) %>% 
                    apply(., 1, function(x) x[max(which(!is.na(x)))]), 
                  drinking9_cat)) %>% 
  mutate("drinking9_cat_impute" = 
           ifelse(is.na(drinking9_cat_impute), 
                  hrs_samp %>% 
                    dplyr::select(paste0("drinking", seq(10, 13), "_cat")) %>% 
                    apply(., 1, function(x) x[min(which(!is.na(x)))]), 
                  drinking9_cat_impute))

# #Sanity check
# View(hrs_samp %>% dplyr::select(contains("drinking", seq(3, 13))))
# View(hrs_samp %>% dplyr::select(contains("drinking", seq(3, 13))) %>% 
#        filter(is.na(drinking4_cat) | is.na(drinking9_cat)))
# table(hrs_samp$drinking4_cat_impute, useNA = "ifany")
# table(hrs_samp$drinking9_cat_impute, useNA = "ifany")

#Drop the one person who has no information on drinking behavior
hrs_samp %<>% filter(!is.na(drinking4_cat_impute))

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

# #Sanity check
# test <- hrs_samp %>% 
#   dplyr::select(paste0("r", seq(4, 9), "diab")) 
# test %<>% cbind(., rowSums(is.na(.)))  
# table(test[, 7])
# sum(table(test[, 7]))

#---- self-reported health ----
#Variable check
self_reported_health <- hrs_samp %>%
  dplyr::select(paste0("r", seq(4, 9), "shlt"))
colSums(is.na(self_reported_health))
drop <- rowSums(is.na(self_reported_health))
#table(drop)

hrs_samp[, "drop"] <- drop

# #Sanity check
# table(hrs_samp$drop, useNA = "ifany")

hrs_samp %<>% filter(drop == 0)

#---- select variables ----
vars <- c("HHIDPN", paste0("r", c(4, 9), "mstat_cat"), "ed_cat", 
          paste0("drinking", c(4, 9), "_cat_impute"), 
          paste0("r", seq(4, 9), "memrye", "_impute"),
          paste0("r", seq(4, 9), "stroke", "_impute"),
          paste0("r", seq(4, 9), "hearte", "_impute"),
          paste0("r", seq(4, 9), "lunge", "_impute"),
          paste0("r", seq(4, 9), "cancre", "_impute"),
          paste0("r", seq(4, 9), "hibpe", "_impute"),
          paste0("r", seq(4, 9), "diabe", "_impute"),
          paste0("r", seq(4, 9), "conde", "_impute"),
          "smoker", 
          paste0("r", seq(4, 9), "BMI"), "hispanic", "white", "black", "other", 
          "female", paste0("r", seq(4, 9), "age_y_int"), "death2018", 
          paste0("r", seq(4, 9), "cesd"), paste0("r", seq(4, 9), "shlt"), 
          "age_death_y")

hrs_samp %<>% dplyr::select(all_of(vars))
                        
#---- save dataset ----
write_csv(hrs_samp, paste0(path_to_dropbox,
                           "/exposure_trajectories/data/",
                           "hrs_samp_6CESD_waves4-9.csv"))


