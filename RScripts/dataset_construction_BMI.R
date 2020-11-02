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

# years <- c("06", "08", "10", "12", "14") #biomarker sample
# letter_waves <- LETTERS[seq(from = 11, to = 15)] #biomarker sample 
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
##         reports heart problems this wave
##        
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
                    paste0("r", number_waves, "diab"),
                    paste0("r", number_waves, "diabe"),
                    paste0("r", number_waves, "cancr"),
                    paste0("r", number_waves, "strok"), 
                    paste0("r", number_waves, "heart"),
                    paste0("r", number_waves, "smoken"), 
                    paste0("r", seq(3, 13, by = 1), "drinkd"),
                    paste0("r", seq(3, 13, by = 1), "drinkn"),
                    paste0("r", seq(7, 13, by = 1), "vgactx"),
                    paste0("r", seq(7, 13, by = 1), "mdactx"), 
                    paste0("r", seq(7, 13, by = 1), "ltactx"))

RAND <- read_dta(paste0(path_to_box, "/Box/HRS/RAND_longitudinal/STATA/", 
                        "randhrs1992_2016v2.dta"), 
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
hrs_samp <- join_all(list(hrs_tracker, RAND, cSES), 
                     by = "HHIDPN", type = "left") 

#---- remove people who were not sampled in 2018 ----
hrs_samp %<>% filter(!is.na(QALIVE))

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
hrs_samp[, paste0(number_waves, "age_y")] <- age_m/12
#Ages rounded down to nearest year
hrs_samp[, paste0(number_waves, "age_y_int")] <- floor(age_m/12)

# #Sanity check
# View(hrs_samp[, c(paste0(number_waves, "age_y"), 
#                   paste0(number_waves, "age_y_int"))])

#Check those missing age data-- these people have no birthdate data so I am 
# dropping them
still_missing <- 
  which(is.na(rowSums(hrs_samp %>% dplyr::select(contains("age_y")))))
sum(is.na(hrs_samp[still_missing, "Bday"])) == length(still_missing)
hrs_samp <- hrs_samp[-c(still_missing), ]

#Impute data of death for those who are dead in 2018
hrs_samp %<>% 
  mutate("age_death_y" = ifelse((is.na(age_death_y) & death2018 == 1), 
                                `13age_y_int` + 2, age_death_y))

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

#26 people missing race/ethnicity data
hrs_samp %<>% filter(unknown_race_eth == 0) %>% 
  #Drop the RAND rahispan variable (recoded as hispanic) and race variables
  dplyr::select(-c("rahispan", "raracem"))

#---- education ----
# #Sanity check
# table(is.na(hrs_samp$raedyrs))

#There are 114 people missing years of education, so I'm going to drop them
hrs_samp %<>% filter(!is.na(raedyrs))

#---- cSES index ----
# #Sanity check
# table(is.na(hrs_samp$cses_index))

#There are 11,165 people missing the cSES index so I am dropping them
hrs_samp %<>% filter(!is.na(cses_index))

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
#                   "med_height", "self_height", "height_measured", "height")])
  
#Drop RAND's height variables + extra derived variables
hrs_samp %<>% dplyr::select(-c(paste0("r", seq(8, 13, by = 1), "pmhght"), 
                               paste0("r", number_waves, "height"), 
                               "med_height", "self_height"))

#---- weight ----
hrs_samp %<>% 
  cbind(measured_self_report(hrs_samp, 
                             paste0("r", seq(8, 13, by = 1), "pmwght"), 
                             paste0("r", number_waves, "weight"), "weight", 
                             number_waves))

# #Checking weird weight values (from YW's analysis)
# View(hrs_samp %>% filter(HHIDPN %in% c(20480010, 203788021)) %>%
#        dplyr::select(c(paste0(letter_waves, "weight"),
#                        paste0(letter_waves, "weight_measured"))))

#Set the weird measures to NA
hrs_samp[which(hrs_samp$HHIDPN == 20480010), "Mweight"] <- NA
hrs_samp[which(hrs_samp$HHIDPN == 203788021), "Nweight"] <- NA

#Fix measured indicators
hrs_samp[which(hrs_samp$HHIDPN == 20480010), "Mweight_measured"] <- 0
hrs_samp[which(hrs_samp$HHIDPN == 203788021), "Nweight_measured"] <- 0

#Drop RAND's weight variables
hrs_samp %<>% dplyr::select(-c(paste0("r", number_waves, "pmwght"), 
                               paste0("r", number_waves, "weight")))

#---- BMI ----
hrs_samp %<>% 
  cbind(measured_self_report(hrs_samp, paste0("r", number_waves, "pmbmi"), 
                             paste0("r", number_waves, "bmi"), "BMI"))

# #Sanity check-- I had an issue with this observation that led to my finding
# #               a bug in my measured_self_report code
# View(hrs_samp %>% filter(HHIDPN %in% c(164907010)) %>%
#        dplyr::select(c(paste0(letter_waves, "BMI"),
#                        paste0(letter_waves, "BMI_measured"))))
       
# #Checking weird BMI values (from YW's analysis)
# View(hrs_samp %>% filter(HHIDPN %in% c(20480010, 203788021)) %>%
#        dplyr::select(c(paste0(letter_waves, "BMI"),
#                        paste0(letter_waves, "BMI_measured"))))

#Set the weird measures to NA-- consistent with weight 
hrs_samp[which(hrs_samp$HHIDPN == 20480010), "MBMI"] <- NA
hrs_samp[which(hrs_samp$HHIDPN == 203788021), "NBMI"] <- NA

#Fix measured indicators
hrs_samp[which(hrs_samp$HHIDPN == 20480010), "MBMI_measured"] <- 0
hrs_samp[which(hrs_samp$HHIDPN == 203788021), "NBMI_measured"] <- 0

#Drop RAND's BMI variables
hrs_samp %<>% dplyr::select(-c(paste0("r", number_waves, "bmi"), 
                               paste0("r", number_waves, "pmbmi")))

#---- smoking ----
hrs_samp %<>% 
  mutate("smoker" = hrs_samp %>% 
           dplyr::select(paste0("r", number_waves, "smoken")) %>%
           apply(1, function(x) x[min(which(!is.na(x)))]))

# #Sanity check
# View(hrs_samp[, c(paste0("r", number_waves, "smoken"), "smoker")])

#Drop RAND's smoking variables
hrs_samp %<>% dplyr::select(-paste0("r", number_waves, "smoken"))

#---- meds ----
med_vars <- colnames(hrs_samp)[str_detect(colnames(hrs_samp), pattern = "rx")]

reformat <- hrs_samp[, med_vars]

# #Sanity check
# which(reformat == 5)
# which(reformat == 8)
# which(reformat == 9)

reformat[reformat == 5] <- 0
reformat[reformat == 8 | reformat == 9] <- NA

hrs_samp[, med_vars] <- reformat

#---- save dataset ----
#3 Cystatin C measures
write_csv(hrs_samp, paste0(path_to_dropbox,
                           "/exposure_trajectories/data/",
                           "hrs_samp_3cysc.csv"))

# #Survival through age 70 and at least one cystatin c measure in [60, 70)
# write_csv(hrs_samp, paste0("C:/Users/yingyan_wu/Dropbox//",
#                            "exposure_trajectories/data/",
#                            "hrs_samp_alive_70_cysc_60_70.csv"))

# #Survival through age 70
# write_csv(hrs_samp, paste0("C:/Users/yingyan_wu/Dropbox/",
#                            "exposure_trajectories/data/",
#                            "hrs_samp_alive_70.csv"))


