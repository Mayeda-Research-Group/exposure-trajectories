#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "readr", "tidyverse", "magrittr", "plyr", "haven", "labelled", 
       "lubridate")

#No scientific notation
options(scipen = 999)

#---- source scripts ----
source(here::here("RScripts", "read_da_dct.R"))
source(here::here("RScripts", "non_missing.R"))
source(here::here("RScripts", "fu_time.R"))
source(here::here("RScripts", "impute_ages.R"))

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

years <- c("06", "08", "10", "12", "14")
letter_waves <- LETTERS[seq(from = 11, to = 15)]
number_waves <- seq(8, 12, by = 1)

#---- read in data ----
hrs_tracker <- 
  read_da_dct("/Users/CrystalShaw/Box/HRS/2016_tracker/trk2016/TRK2016TR_R.da", 
              "/Users/CrystalShaw/Box/HRS/2016_tracker/trk2016/TRK2016TR_R.dct", 
              HHIDPN = TRUE) %>% 
  #Don't use this when trying to figure out survival through age 70
  # #participated in HRS 2016 wave (wave "P")
  # filter(PIWTYPE %in% c(1, 5, 11, 15)) %>% 
  # filter(PALIVE %in% c(1, 5)) %>%
  select("HHIDPN", "PIWTYPE")  

#Don't use this when trying to figure out survival through age 70
# %>% 
#   #people alive at the beginning of 2014 HRS interview wave
#   filter(KNOWNDECEASEDYR >= 2014 | is.na(KNOWNDECEASEDYR)) %>% 
#   filter(EXDEATHYR >= 2014 | is.na(EXDEATHYR))

#2006-2012 biomarker data
biomarker_list <- vector(mode = "list", length = length(years))

for(i in 1:length(years)){
  year <- years[i]
  
  if(i != length(years)){
    biomarker_list[[i]] <-  
      read_da_dct(paste0("/Users/CrystalShaw/Box/HRS/biomarker_data/biomkr", 
                         year, "/BIOMK", year, "BL_R.da"),
                  paste0("/Users/CrystalShaw/Box/HRS/biomarker_data/biomkr", 
                         year, "/BIOMK", year, "BL_R.dct"), HHIDPN = TRUE)
  } else{
    #2014 early release biomarker data
    biomarker_list[[i]] <-  read_da_dct(
      "/Users/CrystalShaw/Box/HRS/biomarker_data/BIOMK14BL/BIOMK14BL.da",
      "/Users/CrystalShaw/Box/HRS/biomarker_data/BIOMK14BL/BIOMK14BL.dct", 
      HHIDPN = TRUE)
  }
}

#RAND longitudinal file-- reading in STATA file because SAS file wouldn't load
#Variables of interest: 
## Demographics: HHIDPN, gender, race, hispanic, birth month, 
##               birth year, birth date, death month, death year, death date, 
##               age data in months (biomarker waves), years of education, 
##               highest degree (masked)
## Health Behaviors: current smoker
## 
# Note: Dates are formatted as SAS dates (days from January 1, 1960)

rand_variables <- c("hhidpn", "ragender", "raracem", "rahispan", "rabmonth", 
                    "rabyear", "rabdate", "radmonth", "radyear", "raddate", 
                    paste0("r", number_waves, "agem_e"), "raedyrs",
                    "raedegrm",
                    paste0("r", number_waves, "smoken"))

RAND <- read_dta("~/Box/HRS/RAND_longitudinal/STATA/randhrs1992_2016v2.dta", 
                 col_select = all_of(rand_variables)) %>% 
  mutate_at("hhidpn", as.factor)

colnames(RAND)[1] <- "HHIDPN" #For merging

#Remove labeled data format
val_labels(RAND) <- NULL

#HRS Core files-- Need for physical measures
#2006-2014 core files
core_list <- vector(mode = "list", length = length(years))

for(i in 1:length(years)){
  year <- years[i]
  core_list[[i]] <- 
    read_da_dct(paste0("/Users/CrystalShaw/Box/HRS/core_files/h", year, 
                       "core/h", year, "da/H", year, "I_R.da"),
                paste0("/Users/CrystalShaw/Box/HRS/core_files/h", year, 
                       "core/h", year, "sta/H", year, "I_R.dct"), 
                HHIDPN = TRUE)
}

#---- merge core and biomarker data across waves ----
core_merge <- join_all(core_list, by = "HHIDPN", type = "left") %>% 
  #Select variables of interest: ID, Weight (pounds), Height (inches), 
  #                              systolic bp (3 times), 
  #                              diastolic bp (3 times)
  dplyr::select(HHIDPN, contains("I841"), contains("I834"), 
                contains("I859"), contains("I864"), contains("I869"), 
                contains("I860"), contains("I865"), contains("I870")) %>% 
  set_colnames(c("HHIDPN", paste0(letter_waves, "wt"), 
                 paste0(letter_waves, "ht"), 
                 paste0(letter_waves, "sbp", rep(seq(1, 3), each = 5)), 
                 paste0(letter_waves, "dbp", rep(seq(1, 3), each = 5))))

biomarker_merge <- join_all(biomarker_list, by = "HHIDPN", type = "left") %>% 
  #Select variables of interest: ID, adjusted Cystatin C, adjusted HbA1c,
  #                              adjusted total cholesterol, adjusted HDL
  dplyr::select(HHIDPN, contains("CYSC_ADJ"), contains("A1C_ADJ"), 
                contains("TC_ADJ"), contains("HDL_ADJ"))

#---- merge datasets ----
#Use this to subset RAND data
hrs_samp <- join_all(list(hrs_tracker, core_merge, RAND, biomarker_merge), 
                     by = "HHIDPN", type = "left")
  
#---- death ----
#death indicator
hrs_samp %<>% mutate("death" = ifelse(is.na(raddate), 0, 1))

#format RAND dates with lubridate
hrs_samp %<>% mutate("DOD" = as.Date(hrs_samp$raddate, origin = "1960-01-01"), 
                     "Bday" = as.Date(hrs_samp$rabdate, origin = "1960-01-01"))

# #Sanity check
# View(hrs_samp[, c("Bday", "rabmonth", "rabyear")] %>% na.omit())

#age at death
hrs_samp %<>% 
  mutate("age_death_d" = difftime(DOD, Bday, units = "days"), 
         "age_death_y" = as.numeric(age_death_d/365.25))

# #Sanity check
# View(hrs_samp[, c("Bday", "DOD", "age_death_y")] %>% na.omit())

#---- gender ----
hrs_samp %<>% 
  mutate("female" = ifelse(ragender == 2, 1, 0))

# #sanity check
# table(hrs_samp$female, hrs_samp$ragender)

#---- race-eth ----
#Code any hispanic as 1, else 0
hrs_samp %<>% mutate("hispanic" = ifelse(HISPANIC %in% c(1, 2, 3), 1, 0)) %>% 
  mutate("black" = ifelse(RACE == 2 & hispanic == 0, 1, 0)) %>% 
  mutate("other" = ifelse(RACE == 7 & hispanic == 0, 1, 0)) %>% 
  mutate("unknown_race_eth" = ifelse(RACE == 0 & hispanic == 0, 1, 0))

# #sanity check
# table(hrs_samp$hispanic, hrs_samp$HISPANIC)
# table(hrs_samp$hispanic, hrs_samp$RACE, hrs_samp$black)
# table(hrs_samp$hispanic, hrs_samp$RACE, hrs_samp$other)
# table(hrs_samp$hispanic, hrs_samp$RACE, hrs_samp$unknown_race_eth)

#There are 13 people missing race/ethnicity data so I am dropping them
hrs_samp %<>% filter(unknown_race_eth == 0) %>% 
  #Drop the HRS HISPANIC variable (recoded as hispanic)
  dplyr::select(-one_of("HISPANIC"))

#---- age ----
ages <- analytic_df %>% dplyr::select(contains("AGE")) %>% 
  apply(., 1, impute_ages)
ages[ages > 900] <- NA

analytic_df[, paste0(LETTERS[seq( from = 11, to = 15)], "AGE")] <- t(ages)

#Flag observations with observed ages at least 70yo
analytic_df %<>% 
  mutate("alive_70" = analytic_df %>% dplyr::select(contains("AGE")) %>% 
  apply(., 1, detect_70))

# #sanity Check
# vars <- paste0(LETTERS[seq( from = 11, to = 15)], "AGE")
# 
# for(var in vars){
#   sum(hrs_samp[, var])
#   hist(hrs_samp[, var])
# }

#---- CysC measures ----
#average of all available CysC measures
analytic_df[, "avg_CYSC"] <- 
  analytic_df %>% 
  dplyr::select(contains("CYSC_ADJ")) %>% 
  rowMeans(na.rm = TRUE)
#change NaN to NA for individuals with no measures of CYSC
analytic_df$avg_CYSC[is.nan(analytic_df$avg_CYSC)] <- NA
#z score the measures
analytic_df[, "avg_CYSC_zscore"] <- 
  (analytic_df$avg_CYSC - mean(analytic_df$avg_CYSC, na.rm = TRUE))/
  sd(analytic_df$avg_CYSC, na.rm = TRUE)

#last CysC measure
analytic_df[, "last_CYSC"] <- 
  analytic_df %>% 
  dplyr::select(contains("CYSC_ADJ")) %>% 
  apply(., 1, function(x) non_missing(x, first = FALSE))
#z score the measures
analytic_df[, "last_CYSC_zscore"] <- 
  (analytic_df$last_CYSC - mean(analytic_df$last_CYSC, na.rm = TRUE))/
  sd(analytic_df$last_CYSC, na.rm = TRUE)

#wave of last CysC measure
waves <- paste0(LETTERS[seq( from = 11, to = 15)])
analytic_df[, "last_CYSC_wave"] <-
  analytic_df %>%
  dplyr::select(contains("CYSC_ADJ")) %>%
  apply(., 1, function(x) waves[max(which(!is.na(x)))])

#age at last CysC measure
analytic_df %<>% 
  mutate("last_CYSC_age" = case_when(last_CYSC_wave == "K" ~ KAGE,
                                     last_CYSC_wave == "L" ~ LAGE, 
                                     last_CYSC_wave == "M" ~ MAGE,
                                     last_CYSC_wave == "N" ~ NAGE, 
                                     last_CYSC_wave == "O" ~ OAGE))

# #sanity Check
# View(analytic_df %>% 
#        dplyr::select(c(contains("CYSC_ADJ"), "avg_CYSC", "last_CYSC")))

#---- restrict to those with Cystatin C measures ----
analytic_df %<>% filter(!is.na(last_CYSC_age))

#---- follow-up time ----
analytic_df[, "fu_time"] <- 
  analytic_df %>% 
  dplyr::select(contains(c("CYSC_ADJ", "AGE"))) %>% 
  apply(., 1, fu_time)

analytic_df$fu_time[is.na(analytic_df$fu_time)] <- 0

# #Sanity Check
# View(analytic_df %>% dplyr::select(contains("CYSC_ADJ"), "fu_time"))

#---- save datasets ----
#For mortality between 2014-2016
write_csv(analytic_df, here::here("Data", "analytic_df.csv"))
write_csv(hrs_samp, here::here("Data", "hrs_samp.csv"))

#For 1931-1941 birth cohort
write_csv(analytic_df %>% 
            filter(BIRTHYR %in% seq(1931, 1941, by = 1)), 
          here::here("Data", "analytic_df_1931-1941_cohort.csv"))

#Survival through age 70
write_csv(analytic_df %>% 
            filter(alive_70 == 1), 
          here::here("Data", "analytic_df_alive_70.csv"))


