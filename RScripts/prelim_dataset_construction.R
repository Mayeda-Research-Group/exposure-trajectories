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
source(here::here("RScripts", "detect_70.R"))
source(here::here("RScripts", "cysc_before_70.R"))

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

years <- c("06", "08", "10", "12", "14") #biomarker sample
letter_waves <- LETTERS[seq(from = 11, to = 16)] #biomarker sample + 2016 HRS
number_waves <- seq(8, 13, by = 1) #biomarker sample + 2016 HRS

#---- read in data ----
hrs_tracker <- 
  read_da_dct("/Users/CrystalShaw/Box/HRS/2016_tracker/trk2016/TRK2016TR_R.da", 
              "/Users/CrystalShaw/Box/HRS/2016_tracker/trk2016/TRK2016TR_R.dct", 
              HHIDPN = TRUE) %>% 
  select("HHIDPN", "OIWTYPE", "OALIVE", "PIWTYPE", "PALIVE") %>% 
  mutate_at("HHIDPN", as.numeric)

#2006-2012 biomarker data and core data 
dataframes_list <- vector(mode = "list", length = 2*length(years))

for(i in 1:length(years)){
  year <- years[i]
  
  if(i != length(years)){
    dataframes_list[[i]] <-  
      read_da_dct(paste0("/Users/CrystalShaw/Box/HRS/biomarker_data/biomkr", 
                         year, "/BIOMK", year, "BL_R.da"),
                  paste0("/Users/CrystalShaw/Box/HRS/biomarker_data/biomkr", 
                         year, "/BIOMK", year, "BL_R.dct"), HHIDPN = TRUE) %>%
      mutate_at("HHIDPN", as.numeric) %>% 
      #Select variables of interest: ID, adjusted Cystatin C, adjusted HbA1c,
      #                              adjusted total cholesterol, adjusted HDL
      dplyr::select(HHIDPN, contains("CYSC_ADJ"), contains("A1C_ADJ"), 
                    contains("TC_ADJ"), contains("HDL_ADJ"))
  } else{
    #2014 early release biomarker data
    dataframes_list[[i]] <-  read_da_dct(
      "/Users/CrystalShaw/Box/HRS/biomarker_data/BIOMK14BL/BIOMK14BL.da",
      "/Users/CrystalShaw/Box/HRS/biomarker_data/BIOMK14BL/BIOMK14BL.dct", 
      HHIDPN = TRUE) %>%
      mutate_at("HHIDPN", as.numeric) %>% 
      #Select variables of interest: ID, adjusted Cystatin C, adjusted HbA1c,
      #                              adjusted total cholesterol, adjusted HDL
      dplyr::select(HHIDPN, contains("CYSC_ADJ"), contains("A1C_ADJ"), 
                    contains("TC_ADJ"), contains("HDL_ADJ"))
    
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

RAND <- read_dta(paste0("/Users/CrystalShaw/Box/HRS/RAND_longitudinal/STATA/", 
                        "randhrs1992_2016v2.dta"), 
                 col_select = all_of(rand_variables)) 

colnames(RAND)[1] <- "HHIDPN" #For merging

#Remove labeled data format
val_labels(RAND) <- NULL

#HRS Core files-- Need for physical measures
#count continues from biomarker data pull
for(i in (length(years) + 1):length(dataframes_list)){
  year <- years[i - length(years)]
  dataframes_list[[i]] <- 
    read_da_dct(paste0("/Users/CrystalShaw/Box/HRS/core_files/h", year, 
                       "core/h", year, "da/H", year, "I_R.da"),
                paste0("/Users/CrystalShaw/Box/HRS/core_files/h", year, 
                       "core/h", year, "sta/H", year, "I_R.dct"), 
                HHIDPN = TRUE) %>% mutate_at("HHIDPN", as.numeric) %>% 
    #Select variables of interest: ID, Weight (pounds), Height (inches), 
    #                              systolic bp (3 times), 
    #                              diastolic bp (3 times)
    dplyr::select(HHIDPN, contains("I841"), contains("I834"), 
                  contains("I859"), contains("I864"), contains("I869"), 
                  contains("I860"), contains("I865"), contains("I870")) %>% 
    set_colnames(c("HHIDPN", paste0(letter_waves[i - length(years)], "wt"), 
                   paste0(letter_waves[i - length(years)], "ht"), 
                   paste0(letter_waves[i - length(years)], "sbp", seq(1, 3)), 
                   paste0(letter_waves[i - length(years)], "dbp", seq(1, 3))))
}

#Anusha Vable's CSES index
cSES <- read_dta(paste0("~/Dropbox/Projects/exposure_trajectories/data/", 
                        "cSES measures/cses_measures.dta")) %>% 
  unite(col = "HHIDPN", c("hhid", "pn"), sep = "") %>% 
  mutate_at("HHIDPN", as.numeric)

#---- merge datasets ----
#Use this to subset RAND data
hrs_samp <- join_all(c(list(hrs_tracker, RAND, cSES), dataframes_list), 
                     by = "HHIDPN", type = "left") 

#---- at least one CysC measure ----
hrs_samp %<>% 
  mutate("some_cysc" = hrs_samp %>% dplyr::select(contains("CYSC_ADJ")) %>% 
  apply(1, function(x) sum(1 - is.na(x)))) %>% filter(some_cysc != 0)

#---- death ----
#death indicator
hrs_samp %<>% mutate("death" = ifelse(is.na(raddate), 0, 1))

#format RAND dates with lubridate
hrs_samp %<>% mutate("DOD" = as.Date(hrs_samp$raddate, origin = "1960-01-01"), 
                     "Bday" = as.Date(hrs_samp$rabdate, origin = "1960-01-01"))

# #Sanity check
# View(hrs_samp[, c("Bday", "rabmonth", "rabyear")] %>% na.omit())
# View(hrs_samp[, c("DOD", "radmonth", "radyear")] %>% na.omit())

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
# table(hrs_samp$female, hrs_samp$ragender, useNA = "ifany")

#---- race-eth ----
#Code any hispanic as 1, else 0
hrs_samp %<>% 
  mutate("hispanic" = ifelse(rahispan == 0 | is.na(rahispan), 0, 1)) %>% 
  mutate("white" = ifelse(raracem == 1 & hispanic == 0, 1, 0)) %>%
  mutate("black" = ifelse(raracem == 2 & hispanic == 0, 1, 0)) %>% 
  mutate("other" = ifelse(raracem == 3 & hispanic == 0, 1, 0)) %>% 
  mutate("unknown_race_eth" = ifelse(is.na(raracem) & hispanic == 0, 1, 0))

# #sanity check
# table(hrs_samp$hispanic, hrs_samp$rahispan, useNA = "ifany")
# table(hrs_samp$hispanic, hrs_samp$raracem, hrs_samp$black, useNA = "ifany")
# table(hrs_samp$hispanic, hrs_samp$raracem, hrs_samp$other, useNA = "ifany")
# table(hrs_samp$hispanic, hrs_samp$raracem, hrs_samp$unknown_race_eth,
#       useNA = "ifany")
# table(hrs_samp$unknown_race_eth, useNA = "ifany")

#There are 6 people missing race/ethnicity data so I am dropping them
hrs_samp %<>% filter(unknown_race_eth == 0) %>% 
  #Drop the RAND rahispan variable (recoded as hispanic)
  dplyr::select(-one_of("rahispan"))

#---- age ----
age_m <- hrs_samp %>% dplyr::select(contains("agem_e")) %>% 
  apply(., 1, impute_ages) %>% t() 

hrs_samp[, paste0(letter_waves, "age_y")] <- age_m/12

#Flag observations with observed ages at least 70yo
hrs_samp %<>% 
  mutate("alive_70" = hrs_samp %>% 
           dplyr::select(paste0(head(letter_waves, -1), "age_y")) %>% 
  apply(., 1, detect_70))

#Flag CysC measure before 70
hrs_samp[, "cysc_before_70"] <- 
  hrs_samp %>% 
  dplyr::select(contains(c("CYSC_ADJ", "age_y"))) %>% 
  apply(., 1, cysc_before_70)

#Flag CysC measures between ages [60-70)
hrs_samp[, "cysc_between_60_70"] <- 
  hrs_samp %>% 
  dplyr::select(contains(c("CYSC_ADJ", "age_y"))) %>% 
  apply(., 1, cysc_between_60_70)

#---- CysC measures ----
#average of all available CysC measures
hrs_samp[, "avg_CYSC"] <- 
  hrs_samp %>% 
  dplyr::select(contains("CYSC_ADJ")) %>% 
  rowMeans(na.rm = TRUE)
#change NaN to NA for individuals with no measures of CYSC
hrs_samp$avg_CYSC[is.nan(hrs_samp$avg_CYSC)] <- NA
#z score the measures
hrs_samp[, "avg_CYSC_zscore"] <- 
  (hrs_samp$avg_CYSC - mean(hrs_samp$avg_CYSC, na.rm = TRUE))/
  sd(hrs_samp$avg_CYSC, na.rm = TRUE)

#last CysC measure
hrs_samp[, "last_CYSC"] <- 
  hrs_samp %>% 
  dplyr::select(contains("CYSC_ADJ")) %>% 
  apply(., 1, function(x) non_missing(x, first = FALSE))

#z score the measures
hrs_samp[, "last_CYSC_zscore"] <- 
  (hrs_samp$last_CYSC - mean(hrs_samp$last_CYSC, na.rm = TRUE))/
  sd(hrs_samp$last_CYSC, na.rm = TRUE)

#wave of first CysC measure
hrs_samp[, "first_CYSC_wave"] <-
  hrs_samp %>%
  dplyr::select(contains("CYSC_ADJ")) %>%
  apply(., 1, function(x) letter_waves[min(which(!is.na(x)))])

#wave of last CysC measure
hrs_samp[, "last_CYSC_wave"] <-
  hrs_samp %>%
  dplyr::select(contains("CYSC_ADJ")) %>%
  apply(., 1, function(x) letter_waves[max(which(!is.na(x)))])

#age at first CysC measure
#restrict sample to those age-eligible at their first Cystatin C measure
hrs_samp %<>% 
  mutate("first_CYSC_age" = case_when(first_CYSC_wave == "K" ~ Kage_y,
                                      first_CYSC_wave == "L" ~ Lage_y, 
                                      first_CYSC_wave == "M" ~ Mage_y,
                                      first_CYSC_wave == "N" ~ Nage_y, 
                                      first_CYSC_wave == "O" ~ Oage_y)) %>% 
  filter(first_CYSC_age >= 50)

#age at last CysC measure
hrs_samp %<>% 
  mutate("last_CYSC_age" = case_when(last_CYSC_wave == "K" ~ Kage_y,
                                     last_CYSC_wave == "L" ~ Lage_y, 
                                     last_CYSC_wave == "M" ~ Mage_y,
                                     last_CYSC_wave == "N" ~ Nage_y, 
                                     last_CYSC_wave == "O" ~ Oage_y))

# #sanity Check
# View(hrs_samp %>% 
#        dplyr::select(c(contains("CYSC_ADJ"), "avg_CYSC", "last_CYSC")))

#---- double check restricting to those with Cystatin C measures ----
# #There are 2 people with missing age data at every wave-- will fix this
# hrs_samp %<>% filter(!is.na(last_CYSC_age))

#---- follow-up time ----
hrs_samp[, "fu_time"] <- 
  hrs_samp %>% 
  dplyr::select(contains(c("CYSC_ADJ", "age_y"))) %>% 
  apply(., 1, fu_time)

hrs_samp$fu_time[is.na(hrs_samp$fu_time)] <- 0

# #Sanity Check
# View(hrs_samp %>% dplyr::select(contains("CYSC_ADJ"), "fu_time"))

#---- checking analytical sample ----
analytical_sample <- hrs_samp %>% 
  filter(alive_70 == 1) %>% 
  filter(cysc_between_60_70 == 1)

nrow(analytical_sample)
table(analytical_sample$alive_70, useNA = "ifany")
table(analytical_sample$cysc_between_60_70, useNA = "ifany")
table(hrs_samp$alive_70, hrs_samp$cysc_between_60_70, 
      useNA = "ifany")

#Do we have cSES measures for everyone?-- missing 19 people
View(analytical_sample[which(is.na(analytical_sample$cses_index)), ])

#---- save datasets ----
#Participation in 2016 HRS (core or exit: mortality between 2014-2016 only)
write_csv(hrs_samp %>% 
            filter(PIWTYPE %in% c(1, 5, 11, 15)) %>% 
            filter(PALIVE %in% c(1, 5)), here::here("Data", "hrs_samp.csv"))

#For 1931-1941 birth cohort
write_csv(hrs_samp %>% 
            filter(PIWTYPE %in% c(1, 5, 11, 15)) %>% 
            filter(PALIVE %in% c(1, 5)) %>%
            filter(rabyear %in% seq(1931, 1941, by = 1)), 
          here::here("Data", "hrs_samp_1931-1941_cohort.csv"))

#Survival through age 70 and at least one cystatin c measure in [60, 70)
write_csv(hrs_samp %>% 
            filter(alive_70 == 1) %>% 
            filter(cysc_between_60_70 == 1),
          here::here("Data", "hrs_samp_alive_70.csv"))


