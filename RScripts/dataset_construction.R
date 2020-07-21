#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "readr", "tidyverse", "magrittr", "plyr", "haven")

#No scientific notation
options(scipen = 999)

#---- source scripts ----
source(here::here("RScripts", "read_da_dct.R"))
source(here::here("RScripts", "non_missing.R"))
source(here::here("RScripts", "fu_time.R"))
source(here::here("RScripts", "impute_ages.R"))

#---- Wave mapping between HRS and RAND ----
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

#---- read in data ----
#2016 tracker file-- for demographics and mortality
#We might want to consider condensing this to only using the RAND file if possible
hrs_tracker <- 
  read_da_dct("/Users/CrystalShaw/Box/HRS/2016_tracker/trk2016/TRK2016TR_R.da", 
              "/Users/CrystalShaw/Box/HRS/2016_tracker/trk2016/TRK2016TR_R.dct", 
              HHIDPN = TRUE)

#2006-2012 biomarker data
years <- c("06", "08", "10", "12")

for(year in years){
  assign(paste0("biomarker_", year), 
         read_da_dct(paste0("/Users/CrystalShaw/Box/HRS/biomarker_data/biomkr", 
                            year, "/BIOMK", year, "BL_R.da"),
                     paste0("/Users/CrystalShaw/Box/HRS/biomarker_data/biomkr", 
                            year, "/BIOMK", year, "BL_R.dct"), HHIDPN = TRUE))
}

#2014 early release biomarker data 
biomarker_14 <- 
  read_da_dct(
    "/Users/CrystalShaw/Box/HRS/biomarker_data/BIOMK14BL/BIOMK14BL.da",
    "/Users/CrystalShaw/Box/HRS/biomarker_data/BIOMK14BL/BIOMK14BL.dct", 
    HHIDPN = TRUE)

#RAND longitudinal file-- for health variables
# RAND <- read_csv(paste0("/Users/CrystalShaw/Box/HRS/RAND_longitudinal/", 
#                         "rndhrs_p.csv")) 
start <- Sys.time()
RAND <- haven::read_sas(paste0("/Users/CrystalShaw/Box/HRS/RAND_longitudinal/", 
                               "randhrs1992_2016v2.sas7bdat"))
end <- Sys.time() - start



#---- pulling variables ----
#We also want their age, sex, race/ethnicity, data to fill in mortality
hrs_samp <- hrs_tracker %>% 
  #participated in HRS 2016 wave (wave "P")
  filter(PIWTYPE %in% c(1, 5, 11, 15)) %>% 
  filter(PALIVE %in% c(1, 5)) %>%
  select("HHIDPN", "PIWTYPE", paste0(LETTERS[seq( from = 11, to = 15)], "AGE"), 
         "GENDER", "RACE", "HISPANIC", "KNOWNDECEASEDMO", "KNOWNDECEASEDYR", 
         "EXDEATHMO", "EXDEATHYR", "PALIVE") %>% 
  #people alive at the beginning of 2014 HRS interview wave
  filter(KNOWNDECEASEDYR >= 2014 | is.na(KNOWNDECEASEDYR)) %>% 
  filter(EXDEATHYR >= 2014 | is.na(EXDEATHYR))
  

#---- DOD ----
#Deriving Date of Death (DOD)
#Looking at values for each variable
table(hrs_samp$KNOWNDECEASEDMO, useNA = "ifany")
table(hrs_samp$KNOWNDECEASEDYR, useNA = "ifany")
table(hrs_samp$EXDEATHYR, useNA = "ifany")
table(hrs_samp$PALIVE, useNA = "ifany")

#Translated from TMM SAS code
hrs_samp %<>% 
  mutate("DOD" = 
           case_when(!(KNOWNDECEASEDMO %in% c(0, 98, NA)) & 
                       !(KNOWNDECEASEDYR %in% c(0, 98, NA)) ~ 
                       paste0(KNOWNDECEASEDMO, "1", KNOWNDECEASEDYR), 
                     TRUE & !is.na(EXDEATHMO) & !is.na(EXDEATHYR) ~ 
                       paste0(EXDEATHMO, "1", EXDEATHYR), 
                     TRUE & PALIVE == 5 ~ "612017")) %>% 
  #died within wave
  mutate("death" = ifelse(is.na(DOD), 0, 1))

# #sanity check 
# #Checking that assignments make sense
# #Checking assignments under first condition: No missing values in this group
# test_1 <- hrs_samp %>% 
#   filter(!(KNOWNDECEASEDMO %in% c(0, 98, NA)) & 
#            !(KNOWNDECEASEDYR %in% c(0, 98, NA))) 
# #There is data in this group
# nrow(test_1) 
# #none of these are empty 
# sum(!is.na(test_1$DOD))
# 
# #Checking assignments under second condition: No data in this group
# test_2 <- hrs_samp %>% 
#   filter(KNOWNDECEASEDMO %in% c(0, 98, NA) & 
#            KNOWNDECEASEDYR %in% c(0, 98, NA)) %>%
#   filter(!is.na(EXDEATHMO) & !is.na(EXDEATHYR)) 
# #There is no data in this group
# nrow(test_2) 
# 
# #Checking assignments under third condition: 
# test_oalive_1 <- hrs_samp %>% 
#   filter(KNOWNDECEASEDMO %in% c(0, 98, NA) & 
#            KNOWNDECEASEDYR %in% c(0, 98, NA)) %>%
#   filter(is.na(EXDEATHMO) | is.na(EXDEATHYR)) %>%
#   filter(OALIVE == 1) 
# #There is data in this group
# nrow(test_oalive_1) 
# #all of these are empty (OALIVE == 1 --> alive at wave)
# sum(!is.na(test_oalive_1$DOD))
# 
# #none of these are empty (OALIVE == 5 --> death reported in wave)
# test_oalive_5 <- hrs_samp %>% 
#   filter(KNOWNDECEASEDMO %in% c(0, 98, NA) & 
#            KNOWNDECEASEDYR %in% c(0, 98, NA)) %>%
#   filter(is.na(EXDEATHMO) | is.na(EXDEATHYR)) %>%
#   filter(OALIVE == 5)
# #There is no data in this group
# nrow(test_oalive_5) 


#---- gender ----
hrs_samp %<>% 
  mutate("female" = ifelse(GENDER == 2, 1, 0))

# #sanity check
# table(hrs_samp$female, hrs_samp$GENDER)

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

#---- merge biomarker data ----
analytic_df <- join_all(list(hrs_samp, biomarker_06, biomarker_08, 
                             biomarker_10, biomarker_12, biomarker_14), 
                        by = "HHIDPN", type = "left") %>% 
  #Select variables of interest
  dplyr::select(colnames(hrs_samp), contains("CYSC"))

#---- age ----
ages <- analytic_df %>% dplyr::select(contains("AGE")) %>% 
  apply(., 1, impute_ages)

analytic_df[, paste0(LETTERS[seq( from = 11, to = 15)], "AGE")] <- t(ages)

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
write_csv(analytic_df, here::here("Data", "analytic_df.csv"))
write_csv(hrs_samp, here::here("Data", "hrs_samp.csv"))

