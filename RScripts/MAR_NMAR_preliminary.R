#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr")

#---- **Expit function ----
expit <- function(x) {
  output <- (exp(x)/(1+exp(x)))
  return(output)
}

logit <- function(x){
  output <- log(x/(1-x))
  return(output)
}

#---- note ----
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
                  "CESD_data_wide.csv"))

#---- Getting the balancing intercept ----

# colnames(CESD_data_wide)
# Age
overall_ages <- CESD_data_wide %>% 
  dplyr::select(paste0("r", seq(4, 9, by = 1), "age_y_int")) %>% 
  pivot_longer(everything(), names_to = "orig_varname", 
               values_to = "age_int")

# CESD previous wave (wave 3-8)
CESD_3_8_long <- CESD_data_wide %>%
  dplyr::select(paste0("r", seq(3, 8, by = 1), "cesd")) %>%
  pivot_longer(everything(), names_to = "orig_varname", 
               values_to = "cesd")

# conde at previous wave
conde_3_8_long <- CESD_data_wide %>%
  dplyr::select(paste0("r", seq(3, 8, by = 1), "conde_impute")) %>%
  pivot_longer(everything(), names_to = "orig_varname", 
               values_to = "conde_impute")

# CESD this wave (wave 4-9)
CESD_4_9_long <- CESD_data_wide %>%
  dplyr::select(paste0("r", seq(4, 9, by = 1), "cesd")) %>%
  pivot_longer(everything(), names_to = "orig_varname", 
               values_to = "cesd")

#---- **E(X)s ----
e_age <- mean(overall_ages$age_int)
e_CESD_3_8 <- mean(CESD_3_8_long$cesd, na.rm = T)
e_CESD_4_9 <- mean(CESD_4_9_long$cesd)
e_conde <- mean(conde_3_8_long$conde_impute, na.rm = T)

#---- ** beta 0 for MAR ----  
beta_0_MAR_10 <- logit(0.1) - (
  log(1.05)*e_age + log(1.10)*e_CESD_3_8 + log(1.25)*e_conde)
beta_0_MAR_25 <- logit(0.25) - (
  log(1.05)*e_age + log(1.10)*e_CESD_3_8 + log(1.25)*e_conde)
beta_0_MAR_50 <- logit(0.50) - (
  log(1.05)*e_age + log(1.10)*e_CESD_3_8 + log(1.25)*e_conde)

#---- ** beta 0 for NMAR ---- 
beta_0_NMAR_10 <- logit(0.1) - (
  log(1.05)*e_age + log(1.10)*e_CESD_3_8 + log(1.25)*e_conde + 
    log(1.15)*e_CESD_4_9)
beta_0_NMAR_25 <- logit(0.25) - (
  log(1.05)*e_age + log(1.10)*e_CESD_3_8 + log(1.25)*e_conde + 
    log(1.15)*e_CESD_4_9)
beta_0_NMAR_50 <- logit(0.50) - (
  log(1.05)*e_age + log(1.10)*e_CESD_3_8 + log(1.25)*e_conde +
    log(1.15)*e_CESD_4_9)

#----  MAR and NMAR function ----
# TBC!!!
MAR_NMAR_func <- function (beta_0, dataset){ # x needs to be changed
  
  n <- nrow(dataset)
  subset <- dataset %>%
    mutate( 
      #MAR
      r4pcesd_MAR = expit(beta_0 + log(1.05) * r4age_y_int + 
                            log(1.10) * r3cesd
                          + log(1.25) * r4conde_impute),
      r5pcesd_MAR = expit(beta_0 + log(1.05) * r5age_y_int + 
                            log(1.10) * r4cesd
                          + log(1.25) * r5conde_impute),
      r6pcesd_MAR = expit(beta_0 + log(1.05) * r6age_y_int + 
                            log(1.10) * r5cesd
                          + log(1.25) * r6conde_impute),
      r7pcesd_MAR = expit(beta_0 + log(1.05) * r7age_y_int + 
                            log(1.10) * r6cesd
                          + log(1.25) * r7conde_impute),
      r8pcesd_MAR = expit(beta_0 + log(1.05) * r8age_y_int + 
                            log(1.10) * r7cesd
                          + log(1.25) * r8conde_impute),
      r9pcesd_MAR = expit(beta_0 + log(1.05) * r9age_y_int + 
                            log(1.10) * r8cesd
                          + log(1.25) * r9conde_impute),
      # NMAR
      r4pcesd_NMAR = expit(beta_0 + log(1.05) * r4age_y_int + 
                             log(1.10) * r3cesd + log(1.15) * r4cesd +
                           + log(1.25) * r4conde_impute),
      r5pcesd_NMAR = expit(beta_0 + log(1.05) * r5age_y_int + 
                             log(1.10) * r4cesd + log(1.15) * r5cesd +
                           + log(1.25) * r5conde_impute),
      r6pcesd_NMAR = expit(beta_0 + log(1.05) * r6age_y_int + 
                             log(1.10) * r5cesd + log(1.15) * r6cesd +
                           + log(1.25) * r6conde_impute),
      r7pcesd_NMAR = expit(beta_0 + log(1.05) * r7age_y_int + 
                             log(1.10) * r6cesd + log(1.15) * r7cesd +
                           + log(1.25) * r7conde_impute),
      r8pcesd_NMAR = expit(beta_0 + log(1.05) * r8age_y_int + 
                             log(1.10) * r7cesd + log(1.15) * r8cesd +
                           + log(1.25) * r8conde_impute),
      r9pcesd_NMAR = expit(beta_0 + log(1.05) * r9age_y_int + 
                             log(1.10) * r8cesd + log(1.15) * r9cesd +
                           + log(1.25) * r9conde_impute)
    ) %>%
    select(contains("MAR"))
  
  # Flag the score based on bernoulli distribution, prob = p_wave_MAR/NMAR
  for (j in 1:ncol(subset)){
    if (j <= 6){
      subset[, paste0("r", j + 3, "cesd_missing_MAR")] <- 
        rbinom(nrow(subset), size = 1, prob = subset[[j]])
    } 
    else{
      subset[, paste0("r", j - 3, "cesd_missing_NMAR")] <-
        rbinom(nrow(subset), size = 1, prob = subset[[j]])
    }
  }

  # Calculate the missing proportion for MAR and NMAR
MAR_long <- subset %>%
    select(contains("cesd_missing_MAR")) %>%
    pivot_longer(
      everything(),
      names_to = "Variables",
      values_to = "Missingness"
    )

NMAR_long <- subset %>%
  select(contains("cesd_missing_NMAR")) %>%
  pivot_longer(
    everything(),
    names_to = "Variables",
    values_to = "Missingness"
  )

return(
  Missing_prop_results <- tibble(
    MAR_missing_prop = round(mean(MAR_long$Missingness, na.rm = T), 4),
    NMAR_missing_prop = round(mean(NMAR_long$Missingness, na.rm = T), 4)
  )
)
}

# Bootstrap for 1000 times
bootsize = 1000
set.seed(123)
missing_prop_bootresults <- 
  map_dfr(1:bootsize, ~MAR_NMAR_func(beta_0_MAR_10, CESD_data_wide),
          .id = "replication")
  
summary(missing_prop_bootresults$MAR_missing_prop)
summary(missing_prop_bootresults$NMAR_missing_prop)
