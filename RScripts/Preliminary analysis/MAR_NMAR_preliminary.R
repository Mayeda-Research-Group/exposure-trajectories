#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "kableExtra")

#---- **Expit function ----
expit <- function(x) {
  output <- (exp(x)/(1+exp(x)))
  return(output)
}

logit <- function(x){
  output <- log(x/(1-x))
  return(output)
}

round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

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

# self-reported health at previous wave
shlt_3_8_long <- CESD_data_wide %>%
  dplyr::select(paste0("r", seq(3, 8, by = 1), "shlt")) %>%
  pivot_longer(everything(), names_to = "orig_varname", 
               values_to = "shlt")

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
e_shlt <- mean(shlt_3_8_long$shlt, na.rm = T)
e_death2018 <- mean(CESD_data_wide$death2018)

#---- ** Coefficients ----
beta_age <- log(1.03)
beta_cesdpre <- log(1.05)
beta_condepre <- log(1.20)
beta_shltpre <- log(1.20)
beta_death2018 <- log(2.00)
beta_cesdcurrent <- log(1.15)

#---- ** beta 0 for MAR ----  
beta_0_MAR_10 <- logit(0.10) - (
  beta_age*e_age + beta_cesdpre*e_CESD_3_8 + beta_condepre*e_conde + 
    beta_shltpre*e_shlt + beta_death2018*e_death2018)
beta_0_MAR_20 <- logit(0.20) - (
  beta_age*e_age + beta_cesdpre*e_CESD_3_8 + beta_condepre*e_conde + 
    beta_shltpre*e_shlt + beta_death2018*e_death2018)
beta_0_MAR_30 <- logit(0.30) - (
  beta_age*e_age + beta_cesdpre*e_CESD_3_8 + beta_condepre*e_conde + 
    beta_shltpre*e_shlt + beta_death2018*e_death2018)

#---- ** beta 0 for MNAR ---- 
beta_0_MNAR_10 <- logit(0.10) - (
  beta_age*e_age + beta_cesdpre*e_CESD_3_8 + beta_condepre*e_conde + 
    beta_shltpre*e_shlt + beta_death2018*e_death2018 + 
    beta_cesdcurrent*e_CESD_4_9)
beta_0_MNAR_20 <- logit(0.20) - (
  beta_age*e_age + beta_cesdpre*e_CESD_3_8 + beta_condepre*e_conde + 
    beta_shltpre*e_shlt + beta_death2018*e_death2018 + 
    beta_cesdcurrent*e_CESD_4_9)
beta_0_MNAR_30 <- logit(0.30) - (
  beta_age*e_age + beta_cesdpre*e_CESD_3_8 + beta_condepre*e_conde + 
    beta_shltpre*e_shlt + beta_death2018*e_death2018 +
    beta_cesdcurrent*e_CESD_4_9)

#----  MAR and MNAR function ----
MAR_MNAR_func <- function (beta_0_MAR, beta_0_MNAR, dataset){
  
  n <- nrow(dataset)
  subset <- dataset %>%
    mutate( 
      #MAR
      r4pcesd_MAR = expit(beta_0_MAR + beta_age * r4age_y_int + 
                            beta_cesdpre * r3cesd + 
                            beta_condepre * r3conde_impute + 
                            beta_shltpre * r3shlt + beta_death2018 * death2018),
      r5pcesd_MAR = expit(beta_0_MAR + beta_age * r5age_y_int + 
                            beta_cesdpre * r4cesd + 
                            beta_condepre * r4conde_impute + 
                            beta_shltpre * r4shlt + beta_death2018 * death2018),
      r6pcesd_MAR = expit(beta_0_MAR + beta_age * r6age_y_int + 
                            beta_cesdpre * r5cesd + 
                            beta_condepre * r5conde_impute + 
                            beta_shltpre * r5shlt + beta_death2018 * death2018),
      r7pcesd_MAR = expit(beta_0_MAR + beta_age * r7age_y_int + 
                            beta_cesdpre * r6cesd + 
                            beta_condepre * r6conde_impute + 
                            beta_shltpre * r6shlt + beta_death2018 * death2018),
      r8pcesd_MAR = expit(beta_0_MAR + beta_age * r8age_y_int + 
                            beta_cesdpre * r7cesd + 
                            beta_condepre * r7conde_impute + 
                            beta_shltpre * r7shlt + beta_death2018 * death2018),
      r9pcesd_MAR = expit(beta_0_MAR + beta_age * r9age_y_int + 
                            beta_cesdpre * r8cesd + 
                            beta_condepre * r8conde_impute + 
                            beta_shltpre * r8shlt + beta_death2018 * death2018),
      # MNAR
      r4pcesd_MNAR = expit(beta_0_MNAR + beta_age * r4age_y_int + 
                             beta_cesdpre * r3cesd + beta_cesdcurrent * r4cesd +
                             beta_condepre * r3conde_impute + 
                             beta_shltpre * r3shlt + beta_death2018 * death2018),
      r5pcesd_MNAR = expit(beta_0_MNAR + beta_age * r5age_y_int + 
                             beta_cesdpre * r4cesd + beta_cesdcurrent * r5cesd +
                             beta_condepre * r4conde_impute + 
                             beta_shltpre * r4shlt + beta_death2018 * death2018),
      r6pcesd_MNAR = expit(beta_0_MNAR + beta_age * r6age_y_int + 
                             beta_cesdpre * r5cesd + beta_cesdcurrent * r6cesd +
                             beta_condepre * r5conde_impute + 
                             beta_shltpre * r5shlt + beta_death2018 * death2018),
      r7pcesd_MNAR = expit(beta_0_MNAR + beta_age * r7age_y_int + 
                             beta_cesdpre * r6cesd + beta_cesdcurrent * r7cesd +
                             beta_condepre * r6conde_impute + 
                             beta_shltpre * r6shlt + beta_death2018 * death2018),
      r8pcesd_MNAR = expit(beta_0_MNAR + beta_age * r8age_y_int + 
                             beta_cesdpre * r7cesd + beta_cesdcurrent * r8cesd +
                             beta_condepre * r7conde_impute + 
                             beta_shltpre * r7shlt + beta_death2018 * death2018),
      r9pcesd_MNAR = expit(beta_0_MNAR + beta_age * r9age_y_int + 
                             beta_cesdpre * r8cesd + beta_cesdcurrent * r9cesd +
                             beta_condepre * r8conde_impute + 
                             beta_shltpre * r8shlt + beta_death2018 * death2018)
    ) %>%
    select(contains(c("MAR", "MNAR")))
  
  # Flag the score based on bernoulli distribution, prob = p_wave_MAR/MNAR
  for (j in 1:ncol(subset)){
    if (j <= 6){
      subset[, paste0("r", j + 3, "cesd_missing_MAR")] <- 
        rbinom(nrow(subset), size = 1, prob = subset[[j]])
    } 
    else{
      subset[, paste0("r", j - 3, "cesd_missing_MNAR")] <-
        rbinom(nrow(subset), size = 1, prob = subset[[j]])
    }
  }
  
  # Calculate the missing proportion for MAR and MNAR
  MAR_long <- subset %>%
    select(contains("cesd_missing_MAR")) %>%
    pivot_longer(
      everything(),
      names_to = "Variables",
      values_to = "Missingness"
    )
  
  MNAR_long <- subset %>%
    select(contains("cesd_missing_MNAR")) %>%
    pivot_longer(
      everything(),
      names_to = "Variables",
      values_to = "Missingness"
    )
  
  return(
    Missing_prop_results <- tibble(
      MAR_missing_prop = mean(MAR_long$Missingness, na.rm = T),
      MNAR_missing_prop = mean(MNAR_long$Missingness, na.rm = T))
  )
}

# Repeat for 1000 times
replicate = 1000
set.seed(20210226)

#---- Missing proportion results ----
system.time({
  # 10%
  missing_prop_represults_10 <- 
    map_dfr(1:replicate, ~ MAR_MNAR_func(beta_0_MAR_10, beta_0_MNAR_10, 
                                         CESD_data_wide),
            .id = "replication")
  
  #25%
  missing_prop_represults_25 <- 
    map_dfr(1:replicate, ~ MAR_MNAR_func(beta_0_MAR_25, beta_0_MNAR_25,
                                         CESD_data_wide),
            .id = "replication")
  
  #50%
  missing_prop_represults_50 <-
    map_dfr(1:replicate, ~ MAR_MNAR_func(beta_0_MAR_50, beta_0_MNAR_50,
                                         CESD_data_wide),
            .id = "replication")
})

#---- ** tibble results ----
#MAR
{
  missing_prop_MAR_summary <- tibble(
    "Prop_Missingness" = c("10%","25%", "50%"),
    "Mean" = c(mean(missing_prop_represults_10$MAR_missing_prop),
               mean(missing_prop_represults_25$MAR_missing_prop),
               mean(missing_prop_represults_50$MAR_missing_prop)),
    "p2.5th" = c(quantile(missing_prop_represults_10$MAR_missing_prop, 0.025),
                 quantile(missing_prop_represults_25$MAR_missing_prop, 0.025),
                 quantile(missing_prop_represults_50$MAR_missing_prop, 0.025)),
    "p97.5th" = c(quantile(missing_prop_represults_10$MAR_missing_prop, 0.975),
                  quantile(missing_prop_represults_25$MAR_missing_prop, 0.975),
                  quantile(missing_prop_represults_50$MAR_missing_prop, 0.975))
  ) %>%
    round_df(digits = 4)
  
  missing_prop_MAR_summary %>%
    kbl(caption = "MAR missingness") %>%
    kable_classic(full_width = F, html_font = "Arial")
}

# MNAR
{
  missing_prop_MNAR_summary <- tibble(
    "Prop_Missingness" = c("10%","25%", "50%"),
    "Mean" = c(mean(missing_prop_represults_10$MNAR_missing_prop),
               mean(missing_prop_represults_25$MNAR_missing_prop),
               mean(missing_prop_represults_50$MNAR_missing_prop)),
    "p2.5th" = c(quantile(missing_prop_represults_10$MNAR_missing_prop, 0.025),
                 quantile(missing_prop_represults_25$MNAR_missing_prop, 0.025),
                 quantile(missing_prop_represults_50$MNAR_missing_prop, 0.025)),
    "p97.5th" = c(quantile(missing_prop_represults_10$MNAR_missing_prop, 0.975),
                  quantile(missing_prop_represults_25$MNAR_missing_prop, 0.975),
                  quantile(missing_prop_represults_50$MNAR_missing_prop, 0.975))
  ) %>%
    round_df(digits = 4)
  
  missing_prop_MNAR_summary %>%
    kbl(caption = "MNAR missingness") %>%
    kable_classic(full_width = F, html_font = "Arial")
}

# Test
# p <- 0.5; (p_new <- 1/((1/(p/(1-p))^2) + 1))


# #10%
# summary(missing_prop_represults_10$MAR_missing_prop)
# summary(missing_prop_represults_10$MNAR_missing_prop)
# #25%
# summary(missing_prop_represults_25$MAR_missing_prop)
# summary(missing_prop_represults_25$MNAR_missing_prop)
# #50%
# summary(missing_prop_represults_50$MAR_missing_prop)
# summary(missing_prop_represults_50$MNAR_missing_prop)


# system.time(
#   missing_prop_represults_10 <-
#     map_dfr(1:replicate, ~MAR_MNAR_func(beta_0_MAR_10, CESD_data_wide),
#             .id = "replication")
# )

