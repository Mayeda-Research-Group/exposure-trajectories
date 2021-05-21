# Aim: to simulate many missing datasets to ensure that we are achieving our 
#   desired levels of missingness
# Created by: Yingyan Wu
# 05.07.2021
# Edited by: Crystal Shaw
# 05.21.2021

#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "kableExtra")

#---- functions ----
expit <- function(x) {
  output <- (exp(x)/(1+exp(x)))
  return(output)
}

logit <- function(x){
  output <- log(x/(1-x))
  return(output)
}

estimate_df <- function(dataframe){
  estimate <- function(x){
    stats <- c(length(x), 
               mean(x,na.rm = T), 
               quantile(x,c(0.025, 0.975),na.rm = T)
    )
    names(stats) <- c("N", "Mean", 
                      "p2.5th", "p97.5th"
    )
    return(round(stats, 4))
  }
  
  stats_tib <- dataframe %>%
    select(where(is.numeric)) %>%
    purrr::map_dfr(estimate) %>%
    as_tibble() %>%
    mutate("Effect" = colnames(dataframe)[-1]) %>%
    select(Effect, everything())
  
  return(stats_tib)
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
data_wide <- 
  read_csv(paste0(path_to_dropbox, 
                  "/exposure_trajectories/data/", 
                  "CESD_data_wide.csv"))

#---- expected value for predictors of MAR & MNAR missingness ----
#---- **E(X)s ----
e_age <- 
  mean(unlist(data_wide[, paste0("r", seq(4, 9, by = 1), "age_y_int")]))
e_CESD_3_8 <- 
  mean(unlist(data_wide[, paste0("r", seq(3, 8, by = 1), "cesd")]), 
       na.rm = TRUE)
e_CESD_4_9 <- 
  mean(unlist(data_wide[, paste0("r", seq(4, 9, by = 1), "cesd")]))
e_conde <- 
  mean(unlist(data_wide[, paste0("r", seq(3, 8, by = 1), "conde_impute")]), 
       na.rm = TRUE)
e_shlt <- mean(unlist(data_wide[, paste0("r", seq(3, 8, by = 1), "shlt")]), 
               na.rm = TRUE)
e_death2018 <- mean(data_wide$death2018)

# #---- ** MAR Coefficients ----
# #beta_age <- log(0.97*scale)
# beta_cesdpre <- log(1.04*scale)
# beta_condepre <- log(1.30*scale)
# #beta_shltpre <- log(1.100)
# beta_death2018 <- log(3*scale)

#---- ** MNAR Coefficients ----
#beta_age <- log(0.955*scale)
#beta_cesdpre <- log(1.04*scale)
#beta_condepre <- log(1.30*scale)
#beta_shltpre <- log(1.100)
beta_death2018 <- log(1.5)
beta_cesdcurrent <- log(4)
beta_death2018*cesdcurrent <- log(1.5)

#----  scale coef function ----
coef_func <- function (dataset, mechanism){
  
  n <- nrow(dataset)
  beta_0 <- -6
  #beta_age_s <- beta_age * scale
  beta_cesdpre_s <- beta_cesdpre * scale
  beta_condepre_s <- beta_condepre * scale
  #beta_shltpre_s <- beta_shltpre * scale
  beta_death2018_s <- beta_death2018 * scale
  beta_cesdcurrent_s <- beta_cesdcurrent * scale
  
  #  # balancing intercept (TEMP)
  #  #---- **E(X)s ----
  #  e_age <- 
  #    mean(unlist(dataset[, paste0("r", seq(4, 9, by = 1), "age_y_int")]))
  #  e_CESD_3_8 <- 
  #    mean(unlist(dataset[, paste0("r", seq(3, 8, by = 1), "cesd")]), 
  #         na.rm = TRUE)
  #  e_CESD_4_9 <- 
  #    mean(unlist(dataset[, paste0("r", seq(4, 9, by = 1), "cesd")]))
  #  e_conde <- 
  #    mean(unlist(dataset[, paste0("r", seq(3, 8, by = 1), "conde_impute")]), 
  #         na.rm = TRUE)
  #  e_shlt <- mean(unlist(dataset[, paste0("r", seq(3, 8, by = 1), "shlt")]), 
  #                 na.rm = TRUE)
  #  e_death2018 <- mean(dataset$death2018)
  #  
  # ( beta_0_MAR <- logit(0.40) - (
  #    beta_age_s*e_age + beta_cesdpre_s*e_CESD_3_8 + beta_condepre_s*e_conde + 
  #      beta_shltpre_s*e_shlt + beta_death2018_s*e_death2018))
  #  exp(beta_0_MAR_10)
  
  if (mechanism == "MAR"){
    #---- MAR ----
    subset <- dataset %>%
      mutate( 
        r4pcesd_MAR = expit(#beta_age_s * r4age_y_int + 
          beta_cesdpre_s * r3cesd +
            beta_condepre_s * r3conde_impute +
            #beta_shltpre_s * r3shlt + 
            beta_death2018_s * death2018 + 
            beta_0),
        r5pcesd_MAR = expit(#beta_age_s * r5age_y_int + 
          beta_cesdpre_s * r4cesd +
            beta_condepre_s * r4conde_impute +
            #beta_shltpre_s * r4shlt + 
            beta_death2018_s * death2018 + 
            beta_0),
        r6pcesd_MAR = expit(#beta_age_s * r6age_y_int + 
          beta_cesdpre_s * r5cesd +
            beta_condepre_s * r5conde_impute +
            #beta_shltpre_s * r5shlt + 
            beta_death2018_s * death2018 + 
            beta_0),
        r7pcesd_MAR = expit(#beta_age_s * r7age_y_int + 
          beta_cesdpre_s * r6cesd +
            beta_condepre_s * r6conde_impute +
            #beta_shltpre_s * r6shlt + 
            beta_death2018_s * death2018 + 
            beta_0),
        r8pcesd_MAR = expit(#beta_age_s * r8age_y_int + 
          beta_cesdpre_s * r7cesd +
            beta_condepre_s * r7conde_impute +
            #beta_shltpre_s * r7shlt + 
            beta_death2018_s * death2018 + 
            beta_0),
        r9pcesd_MAR = expit(#beta_age_s * r9age_y_int + 
          beta_cesdpre_s * r8cesd +
            beta_condepre_s * r8conde_impute +
            #beta_shltpre_s * r8shlt + 
            beta_death2018_s * death2018 + 
            beta_0)) %>%
      dplyr::select(contains("MAR", ignore.case = FALSE))
  } else if(mechanism == "MNAR"){
    #---- MNAR ----
    subset <- dataset %>%
      mutate( 
        r4pcesd_MNAR = expit(#beta_age_s * r4age_y_int + 
          beta_cesdpre_s * r3cesd + 
            beta_cesdcurrent_s * r4cesd +
            beta_condepre_s * r3conde_impute + 
            #beta_shltpre_s * r3shlt + 
            beta_death2018_s * death2018 + beta_0),
        r5pcesd_MNAR = expit(#beta_age_s * r5age_y_int + 
          beta_cesdpre_s * r4cesd + 
            beta_cesdcurrent_s * r5cesd +
            beta_condepre_s * r4conde_impute + 
            #beta_shltpre_s * r4shlt + 
            beta_death2018_s * death2018 + beta_0),
        r6pcesd_MNAR = expit(#beta_age_s * r6age_y_int + 
          beta_cesdpre_s * r5cesd + 
            beta_cesdcurrent_s * r6cesd +
            beta_condepre_s * r5conde_impute + 
            #beta_shltpre_s * r5shlt + 
            beta_death2018_s * death2018 + beta_0),
        r7pcesd_MNAR = expit(#beta_age_s * r7age_y_int + 
          beta_cesdpre_s * r6cesd + 
            beta_cesdcurrent_s * r7cesd +
            beta_condepre_s * r6conde_impute + 
            #beta_shltpre_s * r6shlt + 
            beta_death2018_s * death2018 + beta_0),
        r8pcesd_MNAR = expit(#beta_age_s * r8age_y_int + 
          beta_cesdpre_s * r7cesd + 
            beta_cesdcurrent_s * r8cesd +
            beta_condepre_s * r7conde_impute + 
            #beta_shltpre_s * r7shlt + 
            beta_death2018_s * death2018 + beta_0),
        r9pcesd_MNAR = expit(#beta_age_s * r9age_y_int + 
          beta_cesdpre_s * r8cesd + 
            beta_cesdcurrent_s * r9cesd +
            beta_condepre_s * r8conde_impute + 
            #beta_shltpre_s * r8shlt + 
            beta_death2018_s * death2018 + beta_0)) %>%
      select(contains("MNAR", ignore.case = FALSE))
  }
  
  # Flag the score based on bernoulli distribution, prob = p_wave_MAR/MNAR
  for (j in 1:ncol(subset)){
    subset[, paste0("r", j + 3, "cesd_missing")] <- 
      rbinom(nrow(subset), size = 1, prob = subset[[j]])
  }
  
  # Calculate the missing proportion
  subset_long <- subset %>%
    select(contains("cesd_missing")) %>%
    pivot_longer(
      everything(),
      names_to = "Variables",
      values_to = "Missingness"
    )
  
  return(
    Missing_prop_results <- tibble(
      missing_prop = mean(subset_long$Missingness, na.rm = T)
    )
  )
}

# Repeat for 1000 times
replicate = 1000
set.seed(20210507)

#---- run MAR sim ----
# Missing proportion
{missing_prop_MAR <- 
  map_dfr(1:replicate, ~ scale_coef_func(data_wide, "MAR", 1),
          .id = "replication") %>% estimate_df()

missing_prop_MAR %>%
  kbl(caption = "MAR missingness")%>%
  kable_classic(full_width = F, html_font = "Arial")
}

#---- run MNAR sim ----
{missing_prop_MNAR <- 
  map_dfr(1:replicate, ~ scale_coef_func(data_wide, "MNAR", 1),
          .id = "replication") %>% estimate_df()

missing_prop_MNAR %>%
  kbl(caption = "MNAR missingness")%>%
  kable_classic(full_width = F, html_font = "Arial")
}
