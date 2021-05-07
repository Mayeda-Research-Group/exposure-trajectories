# Aim: to find out  the coefficients which will produce approximately the proportion
# we want for the missingness, and thus reduce the intercepts' effect
# Created by: Yingyan Wu
# 05.07.2021

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
path_to_box <- "C:/Users/yingyan_wu"
path_to_dropbox <- "C:/Users/yingyan_wu/Dropbox"

#---- read in analytical sample ----
CESD_data_wide <- 
  read_csv(paste0(path_to_dropbox, 
                  "/exposure_trajectories/data/", 
                  "CESD_data_wide.csv"))

#---- scaling the coefficients ----

#---- ** Coefficients ----
beta_age <- log(1.040)
beta_cesdpre <- log(1.050)
beta_condepre <- log(1.250)
beta_shltpre <- log(1.200)
beta_death2018 <- log(2.50)
beta_cesdcurrent <- log(1.15)

#----  scale coef function ----
scale_coef_func <- function (dataset, mechanism, scale){
  
  n <- nrow(dataset)
  beta_age_s <- scale * beta_age
  beta_cesdpre_s <- scale * beta_cesdpre
  beta_condepre_s <- scale * beta_condepre
  beta_shltpre_s <- scale * beta_shltpre
  beta_death2018_s <- scale * beta_death2018
  beta_cesdcurrent_s <- scale * beta_cesdcurrent
  
  if (mechanism == "MAR"){
    #---- MAR ----
    subset <- dataset %>%
      mutate( 
        r4pcesd_MAR = expit(beta_age_s * r4age_y_int + 
                              beta_cesdpre_s * r3cesd + 
                              beta_condepre_s * r3conde_impute + 
                              beta_shltpre_s * r3shlt + 
                              beta_death2018_s * death2018),
        r5pcesd_MAR = expit(beta_age_s * r5age_y_int + 
                              beta_cesdpre_s * r4cesd + 
                              beta_condepre_s * r4conde_impute + 
                              beta_shltpre_s * r4shlt + 
                              beta_death2018_s * death2018),
        r6pcesd_MAR = expit(beta_age_s * r6age_y_int + 
                              beta_cesdpre_s * r5cesd + 
                              beta_condepre_s * r5conde_impute + 
                              beta_shltpre_s * r5shlt + 
                              beta_death2018_s * death2018),
        r7pcesd_MAR = expit(beta_age_s * r7age_y_int + 
                              beta_cesdpre_s * r6cesd + 
                              beta_condepre_s * r6conde_impute + 
                              beta_shltpre_s * r6shlt + 
                              beta_death2018_s * death2018),
        r8pcesd_MAR = expit(beta_age_s * r8age_y_int + 
                              beta_cesdpre_s * r7cesd + 
                              beta_condepre_s * r7conde_impute + 
                              beta_shltpre_s * r7shlt + 
                              beta_death2018_s * death2018),
        r9pcesd_MAR = expit(beta_age_s * r9age_y_int + 
                              beta_cesdpre_s * r8cesd + 
                              beta_condepre_s * r8conde_impute + 
                              beta_shltpre_s * r8shlt + 
                              beta_death2018_s * death2018)
      ) %>%
      select(contains("MAR"))
  } else if(mechanism == "MNAR"){
    #---- MNAR ----
    subset <- dataset %>%
      mutate( 
        r4pcesd_MNAR = expit(beta_age_s * r4age_y_int + 
                               beta_cesdpre_s * r3cesd + 
                               beta_cesdcurrent_s * r4cesd +
                               beta_condepre_s * r3conde_impute + 
                               beta_shltpre_s * r3shlt + 
                               beta_death2018_s * death2018),
        r5pcesd_MNAR = expit(beta_age_s * r5age_y_int + 
                               beta_cesdpre_s * r4cesd + 
                               beta_cesdcurrent_s * r5cesd +
                               beta_condepre_s * r4conde_impute + 
                               beta_shltpre_s * r4shlt + 
                               beta_death2018_s * death2018),
        r6pcesd_MNAR = expit(beta_age_s * r6age_y_int + 
                               beta_cesdpre_s * r5cesd + 
                               beta_cesdcurrent_s * r6cesd +
                               beta_condepre_s * r5conde_impute + 
                               beta_shltpre_s * r5shlt + 
                               beta_death2018_s * death2018),
        r7pcesd_MNAR = expit(beta_age_s * r7age_y_int + 
                               beta_cesdpre_s * r6cesd + 
                               beta_cesdcurrent_s * r7cesd +
                               beta_condepre_s * r6conde_impute + 
                               beta_shltpre_s * r6shlt + 
                               beta_death2018_s * death2018),
        r8pcesd_MNAR = expit(beta_age_s * r8age_y_int + 
                               beta_cesdpre_s * r7cesd + 
                               beta_cesdcurrent_s * r8cesd +
                               beta_condepre_s * r7conde_impute + 
                               beta_shltpre_s * r7shlt + 
                               beta_death2018_s * death2018),
        r9pcesd_MNAR = expit(beta_age_s * r9age_y_int + 
                               beta_cesdpre_s * r8cesd + 
                               beta_cesdcurrent_s * r9cesd +
                               beta_condepre_s * r8conde_impute + 
                               beta_shltpre_s * r8shlt + 
                               beta_death2018_s * death2018)
      ) %>%
      select(contains(c("MNAR")))
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
replicate = 10
set.seed(20210507)

# Missing proportion grid
missing_prop_MAR <- 
  map_dfr(1:replicate, ~ scale_coef_func(CESD_data_wide, "MAR", 3),
          .id = "replication") %>% estimate_df()

missing_prop_MAR %>%
  kbl(caption = "MAR missingness")%>%
  kable_classic(full_width = F, html_font = "Arial")
