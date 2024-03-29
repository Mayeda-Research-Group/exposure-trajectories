# Aim: to optimize the value of beta_0 for fixed values of other coefficients
#   in our masking models so that we get the desired level of missingness from
#   our masking function
# Created by: Yingyan Wu
# Date created: 05.07.2021
# Edited by: Crystal Shaw (05.21.2021)
# 
# Data input: CESD_data_wide.csv
# Data output: beta_0_table.csv, beta_mat.csv

#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "kableExtra")

#---- custom functions ----
expit <- function(x) {
  output <- (exp(x)/(1 + exp(x)))
  return(output)
}

logit <- function(x){
  output <- log(x/(1 - x))
  return(output)
}

estimate_df <- function(dataframe){
  estimate <- function(x){
    stats <- c(length(x), mean(x, na.rm = TRUE), 
               quantile(x, c(0.025, 0.975), na.rm = TRUE))
    names(stats) <- c("N", "Mean", "p2.5th", "p97.5th")
    return(round(stats, 4))
  }
  
  stats_tib <- dataframe %>% select(where(is.numeric)) %>% 
    purrr::map_dfr(estimate) %>% as_tibble() %>%
    mutate("Effect" = colnames(dataframe)[-1]) %>% select(Effect, everything())
  
  return(stats_tib)
}

#---- note ----
# Since the difference between win and OS, put substituted directory here
# Yingyan's directory: C:/Users/yingyan_wu/Dropbox
#                      
# Crystal's directory: ~/Dropbox/Projects 

#Changing directories here will change them throughout the script
path_to_dropbox <- "~/Dropbox/Projects"

#---- read in analytical sample ----
data_wide <- read_csv(paste0(path_to_dropbox, "/exposure_trajectories/data/", 
                             "CESD_data_wide.csv"))

for(wave in seq(4, 9)){
  data_wide %<>% 
    mutate(!!paste0("r", wave, "cesd_death2018") := 
             !!sym(paste0("r", wave, "cesd"))*death2018, 
           !!paste0("r", wave - 1, "cesd_conde_impute") := 
             !!sym(paste0("r", wave - 1, "cesd"))*
             !!sym(paste0("r", wave - 1, "conde_impute")))
}

#---- beta and expected values matrix ----
beta_mat <- #effect sizes: fixed at what we thought were reasonable values
  matrix(c(log(1.1), log(1.05), log(1.05), log(1.25), log(1.1), log(1.25),
           #expected values: used to calculate balancing intercept for warm start
           #  for optimizer
           #MAR variables
           mean(unlist(data_wide[, paste0("r", seq(3, 8, by = 1), "cesd")]), 
                na.rm = TRUE), 
           mean(unlist(data_wide[, paste0("r", seq(3, 8, by = 1), 
                                          "conde_impute")]), na.rm = TRUE), 
           mean(unlist(data_wide[, paste0("r", seq(3, 8, by = 1), 
                                          "cesd_conde_impute")]), na.rm = TRUE),
           #MNAR variables
           mean(data_wide$death2018), 
           mean(unlist(data_wide[, paste0("r", seq(4, 9, by = 1), "cesd")])), 
           mean(unlist(data_wide[, paste0("r", seq(4, 9, by = 1), 
                                          "cesd_death2018")]))), 
         nrow = 2, byrow = TRUE) %>% 
  set_colnames(c("cesdpre", "condepre", "cesdpre_condepre",
                 "death2018", "cesdcurrent", "death2018_cesdcurrent")) %>% 
  set_rownames(c("beta", "expected"))

#----  prop missing function: to be optimized ----
missing_prop <- function(BETA_0, dataset, mechanism, mask_prop, beta_mat, 
                         optimize = "yes"){
  subset <- dataset
  
  if (mechanism == "MNAR"){
    #---- MNAR ----
    for(wave in seq(4, 9)){
      subset %<>% 
        mutate(!!paste0("r", wave, "pcesd") := 
                 expit(beta_mat["beta", "death2018"]*death2018 + 
                         beta_mat["beta", "cesdcurrent"]*
                         !!sym(paste0("r", wave, "cesd")) + 
                         beta_mat["beta", "death2018_cesdcurrent"]*
                         !!sym(paste0("r", wave, "cesd_death2018")) + BETA_0))
    }
  } else if(mechanism == "MAR"){
    #---- MAR ----
    for(wave in seq(4, 9)){
      subset %<>% 
        mutate(!!paste0("r", wave, "pcesd") := 
                 expit(beta_mat["beta", "cesdpre"]*
                         !!sym(paste0("r", wave - 1, "cesd")) + 
                         beta_mat["beta", "condepre"]*
                         !!sym(paste0("r", wave - 1, "conde_impute")) + 
                         beta_mat["beta", "cesdpre_condepre"]*
                         !!sym(paste0("r", wave - 1, "cesd_conde_impute")) + 
                         BETA_0))
    }
  } 
  
  subset %<>% dplyr::select(contains("pcesd", ignore.case = FALSE))
  subset[is.na(subset)] <- 0
  
  # Flag the score based on bernoulli distribution, prob = p_wave_MAR/MNAR
  for (j in 1:ncol(subset)){
    subset[, paste0("r", j + 3, "cesd_missing")] <- 
      rbinom(nrow(subset), size = 1, prob = subset[[j]])
  }
  
  if(optimize == "yes"){
    return(abs(
      mean(subset %>% dplyr::select(contains("cesd_missing")) %>% as.matrix()) - 
        mask_prop)) 
  } else{
    #This is to double check the results missing probabilities from the optimized
    # parameters
    return(Missing_prop_results <- 
             tibble(missing_prop = 
                      mean(subset %>% 
                             dplyr::select(contains("cesd_missing")) %>% 
                             as.matrix(), na.rm = T)))
  }
}

#---- optimizer ----
#Balancing intercepts (Rudolph et.al., 2021) are used for warm starts 
warm_start <- function(mask_prop, mechanism, beta_mat){
  if(mechanism == "MNAR"){
    return(-log(1/(mask_prop) - 1) -
             (beta_mat["beta", "death2018"]*beta_mat["expected", "death2018"] + 
                beta_mat["beta", "cesdcurrent"]*
                beta_mat["expected", "cesdcurrent"] +
                beta_mat["beta", "death2018_cesdcurrent"]*
                beta_mat["expected", "death2018_cesdcurrent"]))
  } else{
    return(-log(1/(mask_prop) - 1) -
             (beta_mat["beta", "cesdpre"]*beta_mat["expected", "cesdpre"] + 
                beta_mat["beta", "condepre"]*beta_mat["expected", "condepre"] +
                beta_mat["beta", "cesdpre_condepre"]*
                beta_mat["expected", "cesdpre_condepre"]))
  }
} 

for(mechanism in c("MAR", "MNAR")){
  for(mask_prop in c(0.10, 0.20, 0.30)){
    start <- warm_start(mask_prop, mechanism, beta_mat)
    
    # Assign dataset name to the datasets generated by `optimize` function
    assign(paste0("optim_", mechanism, 100*mask_prop), 
           optimize(missing_prop, lower = start + 3*start, upper = 0, 
                    maximum = FALSE, dataset = data_wide, mechanism = mechanism, 
                    mask_prop = mask_prop, beta_mat = beta_mat, 
                    optimize = "yes"))
  }
}

#---- run sims ----
#Check the optimized values using simulation
replicate = 1000
set.seed(20210507)
# Missing proportion
{test_run <- 
    map_dfr(1:replicate, 
            ~ missing_prop(BETA_0 = optim_MNAR30$minimum, dataset = data_wide, 
                           mechanism = "MNAR", mask_prop = 0.30, 
                           beta_mat = beta_mat, optimize = "No"),
            .id = "replication") %>% estimate_df()
  
  test_run %>%
    kbl(caption = "MAR missingness")%>%
    kable_classic(full_width = F, html_font = "Arial")
}

#---- save results ----
#---- **optimized betas table ----
#don't need MCAR
mechanisms <- c("MAR", "MNAR")
percents <- c(10, 20, 30)

beta_0_table <- expand_grid(mechanisms, percents) %>% 
  mutate("beta0" = 0)

for(mechanism in mechanisms){
  for(percent in percents){
    optimized <- get(paste0("optim_", mechanism, percent))
    
    beta_0_table[which(beta_0_table$mechanisms == mechanism & 
                         beta_0_table$percents == percent), "beta0"] <- 
      optimized$minimum
  }
}

write_csv(beta_0_table, paste0(path_to_dropbox, "/exposure_trajectories/data/", 
                               "beta_0_table.csv"))

#---- **beta matrix ----
write_csv(as.data.frame(t(beta_mat["beta", ])), 
          paste0(path_to_dropbox, "/exposure_trajectories/data/beta_mat.csv"))

