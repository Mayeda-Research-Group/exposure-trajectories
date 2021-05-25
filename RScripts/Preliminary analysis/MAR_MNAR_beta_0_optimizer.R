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
    stats <- c(length(x), mean(x,na.rm = T), 
               quantile(x,c(0.025, 0.975),na.rm = T))
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
# Yingyan's directory: C:/Users/yingyan_wu
#                      C:/Users/yingyan_wu/Dropbox
# Crystal's directory: /Users/CrystalShaw
#                     ~/Dropbox/Projects

#Changing directories here will change them throughout the script
path_to_box <- "/Users/CrystalShaw"
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
beta_mat <- #effect sizes
  matrix(c(log(1.1), log(1.15), log(1.15), log(1.25), log(1.1), log(1.25),
           #expected values
           mean(unlist(data_wide[, paste0("r", seq(3, 8, by = 1), "cesd")]), 
                na.rm = TRUE), 
           mean(unlist(data_wide[, paste0("r", seq(3, 8, by = 1), 
                                          "conde_impute")]), na.rm = TRUE), 
           mean(unlist(data_wide[, paste0("r", seq(3, 8, by = 1), 
                                          "cesd_conde_impute")]), na.rm = TRUE), 
           mean(data_wide$death2018), 
           mean(unlist(data_wide[, paste0("r", seq(4, 9, by = 1), "cesd")])), 
           mean(unlist(data_wide[, paste0("r", seq(4, 9, by = 1), 
                                          "cesd_death2018")]))), 
         nrow = 2, byrow = TRUE) %>% 
  #MAR
  set_colnames(c("cesdpre", "condepre", "cesdpre_condepre",
                 #MNAR
                 "death2018", "cesdcurrent", "death2018_cesdcurrent")) %>% 
  set_rownames(c("beta", "expected"))

#----  prop missing function ----
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
    return(Missing_prop_results <- 
             tibble(missing_prop = 
                      mean(subset %>% 
                             dplyr::select(contains("cesd_missing")) %>% 
                             as.matrix(), na.rm = T)))
  }
}

#---- optimizer ----
#---- **MNAR ----
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
    
    assign(paste0("optim_", mechanism, 100*mask_prop), 
           optimize(missing_prop, lower = start + 3*start, upper = 0, 
                    maximum = FALSE, dataset = data_wide, mechanism = mechanism, 
                    mask_prop = mask_prop, beta_mat = beta_mat, 
                    optimize = "yes"))
  }
}

# #---- test masking ----
# test_mask <- function (dataset, mechanism, beta_0){
#   subset <- dataset
#   
#   if (mechanism == "MNAR"){
#     #---- MNAR ----
#     for(wave in seq(4, 9)){
#       subset %<>% 
#         mutate(!!paste0("r", wave, "pcesd") := 
#                  expit(beta_death2018*death2018 + 
#                          beta_cesdcurrent*
#                          !!sym(paste0("r", wave, "cesd")) + 
#                          beta_death2018_cesdcurrent*
#                          death2018*!!sym(paste0("r", wave, "cesd")) + beta_0))
#     }
#   } else if(mechanism == "MAR"){
#     #---- MAR ----
#     for(wave in seq(4, 9)){
#       subset %<>% 
#         mutate(!!paste0("r", wave, "pcesd") := 
#                  expit(beta_cesdpre*
#                          !!sym(paste0("r", wave - 1, "cesd")) + 
#                          beta_condepre*
#                          !!sym(paste0("r", wave - 1, "conde_impute")) + 
#                          beta_cesdpre_condepre*
#                          !!sym(paste0("r", wave - 1, "cesd"))*
#                          !!sym(paste0("r", wave - 1, "conde_impute")) + beta_0))
#     }
#   }
#   
#   subset %<>% dplyr::select(contains("pcesd", ignore.case = FALSE))
#   
#   # Flag the score based on bernoulli distribution, prob = p_wave_MAR/MNAR
#   for (j in 1:ncol(subset)){
#     subset[, paste0("r", j + 3, "cesd_missing")] <- 
#       rbinom(nrow(subset), size = 1, prob = subset[[j]])
#   }
#   
#   return(Missing_prop_results <- 
#            tibble(missing_prop = 
#                     mean(subset %>% dplyr::select(contains("cesd_missing")) %>% 
#                            as.matrix(), na.rm = T)))
# }

# #---- run MNAR sim ----
# # Repeat for 1000 times
# replicate = 1000
# set.seed(20210507)
# 
# {missing_prop_MNAR <- 
#     map_dfr(1:replicate, ~ test_mask(dataset = data_wide, mechanism = "MNAR", 
#                                      beta_0 = optim_MNAR30$minimum), 
#             .id = "replication") %>% estimate_df()
#   
#   missing_prop_MNAR %>%
#     kbl(caption = "MNAR missingness")%>%
#     kable_classic(full_width = F, html_font = "Arial")
# }
# 
# #---- **save optimized beta_0 ----
# for(percent in c(10, 20, 30)){
#   write_rds(get(paste0("optim_MNAR", percent)), 
#             file = paste0(path_to_dropbox, "/exposure_trajectories/data/", 
#                           "optimized_masking_intercepts/optim_MNAR", percent, 
#                           ".RDS"))
# }

#---- run MAR sim ----
# Missing proportion
{missing_prop_MAR <-
  map_dfr(1:replicate, ~ test_mask(dataset = data_wide, mechanism = "MAR", 
                                   beta_0 = optim_MAR30$minimum),
          .id = "replication") %>% estimate_df()

missing_prop_MAR %>%
  kbl(caption = "MAR missingness")%>%
  kable_classic(full_width = F, html_font = "Arial")
}

#---- **save optimized beta_0 ----
for(percent in c(10, 20, 30)){
  write_rds(get(paste0("optim_MAR", percent)), 
            file = paste0(path_to_dropbox, "/exposure_trajectories/data/", 
                          "optimized_masking_intercepts/optim_MAR", percent, 
                          ".RDS"))
}

