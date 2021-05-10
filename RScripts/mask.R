#---- Expit function ----
expit <- function(x) {
  output <- (exp(x)/(1+exp(x)))
  return(output)
}

logit <- function(x){
  output <- log(x/(1-x))
  return(output)
}

#---- Mask function ----
mask <- function(data_wide, mechanism, mask_percent){
  
  mask_prop <- as.numeric(sub("%","", mask_percent))/100
  
  #---- create incomplete data ----
  if(mechanism == "MCAR"){
    #---- MCAR ----
    #it's easier to do this with my own code than the ampute function in MICE, 
    # which requires specifying all possible missing patterns you'd like it to 
    # consider
    #   Need to add the following: wave-updated marital status, wave-updated 
    #   drinking behavior, wave-updated chronic conditions,
    total_indices <- nrow(data_wide)*6 #6 waves of data per person
    mask_index <- sample(seq(1, total_indices), 
                         size = floor(mask_prop*total_indices), 
                         replace = FALSE)
  } else{
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
    
    if (mechanism == "MAR"){
      #---- ** MAR Coefficients ----
      beta_age_10 <- log(0.96)
      beta_cesdpre_10 <- log(1.04)
      beta_condepre_10 <- log(1.30)
      beta_death2018 <- log(2.25)
      
      beta_0_10 <- logit(0.1) -
        (beta_age_10*e_age + beta_cesdpre_10*e_CESD_3_8 + 
           beta_condepre_10*e_conde + beta_death2018_10*e_death2018)
      
      if (mask_prop == 0.1){
        beta_0 <- beta_0_10
        beta_age <- beta_age_10
        beta_cesdpre <- beta_cesdpre_10
        beta_condepre <- beta_condepre_10
        beta_death2018 <- beta_death2018_10
      } else{
        scaling_coef <- logit(mask_prop)/logit(0.1)
        beta_0 <- scaling_coef * beta_0_10
        beta_age <- scaling_coef * beta_age_10
        beta_cesdpre <- scaling_coef * beta_cesdpre_10
        beta_condepre <- scaling_coef * beta_condepre_10
        beta_shlt <- scaling_coef * beta_shltpre_10
        beta_death2018 <- scaling_coef * beta_death2018_10
      }
      
      subset <- data_wide %>%
        mutate(
          r4pcesd = ifelse(is.na(expit(beta_0 + beta_age * r4age_y_int + 
                                         beta_cesdpre * r3cesd + 
                                         beta_condepre * r3conde_impute + 
                                         beta_death2018 * death2018)),
                           0, expit(beta_0 + beta_age * r4age_y_int +
                                      beta_cesdpre * r3cesd + 
                                      beta_condepre * r3conde_impute + 
                                      beta_death2018 * death2018)),
          r5pcesd = expit(beta_0 + beta_age * r5age_y_int + 
                            beta_cesdpre * r4cesd + 
                            beta_condepre * r4conde_impute + 
                            beta_death2018 * death2018),
          r6pcesd = expit(beta_0 + beta_age * r6age_y_int + 
                            beta_cesdpre * r5cesd +
                            beta_condepre * r5conde_impute + 
                            beta_death2018 * death2018),
          r7pcesd = expit(beta_0 + beta_age * r7age_y_int + 
                            beta_cesdpre * r6cesd + 
                            beta_condepre * r6conde_impute + 
                            beta_death2018 * death2018),
          r8pcesd = expit(beta_0 + beta_age * r8age_y_int + 
                            beta_cesdpre * r7cesd +
                            beta_condepre * r7conde_impute + 
                            beta_death2018 * death2018),
          r9pcesd = expit(beta_0 + beta_age * r9age_y_int + 
                            beta_cesdpre * r8cesd +
                            beta_condepre * r8conde_impute + 
                            beta_death2018 * death2018)) %>%
        select(contains("pcesd"))
      
      for (j in 1:ncol(subset)){
        subset[, paste0("r", j + 3, "cesd_missing")] <- 
          rbinom(nrow(subset), size = 1, prob = subset[[j]])
      }
      
      subset_long <- subset %>%
        select(contains("cesd_missing")) %>%
        pivot_longer(everything(),
                     names_to = "orig_varname",
                     values_to = "cesd_missing")
      
      mask_index <- which(subset_long$cesd_missing == 1)
      
    } else if(mechanism == "MNAR"){
      #---- ** MNAR Coefficients ----
      beta_age_10 <- log(0.955)
      beta_cesdpre_10 <- log(1.04)
      beta_condepre_10 <- log(1.20)
      beta_death2018_10 <- log(2.25)
      beta_cesdcurrent_10 <- log(1.25)
      
      beta_0_10 <- logit(0.1) -
        (beta_age_10*e_age + beta_cesdpre_10*e_CESD_3_8 + 
           beta_condepre_10*e_conde + beta_cesdcurrent_10*e_CESD_4_9 + 
           beta_death2018_10*e_death2018)
      
      if (mask_prop == 0.1){
        beta_0 <- beta_0_10
        beta_age <- beta_age_10
        beta_cesdpre <- beta_cesdpre_10
        beta_condepre <- beta_condepre_10
        beta_death2018 <- beta_death2018_10
        beta_cesdcurrent <- beta_cesdcurrent_10
      } else {
        scaling_coef <- logit(mask_prop)/logit(0.1)
        beta_0 <- scaling_coef * beta_0_10
        beta_age <- scaling_coef * beta_age_10
        beta_cesdpre <- scaling_coef * beta_cesdpre_10
        beta_condepre <- scaling_coef * beta_condepre_10
        beta_shlt <- scaling_coef * beta_shltpre_10
        beta_death2018 <- scaling_coef * beta_death2018_10
        beta_cesdcurrent <- scaling_coef * beta_cesdcurrent_10
      }
      
      subset <- data_wide %>%
        mutate(
          r4pcesd = ifelse(is.na(expit(beta_0 + beta_age * r4age_y_int + 
                                         beta_cesdpre * r3cesd + 
                                         beta_cesdcurrent * r4cesd +
                                         beta_condepre * r3conde_impute + 
                                         beta_death2018 * death2018)), 
                           0, expit(beta_0 + beta_age * r4age_y_int + 
                                      beta_cesdpre * r3cesd + 
                                      beta_cesdcurrent * r4cesd +
                                      beta_condepre * r3conde_impute + 
                                      beta_death2018 * death2018)),
          r5pcesd = expit(beta_0 + beta_age * r5age_y_int + 
                            beta_cesdpre * r4cesd + beta_cesdcurrent * r5cesd +
                            beta_condepre * r4conde_impute + 
                            beta_death2018 * death2018),
          r6pcesd = expit(beta_0 + beta_age * r6age_y_int + 
                            beta_cesdpre * r5cesd + beta_cesdcurrent * r6cesd +
                            beta_condepre * r5conde_impute + 
                            beta_death2018 * death2018),
          r7pcesd = expit(beta_0 + beta_age * r7age_y_int + 
                            beta_cesdpre * r6cesd + beta_cesdcurrent * r7cesd +
                            beta_condepre * r6conde_impute + 
                            beta_death2018 * death2018),
          r8pcesd = expit(beta_0 + beta_age * r8age_y_int + 
                            beta_cesdpre * r7cesd + beta_cesdcurrent * r8cesd +
                            beta_condepre * r7conde_impute + 
                            beta_death2018 * death2018),
          r9pcesd = expit(beta_0 + beta_age * r9age_y_int + 
                            beta_cesdpre * r8cesd + beta_cesdcurrent * r9cesd +
                            beta_condepre * r8conde_impute + 
                            beta_death2018 * death2018)) %>%
        select(contains("pcesd"))
      
      for (j in 1:ncol(subset)){
        subset[, paste0("r", j + 3, "cesd_missing")] <- 
          rbinom(nrow(subset), size = 1, prob = subset[[j]])
      }
      
      subset_long <- subset %>%
        select(contains("cesd_missing")) %>%
        pivot_longer(everything(),
                     names_to = "orig_varname",
                     values_to = "cesd_missing")
      
      mask_index <- which(subset_long$cesd_missing == 1)
    }
  }
  
  #---- masking wave-specific values ----
  mask_wave_specific <- c("married_partnered", "not_married_partnered", 
                          "widowed", "drinking_cat", "memrye_impute", 
                          "stroke_impute", "hearte_impute", "lunge_impute", 
                          "cancre_impute", "hibpe_impute", "diabe_impute", 
                          "cesd", "BMI", "shlt")
  
  for(var in mask_wave_specific){
    #mask values
    data_long <- data_wide %>% 
      dplyr::select("HHIDPN", paste0("r", seq(4, 9), var)) %>% 
      pivot_longer(-"HHIDPN")
    
    data_long[mask_index, "value"] <- NA
    
    #need to get rid of rows with NA values to pivot back to wide
    data_long %<>% na.omit() %>%
      pivot_wider(names_from = "name", values_from = "value")
    
    #mask in wide dataset with masked values
    data_wide[which(data_wide$HHIDPN %in% data_long$HHIDPN),
              paste0("r", seq(4, 9), var)] <- data_long[, -1]
    data_wide[which(!data_wide$HHIDPN %in% data_long$HHIDPN),
              paste0("r", seq(4, 9), var)] <- NA
  }
  
  #masking derived values
  data_wide %<>% 
    mutate("r4cesd_elevated" = ifelse(r4cesd > 4, 1, 0), 
           "r9cesd_elevated" = ifelse(r9cesd > 4, 1, 0), 
           "total_elevated_cesd" = 
             rowSums(data_wide %>% 
                       dplyr::select(paste0("r", seq(4, 9), "cesd")) %>% 
                       mutate_all(function(x) x > 4)), 
           "avg_cesd" = 
             rowMeans(data_wide %>% 
                        dplyr::select(paste0("r", seq(4, 9), "cesd"))), 
           "avg_cesd_elevated" = ifelse(avg_cesd > 4, 1, 0))
  
  #---- check missings ----
  # make sure no one is missing every cesd measure
  data_wide %<>%
    mutate("CESD_missing" =
             rowSums(data_wide %>%
                       dplyr::select(paste0("r", seq(4, 9), "cesd")) %>%
                       mutate_all(function(x) is.na(x)))) %>%
    filter(CESD_missing < 6) %>%
    select(-CESD_missing)
  
  # Return the dataset
  # return(mask_index)
  return(data_wide)
}

# # Test
# path_to_box <- "C:/Users/yingyan_wu"
# path_to_dropbox <- "C:/Users/yingyan_wu/Dropbox"
# 
# #---- read in analytical sample ----
# CESD_data_wide <-
#   read_csv(paste0(path_to_dropbox,
#                   "/exposure_trajectories/data/",
#                   "CESD_data_wide.csv"))
# # geting the coefficients in MAR and MNAR models
# coef_df <- data.frame(
#   intercept = exp(beta_0),
#   age_j = exp(beta_age),
#   CESD_j_1 = exp(beta_cesdpre),
#   conde_j_1 = exp(beta_condepre),
#   cesd_j = exp(beta_cesdcurrent)
# ) %>% print()

# MCAR
# test <- mask(CESD_data_wide, "MCAR", "20%")
# {index <- mask(CESD_data_wide, "MCAR", "20%")
#   length(index)/(9445*6)
# }
# test %>%
#   select(contains("cesd"),
#          -contains(c("elevated", "avg", "r3"))) %>%
#   pivot_longer(everything(),
#                names_to = "origvar",
#                values_to = "CESD") %>%
#     dplyr::count(CESD) %>%
#     mutate(prop = prop.table(n))
# 
# view(test[is.na(test$r9cesd), ])
# 
# test_func <- function(data_wide, mechanism, mask_percent){
#   mask_prop <- as.numeric(sub("%","", mask_percent))/100
#   n <- nrow(mask(data_wide, mechanism, mask_percent))
#   return(n)
# }
# 
# size = 1000
# system.time(
#   test_sample_size <- replicate(size, test_func(CESD_data_wide, "MCAR", "30%")))
# summary(test_sample_size)
# 
# MAR
# test1 <- mask(CESD_data_wide, "MAR", "20%")
# test1 %>%
#   select(contains("cesd"),
#          -contains(c("elevated", "avg", "r3"))) %>%
#   pivot_longer(everything(),
#                names_to = "origvar",
#                values_to = "CESD") %>%
#   dplyr::count(CESD) %>%
#   mutate(prop = prop.table(n))
# 
# test1 %>%
#   select(contains("r5")) %>%
#   filter(is.na(r5cesd))
# 
# size = 1000
# system.time(
#   test1_sample_size <- replicate(size, test_func(CESD_data_wide, "MAR", "30%")))
# summary(test1_sample_size)

# # MNAR
# test2 <- mask(CESD_data_wide, "MNAR", "30%")
# test2 %>%
#   select(contains("cesd"),
#          -contains(c("elevated", "avg", "r3"))) %>%
#   pivot_longer(everything(),
#                names_to = "origvar",
#                values_to = "CESD") %>%
#   dplyr::count(CESD) %>%
#   mutate(prop = prop.table(n))
# test2 %>%
#   select(contains("r5")) %>%
#   filter(is.na(r5cesd))
# 
# size = 1000
# system.time(
#   test2_sample_size <- replicate(size, test_func(CESD_data_wide, "MNAR", "30%")))
# summary(test2_sample_size)
