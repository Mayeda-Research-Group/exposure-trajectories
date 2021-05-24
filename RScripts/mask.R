#---- read in optimized betas ----
for(mechanism in c("MNAR")){
  for(percent in c(10, 20, 30)){
    assign(paste0("beta0_", mechanism, percent), 
           read_rds(file = 
                      paste0(path_to_dropbox, "/exposure_trajectories/data/", 
                             "optimized_masking_intercepts/optim_", mechanism, 
                             percent, ".RDS")))
  }
}

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
      MAR_coeff <- tibble("Var" = c("cesdpre", "condepre", "death2018"), 
                          "10%" = c(log(2.028), log(2.535), log(5.85)), 
                          "20%" = c(log(2.86), log(3.575), log(8.25)), 
                          "30%" = c(log(4.16), log(5.2), log(12)))
      
      #beta_age <- MAR_coeff[[which(MAR_coeff$Var == "age"), mask_percent]]
      beta_cesdpre <- 
        MAR_coeff[[which(MAR_coeff$Var == "cesdpre"), mask_percent]]
      beta_condepre <- 
        MAR_coeff[[which(MAR_coeff$Var == "condepre"), mask_percent]]
      beta_death2018 <- 
        MAR_coeff[[which(MAR_coeff$Var == "death2018"), mask_percent]]
      
      beta_0 <- logit(mask_prop) -
        (#beta_age*e_age + 
          beta_cesdpre*e_CESD_3_8 + beta_condepre*e_conde +
            beta_death2018*e_death2018)
      
      subset <- data_wide %>%
        mutate(
          r4pcesd = ifelse(is.na(expit(beta_0 + #beta_age * r4age_y_int + 
                                         beta_cesdpre * r3cesd + 
                                         beta_condepre * r3conde_impute + 
                                         beta_death2018 * death2018)),
                           0, expit(beta_0 + #beta_age * r4age_y_int +
                                      beta_cesdpre * r3cesd + 
                                      beta_condepre * r3conde_impute + 
                                      beta_death2018 * death2018)),
          r5pcesd = expit(beta_0 + #beta_age * r5age_y_int + 
                            beta_cesdpre * r4cesd + 
                            beta_condepre * r4conde_impute + 
                            beta_death2018 * death2018),
          r6pcesd = expit(beta_0 + #beta_age * r6age_y_int + 
                            beta_cesdpre * r5cesd +
                            beta_condepre * r5conde_impute + 
                            beta_death2018 * death2018),
          r7pcesd = expit(beta_0 + #beta_age * r7age_y_int + 
                            beta_cesdpre * r6cesd + 
                            beta_condepre * r6conde_impute + 
                            beta_death2018 * death2018),
          r8pcesd = expit(beta_0 + #beta_age * r8age_y_int + 
                            beta_cesdpre * r7cesd +
                            beta_condepre * r7conde_impute + 
                            beta_death2018 * death2018),
          r9pcesd = expit(beta_0 + #beta_age * r9age_y_int + 
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
      MNAR_coeff <- 
        tibble("Var" = c("age", "cesdpre", "condepre", "death2018", "cesdcurrent"), 
               "10%" = c(log(0.9562), log(1.0413), log(1.1514), log(2.2528), 
                         log(1.3017)), 
               "20%" = c(log(0.967), log(1.053), log(1.164), log(2.278), 
                         log(1.316)), 
               "30%" = c(log(0.9743), log(1.0611), log(1.1733), log(2.2956), 
                         log(1.3263)))
      
      beta_age <- MNAR_coeff[[which(MNAR_coeff$Var == "age"), mask_percent]]
      beta_cesdpre <- 
        MNAR_coeff[[which(MNAR_coeff$Var == "cesdpre"), mask_percent]]
      beta_condepre <- 
        MNAR_coeff[[which(MNAR_coeff$Var == "condepre"), mask_percent]]
      beta_death2018 <- 
        MNAR_coeff[[which(MNAR_coeff$Var == "death2018"), mask_percent]]
      beta_cesdcurrent <- 
        MNAR_coeff[[which(MNAR_coeff$Var == "cesdcurrent"), mask_percent]]
      
      beta_0 <- logit(mask_prop) -
        (beta_age*e_age + beta_cesdpre*e_CESD_3_8 + beta_condepre*e_conde + 
           beta_cesdcurrent*e_CESD_4_9 + beta_death2018*e_death2018)
      
      subset <- data_wide %>%
        mutate(
          r4pcesd = ifelse(is.na(expit(beta_0 + #beta_age * r4age_y_int + 
                                         beta_cesdpre * r3cesd + 
                                         beta_cesdcurrent * r4cesd +
                                         beta_condepre * r3conde_impute + 
                                         beta_death2018 * death2018)), 
                           0, expit(beta_0 + #beta_age * r4age_y_int + 
                                      beta_cesdpre * r3cesd + 
                                      beta_cesdcurrent * r4cesd +
                                      beta_condepre * r3conde_impute + 
                                      beta_death2018 * death2018)),
          r5pcesd = expit(beta_0 + #beta_age * r5age_y_int + 
                            beta_cesdpre * r4cesd + beta_cesdcurrent * r5cesd +
                            beta_condepre * r4conde_impute + 
                            beta_death2018 * death2018),
          r6pcesd = expit(beta_0 + #beta_age * r6age_y_int + 
                            beta_cesdpre * r5cesd + beta_cesdcurrent * r6cesd +
                            beta_condepre * r5conde_impute + 
                            beta_death2018 * death2018),
          r7pcesd = expit(beta_0 + #beta_age * r7age_y_int + 
                            beta_cesdpre * r6cesd + beta_cesdcurrent * r7cesd +
                            beta_condepre * r6conde_impute + 
                            beta_death2018 * death2018),
          r8pcesd = expit(beta_0 + #beta_age * r8age_y_int + 
                            beta_cesdpre * r7cesd + beta_cesdcurrent * r8cesd +
                            beta_condepre * r7conde_impute + 
                            beta_death2018 * death2018),
          r9pcesd = expit(beta_0 + #beta_age * r9age_y_int + 
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
