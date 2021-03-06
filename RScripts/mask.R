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
  
  #---- create incomplete data ----
  if(mechanism == "MCAR"){
    #---- MCAR ----
    #it's easier to do this with my own code than the ampute function in MICE, 
    # which requires specifying all possible missing patterns you'd like it to 
    # consider
    #   Need to add the following: wave-updated marital status, wave-updated 
    #   drinking behavior, wave-updated chronic conditions,
    mask_prop <- as.numeric(sub("%","", mask_percent))/100
    total_indices <- nrow(data_wide)*6 #6 waves of data per person
    mask_index <- sample(seq(1, total_indices), 
                         size = floor(mask_prop*total_indices), 
                         replace = FALSE)
  } else{
    
    mask_prop <- as.numeric(sub("%","", mask_percent))/100
    
    #---- expected value for predictors of missingness ----
    # Age
    overall_ages <- data_wide %>% 
      dplyr::select(paste0("r", seq(4, 9, by = 1), "age_y_int")) %>% 
      pivot_longer(everything(), names_to = "orig_varname", 
                   values_to = "age_int")
    
    # CESD previous wave (wave 3-8)
    CESD_3_8_long <- data_wide %>%
      dplyr::select(paste0("r", seq(3, 8, by = 1), "cesd")) %>%
      pivot_longer(everything(), names_to = "orig_varname", 
                   values_to = "cesd")
    
    # conde at previous wave
    conde_3_8_long <- data_wide %>%
      dplyr::select(paste0("r", seq(3, 8, by = 1), "conde_impute")) %>%
      pivot_longer(everything(), names_to = "orig_varname", 
                   values_to = "conde_impute")
    
    # CESD this wave (wave 4-9)
    CESD_4_9_long <- data_wide %>%
      dplyr::select(paste0("r", seq(4, 9, by = 1), "cesd")) %>%
      pivot_longer(everything(), names_to = "orig_varname", 
                   values_to = "cesd")
    #---- **E(X)s ----
    e_age <- mean(overall_ages$age_int)
    e_CESD_3_8 <- mean(CESD_3_8_long$cesd, na.rm = T)
    e_CESD_4_9 <- mean(CESD_4_9_long$cesd)
    e_conde <- mean(conde_3_8_long$conde_impute, na.rm = T)
    
    #---- betas ----
    beta_age <- log(1.05)
    beta_cesdpre <- log(1.10)
    beta_condepre <- log(1.30)
    beta_cesdcurrent <- log(1.15)
    
    
    if (mechanism == "MAR"){
      #---- MAR ----
      beta_0 <- logit(mask_prop) - 
        (beta_age*e_age + beta_cesdpre*e_CESD_3_8 + beta_condepre*e_conde)
      
      # # Pending to improve the efficiency
      # Not working unquote functions: 
      # for (i in (4:9)){
      #   data_wide[, paste0("r", i, "pcesd")] <-
      #     with(
      #       data_wide,
      #       expit(noquote(paste0("beta_0 + beta_age * ", "r", i, "age_y_int +",
      #                            "beta_cesdpre * ", "r", i-1, "cesd",
      #                            "beta_condepre * ", "r", i, "conde_impute"))))
      # }
      # 
        subset <- data_wide %>%
          mutate(
          r4pcesd = expit(beta_0 + beta_age * r4age_y_int + 
                                beta_cesdpre * r3cesd
                              + beta_condepre * r4conde_impute),
          r5pcesd = expit(beta_0 + beta_age * r5age_y_int + 
                                beta_cesdpre * r4cesd
                              + beta_condepre * r5conde_impute),
          r6pcesd = expit(beta_0 + beta_age * r6age_y_int + 
                                beta_cesdpre * r5cesd
                              + beta_condepre * r6conde_impute),
          r7pcesd = expit(beta_0 + beta_age * r7age_y_int + 
                                beta_cesdpre * r6cesd
                              + beta_condepre * r7conde_impute),
          r8pcesd = expit(beta_0 + beta_age * r8age_y_int + 
                                beta_cesdpre * r7cesd
                              + beta_condepre * r8conde_impute),
          r9pcesd = expit(beta_0 + beta_age * r9age_y_int + 
                                beta_cesdpre * r8cesd
                              + beta_condepre * r9conde_impute)
          ) %>%
          select(contains("pcesd"))
        
        for (j in 1:ncol(subset)){
          if (j <= 6){
            subset[, paste0("r", j + 3, "cesd_missing")] <- 
              rbinom(nrow(subset), size = 1, prob = subset[[j]])
          } 
          
      }
      
     subset_long <- subset %>%
       select(contains("cesd_missing")) %>%
        pivot_longer(everything(),
                     names_to = "orig_varname",
                     values_to = "cesd_missing")
        
      mask_index <- which(subset_long$cesd_missing == 1)
      
    } else if(mechanism == "MNAR"){
      # Pending!
      #---- MNAR ----
      beta_0 <- logit(mask_prop) - 
        (beta_age*e_age + beta_cesdpre*e_CESD_3_8 + beta_condepre*e_conde + 
          beta_cesdcurrent*e_CESD_4_9)
      
      
    }
  }
  
  
  #masking wave-specific values
  mask_wave_specific <- c("mstat_impute", "mstat_cat", "drinking_impute", 
                          "drinking_cat", "memrye_impute", "stroke_impute", 
                          "hearte_impute", "lunge_impute", "cancre_impute", 
                          "hibpe_impute", "diabe_impute", "cesd", "BMI", 
                          "shlt")
  
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
  
  
  
  # Pending
  # return()
}