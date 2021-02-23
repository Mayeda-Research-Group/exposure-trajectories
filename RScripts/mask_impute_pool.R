mask_impute_pool <- 
  function(data_wide, mechanism, method, mask_percent, num_impute, save = "no"){
    #---- create shell for output ----
    exposures <- c("CES-D Wave 4", "CES-D Wave 9", "Elevated CES-D Count", 
                   "Elevated Average CES-D")
    
    pooled_effect_ests <- 
      data.frame("Exposure" = exposures, "beta" = NA, "LCI" = NA, "UCI" = NA, 
                 "Method" = method, "Missingness" = mask_percent, 
                 "Type" = mechanism)
    
    #---- create incomplete data ----
    if(mechanism == "MCAR"){
      #---- **MCAR ----
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
    }
    #making wave-specific values
    #---- EDIT HERE ----
    mask_wave_specific <- c("cesd", "BMI", "shlt")
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
    #make sure no one is missing every cesd measure
    data_wide %<>% 
      mutate("CESD_missing" = 
               rowSums(data_wide %>% 
                         dplyr::select(paste0("r", seq(4, 9), "cesd")) %>% 
                         mutate_all(function(x) is.na(x)))) %>% 
      filter(CESD_missing < 6)
    
    #Sanity check-- table of num missings per masking proportion
    # missings <- as.data.frame(matrix(nrow = length(mask_props), ncol = 8)) %>%
    #   set_colnames(c("Mask Prop", seq(0, 6)))
    # missings[, "Mask Prop"] <- 100*mask_props
    # 
    # for(mask in 100*mask_props){
    #   data <- get(paste0("mask", mask)) %>%
    #     dplyr::select(paste0("r", seq(4, 9), "cesd"))
    #   missing_count <- table(rowSums(is.na(data)))
    #   missings[which(missings$`Mask Prop` == mask), names(missing_count)] <-
    #     missing_count
    # }
    # 
    # write_csv(missings, path = paste0(path_to_dropbox,
    #                                   "/exposure_trajectories/manuscript/",
    #                                   "tables/missing_counts.csv"))
    
    if(method == "JMVN"){
      #---- transformations ----
      #Taking the log for joint MVN models
      for(wave in 4:9){
        data_wide[, paste0("logr", wave, "cesd")] <- 
          log(1 + data_wide[, paste0("r", wave, "cesd")])
      }
    }
    
    #---- imputation ----
    if(method == "JMVN"){
      #---- **JMVN ----
      #Joint multivariate normal
      #---- ***predictor matrix ----
      predict <- 
        matrix(1, nrow = 23, ncol = ncol(data_wide)) %>% 
        set_rownames(c(paste0("logr", seq(4, 9), "cesd"), 
                       apply(expand.grid("r", seq(4, 9), 
                                         mask_wave_specific[2:3]), 1, 
                             paste, collapse = ""),
                       "r4cesd_elevated", "r9cesd_elevated", 
                       "total_elevated_cesd", "avg_cesd", 
                       "avg_cesd_elevated")) %>% 
        set_colnames(colnames(data_wide))
      
      #Don't use these as predictors
      predict[, c("HHIDPN", "conde", "age_death_y", 
                  paste0("r", seq(4, 9), "cesd"), "observed", 
                  "CESD_missing")] <- 0
      
      #---- EDIT HERE ----
      #---- ****logCESD models ----
      predict[paste0("logr", seq(4, 9), "cesd"), 
              paste0("r", seq(4, 9), "BMI")] <- diag(x = 1, nrow = 6, ncol = 6)
      predict[paste0("logr", seq(4, 9), "cesd"), 
              paste0("r", seq(4, 9), "age_y_int")] <- 
        diag(x = 1, nrow = 6, ncol = 6)
      predict[paste0("logr", seq(4, 9), "cesd"), 
              paste0("r", seq(4, 9), "shlt")] <- diag(x = 1, nrow = 6, ncol = 6)
      predict[paste0("logr", seq(4, 9), "cesd"), 
              paste0("logr", seq(4, 9), "cesd")] <- 
        (diag(x = 1, nrow = 6, ncol = 6) == 0)*1
      predict[paste0("logr", seq(4, 9), "cesd"), "r4cesd_elevated"] <- 
        c(1, rep(0, 5))
      predict[paste0("logr", seq(4, 9), "cesd"), "r9cesd_elevated"] <- 
        c(rep(0, 5), 1)
      
      #---- ****BMI models ----
      predict[paste0("r", seq(4, 9), "BMI"), 
              paste0("r", seq(4, 9), "BMI")] <- 
        (diag(x = 1, nrow = 6, ncol = 6) == 0)*1
      predict[paste0("r", seq(4, 9), "BMI"), 
              paste0("r", seq(4, 9), "cesd")] <- 
        diag(x = 1, nrow = 6, ncol = 6)
      predict[paste0("r", seq(4, 9), "BMI"), 
              paste0("r", seq(4, 9), "age_y_int")] <- 
        diag(x = 1, nrow = 6, ncol = 6)
      predict[paste0("r", seq(4, 9), "BMI"), 
              paste0("r", seq(4, 9), "shlt")] <- diag(x = 1, nrow = 6, ncol = 6)
      predict[paste0("r", seq(4, 9), "BMI"), 
              paste0("logr", seq(4, 9), "cesd")] <- 
        diag(x = 1, nrow = 6, ncol = 6)
      predict[paste0("r", seq(4, 9), "BMI"), "r4cesd_elevated"] <- 
        c(1, rep(0, 5))
      predict[paste0("r", seq(4, 9), "BMI"), "r9cesd_elevated"] <- 
        c(rep(0, 5), 1)
      
      #---- ****shlt models ----
      predict[paste0("r", seq(4, 9), "shlt"), 
              paste0("r", seq(4, 9), "BMI")] <- diag(x = 1, nrow = 6, ncol = 6)
      predict[paste0("r", seq(4, 9), "shlt"), 
              paste0("r", seq(4, 9), "age_y_int")] <- 
        diag(x = 1, nrow = 6, ncol = 6)
      predict[paste0("r", seq(4, 9), "shlt"), 
              paste0("r", seq(4, 9), "shlt")] <- 
        (diag(x = 1, nrow = 6, ncol = 6) == 0)*1
      predict[paste0("r", seq(4, 9), "shlt"), 
              paste0("logr", seq(4, 9), "cesd")] <- 
        diag(x = 1, nrow = 6, ncol = 6)
      predict[paste0("r", seq(4, 9), "shlt"), "r4cesd_elevated"] <- 
        c(1, rep(0, 5))
      predict[paste0("r", seq(4, 9), "shlt"), "r9cesd_elevated"] <- 
        c(rep(0, 5), 1)
      
      #---- ****Exposure models ----
      predict[c("r4cesd_elevated", "r9cesd_elevated", "total_elevated_cesd", 
              "avg_cesd", "avg_cesd_elevated"),  
              c("r4cesd_elevated", "r9cesd_elevated", "total_elevated_cesd", 
                "avg_cesd", "avg_cesd_elevated")] <- 
        diag(x = 0, nrow = 5, ncol = 5)
      
      #E1a
      predict["r4cesd_elevated", c(paste0("r", seq(5, 9), "BMI"), 
                                   paste0("r", seq(5, 9), "age_y_int"), 
                                   paste0("r", seq(5, 9), "shlt"), 
                                   paste0("logr", seq(5, 9), "cesd"))] <- 0
      
      #E1a
      predict["r9cesd_elevated", c(paste0("r", seq(4, 8), "BMI"), 
                                   paste0("r", seq(4, 8), "age_y_int"), 
                                   paste0("r", seq(4, 8), "shlt"), 
                                   paste0("logr", seq(4, 8), "cesd"))] <- 0
      
      #---- ***run imputation ----
      for(prop in mask_props){
        data <- get(paste0("mask", 100*prop))
        assign(paste0("jmvn_impute", 100*prop), 
               mice(data = data, m = num_impute, method = "norm", 
                    predictorMatrix = predict, where = is.na(data), 
                    blocks = as.list(paste0("logr", seq(4, 9), "cesd")), 
                    seed = 20210126))
        
        #---- ***save results ----
        if(save == "yes"){
          saveRDS(get(paste0("jmvn_impute", 100*prop)), 
                  file = here::here("MI datasets", 
                                    paste0("jmvn_", tolower(mechanism), 100*prop)))
        }
      }
    }
    
    
    # #---- **FCS ----
    # #---- ***predictor matrix ----
    # predict <- 
    #   matrix(1, nrow = 6, ncol = ncol(get(paste0("mask", 100*mask_props[1])))) %>% 
    #   set_rownames(paste0("r", seq(4, 9), "cesd")) %>% 
    #   set_colnames(colnames(get(paste0("mask", 100*mask_props[1]))))
    # #Don't use these as predictors
    # predict[, c("HHIDPN", "conde", "age_death_y", "r4cesd_elevated", 
    #             paste0("logr", seq(4, 9), "cesd"), "r9cesd_elevated", 
    #             "total_elevated_cesd", "avg_cesd", "avg_cesd_elevated", 
    #             "observed")] <- 0
    # 
    # #Use values at the current wave to predict-- not sure about this yet
    # predict[, paste0("r", seq(4, 9), "BMI")] <- diag(x = 1, nrow = 6, ncol = 6)
    # predict[, paste0("r", seq(4, 9), "age_y_int")] <- diag(x = 1, nrow = 6, 
    #                                                        ncol = 6)
    # predict[, paste0("r", seq(4, 9), "shlt")] <- diag(x = 1, nrow = 6, ncol = 6)
    # 
    # #Exclude values that predict themselves
    # predict[, paste0("r", seq(4, 9), "cesd")] <- 
    #   (diag(x = 1, nrow = 6, ncol = 6) == 0)*1
    # 
    # #---- ***run imputation ----
    # for(prop in mask_props){
    #   data <- get(paste0("mask", 100*prop))
    #   data %<>% mutate_at(paste0("r", seq(4, 9), "cesd"), as.factor)
    #   assign(paste0("fcs_impute", 100*prop), 
    #          mice(data = data, m = num_impute, method = "polr", 
    #               predictorMatrix = predict, where = is.na(data), 
    #               blocks = as.list(paste0("r", seq(4, 9), "cesd")), 
    #               seed = 20210126))
    #   
    #   #---- ***save results ----
    #   if(save == "yes"){
    #     saveRDS(get(paste0("fcs_impute", 100*prop)), 
    #             file = here::here("MI datasets", paste0("fcs_", tolower(mechanism), 
    #                                                     100*prop)))
    #   }
    # }
    # 
    # #---- **JMVN long ----
    # #Longitudinal joint multivariate normal model
    # 
    # #---- **FCS long ----
    # #Longitudinal fully conditional specification
    
    for(method in methods){
      for(mask in 100*mask_props){
        for(exposure in exposures){
          #---- fitting models ----
          if(method == "JMVN"){
            if(exposure == "CES-D Wave 4"){
              fitted_models <- 
                with(get(paste0(tolower(method), "_impute", mask)), 
                     coxph(Surv(survtime, observed) ~ r4age_y_int + female + 
                             hispanic + black + other + ed_cat + r4mstat_cat + 
                             ever_mem + ever_arthritis + ever_stroke + 
                             ever_heart + ever_lung + ever_cancer + ever_hibp + 
                             ever_diabetes + r4BMI + drinking4_cat_impute + 
                             smoker + 
                             ifelse(round((exp(logr4cesd) - 1), 0) > 4, 1, 0)))
            } else if(exposure == "CES-D Wave 9"){
              fitted_models <- 
                with(get(paste0(tolower(method), "_impute", mask)), 
                     coxph(Surv(survtime, observed) ~ r9age_y_int + female + 
                             hispanic + black + other + ed_cat + r9mstat_cat + 
                             ever_mem + ever_arthritis + ever_stroke + 
                             ever_heart + ever_lung + ever_cancer + ever_hibp + 
                             ever_diabetes + r9BMI + drinking9_cat_impute + 
                             smoker + 
                             ifelse(round((exp(logr9cesd) - 1), 0) > 4, 1, 0)))
            } else if(exposure == "Elevated CES-D Count"){
              fitted_models <- 
                with(get(paste0(tolower(method), "_impute", mask)), 
                     coxph(Surv(survtime, observed) ~ r4age_y_int + female + 
                             hispanic + black + other + ed_cat + r4mstat_cat + 
                             ever_mem + ever_arthritis + ever_stroke + 
                             ever_heart + ever_lung + ever_cancer + ever_hibp + 
                             ever_diabetes + r4BMI + drinking4_cat_impute + 
                             smoker + 
                             sum(ifelse(round((exp(logr4cesd) - 1), 0) > 4, 
                                        1, 0), 
                                 ifelse(round((exp(logr5cesd) - 1), 0) > 4, 
                                        1, 0), 
                                 ifelse(round((exp(logr6cesd) - 1), 0) > 4, 
                                        1, 0), 
                                 ifelse(round((exp(logr7cesd) - 1), 0) > 4, 
                                        1, 0), 
                                 ifelse(round((exp(logr8cesd) - 1), 0) > 4, 
                                        1, 0), 
                                 ifelse(round((exp(logr9cesd) - 1), 0) > 4, 
                                        1, 0)))) 
            } else{
              fitted_models <- 
                with(get(paste0(tolower(method), "_impute", mask)), 
                     coxph(Surv(survtime, observed) ~ r4age_y_int + female + 
                             hispanic + black + other + ed_cat + r4mstat_cat + 
                             ever_mem + ever_arthritis + ever_stroke + 
                             ever_heart + ever_lung + ever_cancer + ever_hibp + 
                             ever_diabetes + r4BMI + drinking4_cat_impute + 
                             smoker + 
                             ifelse(sum(ifelse(round((exp(logr4cesd) - 1), 
                                                     0) > 4, 1, 0), 
                                        ifelse(round((exp(logr5cesd) - 1), 
                                                     0) > 4, 1, 0), 
                                        ifelse(round((exp(logr6cesd) - 1), 
                                                     0) > 4, 1, 0), 
                                        ifelse(round((exp(logr7cesd) - 1), 
                                                     0) > 4, 1, 0), 
                                        ifelse(round((exp(logr8cesd) - 1), 
                                                     0) > 4, 1, 0), 
                                        ifelse(round((exp(logr9cesd) - 1), 
                                                     0) > 4, 1, 0))/6 > 4, 
                                    1, 0)))
            }
          }
          #---- pooling models ----
          pooled <- pool(fitted_models)
          
          #---- storing results ----
          pooled_effect_ests[which(
            pooled_effect_ests$Exposure == exposure & 
              pooled_effect_ests$Method == method & 
              pooled_effect_ests$Missingness == paste0(mask, "%")), 
            c("beta", "LCI", "UCI")] <- 
            summary(pooled, conf.int = TRUE, conf.level = 0.95, 
                    exponentiate = TRUE)[nrow(summary(pooled)), 
                                         c("estimate", "2.5 %", "97.5 %")]
        }
      }
    }
    
    #---- return values ----
    return(pooled_effect_ests)
  }

#---- testing ----
test <- mask_impute(CESD_data_wide, mechanism = "MCAR", mask_props, 
                    num_impute = 5, save = "no")
