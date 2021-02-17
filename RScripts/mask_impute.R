mask_impute <- 
  function(data_wide, mechanism, mask_props, num_impute, save = "no"){
    #---- create shell for data ----
    exposures <- c("CES-D Wave 4", "CES-D Wave 9", "Elevated CES-D Count", 
                   "Elevated Average CES-D")
    methods <- c("JMVN")
    
    pooled_effect_ests <- 
      data.frame("Exposure" = rep(exposures, length(mask_props)*length(methods)),
                 "beta" = NA, "LCI" = NA, "UCI" = NA, 
                 "Method" = rep(methods, 
                                each = length(exposures)*length(mask_props)), 
                 "Missingness" = 
                   rep(rep(paste0(mask_props*100, "%"), 
                           each = length(exposures)), length(methods)), 
                 "Type" = mechanism)
    
    #---- create incomplete data ----
    if(mechanism == "MCAR"){
      #---- **MCAR ----
      #it's easier to do this with my own code than the ampute function in MICE, 
      # which requires specifying all possible missing patterns you'd like it to 
      # consider
      data_long <- data_wide %>% 
        dplyr::select("HHIDPN", paste0("r", seq(4, 9), "cesd")) %>% 
        pivot_longer(-c("HHIDPN"), names_to = "wave", values_to = "cesd")
      
      for(prop in mask_props){
        assign(paste0(prop*100, "_mask"), 
               sample.int(n = nrow(data_long), 
                          size = floor(prop*nrow(data_long))))
      }
      
      for(mask in 100*mask_props){
        #mask these values in long data
        data_mask <- data_long
        data_mask[get(paste0(mask, "_mask")), "cesd"] <- NA
        data_mask %<>% na.omit() %>% 
          pivot_wider(names_from = "wave", values_from = "cesd")
        
        wide_data <- data_wide
        wide_data[which(wide_data$HHIDPN %in% data_mask$HHIDPN), 
                  paste0("r", seq(4, 9), "cesd")] <- data_mask[, -1] 
        wide_data[which(!wide_data$HHIDPN %in% data_mask$HHIDPN), 
                  paste0("r", seq(4, 9), "cesd")] <- NA
        
        assign(paste0("mask", mask), wide_data)
      }
    }
    
    #---- check missings ----
    #make sure no one is missing every cesd measure
    for(mask in 100*mask_props){
      data <- get(paste0("mask", mask))
      missing_counts <- 
        data %>% dplyr::select(paste0("r", seq(4, 9), "cesd")) %>% is.na() %>% 
        rowSums()
        
      data %<>% mutate("CESD_missing" =  missing_counts) %>% 
        filter(missing_counts < 6)
      
      assign(paste0("mask", mask), data)
    }
    
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
    
    #---- transformations ----
    #Taking the log for joint MNV models
    for(prop in mask_props){
      data <- get(paste0("mask", 100*prop))
      for(wave in 4:9){
        data[, paste0("logr", wave, "cesd")] <- 
          log(1 + data[, paste0("r", wave, "cesd")])
      }
      assign(paste0("mask", 100*prop), data)
    }
    
    #---- imputation ----
    num_impute <- num_impute
    
    #---- **JMVN ----
    #Joint multivariate normal
    #---- ***predictor matrix ----
    predict <- 
      matrix(1, nrow = 6, 
             ncol = ncol(get(paste0("mask", 100*mask_props[1])))) %>% 
      set_rownames(paste0("logr", seq(4, 9), "cesd")) %>% 
      set_colnames(colnames(get(paste0("mask", 100*mask_props[1]))))
    #Don't use these as predictors
    predict[, c("HHIDPN", "conde", "age_death_y", "r4cesd_elevated", 
                paste0("r", seq(4, 9), "cesd"), "r9cesd_elevated", 
                "total_elevated_cesd", "avg_cesd", "avg_cesd_elevated", 
                "observed")] <- 0
    
    #Use values at the current wave to predict
    predict[, paste0("r", seq(4, 9), "BMI")] <- diag(x = 1, nrow = 6, ncol = 6)
    predict[, paste0("r", seq(4, 9), "age_y_int")] <- diag(x = 1, nrow = 6, 
                                                           ncol = 6)
    predict[, paste0("r", seq(4, 9), "shlt")] <- diag(x = 1, nrow = 6, ncol = 6)
    
    #Exclude values that predict themselves
    predict[, paste0("logr", seq(4, 9), "cesd")] <- 
      (diag(x = 1, nrow = 6, ncol = 6) == 0)*1
    
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
          if(exposure == "CES-D Wave 4" & method == "JMVN"){
            fitted_models <- 
              with(get(paste0(tolower(method), "_impute", mask)), 
                   coxph(Surv(survtime, observed) ~ r4age_y_int + female + 
                           hispanic + black + other + ed_cat + r4mstat_cat + 
                           ever_mem + ever_arthritis + ever_stroke + 
                           ever_heart + ever_lung + ever_cancer + ever_hibp + 
                           ever_diabetes + r4BMI + drinking4_cat_impute + 
                           smoker + 
                           ifelse(round((exp(logr4cesd) - 1), 0) > 4, 1, 0)))
          } elseif(exposure == "CES-D Wave 9" & method == "JMVN"){
            fitted_models <- 
              with(get(paste0(tolower(method), "_impute", mask)), 
                   coxph(Surv(survtime, observed) ~ r9age_y_int + female + 
                           hispanic + black + other + ed_cat + r9mstat_cat + 
                           ever_mem + ever_arthritis + ever_stroke + 
                           ever_heart + ever_lung + ever_cancer + ever_hibp + 
                           ever_diabetes + r4BMI + drinking9_cat_impute + 
                           smoker + 
                           ifelse(round((exp(logr9cesd) - 1), 0) > 4, 1, 0)))
          } elseif(exposure == "Elevated CES-D Count" & method == "JMVN"){
            fitted_models <- 
              with(get(paste0(tolower(method), "_impute", mask)), 
                   coxph(Surv(survtime, observed) ~ r4age_y_int + female + 
                           hispanic + black + other + ed_cat + r4mstat_cat + 
                           ever_mem + ever_arthritis + ever_stroke + 
                           ever_heart + ever_lung + ever_cancer + ever_hibp + 
                           ever_diabetes + r4BMI + drinking4_cat_impute + 
                           smoker + 
                           ifelse(round((exp(logr9cesd) - 1), 0) > 4, 1, 0)))
          } elseif(exposure == "Elevated CES-D Count" & method == "JMVN"){
            fitted_models <- 
              with(get(paste0(tolower(method), "_impute", mask)), 
                   coxph(Surv(survtime, observed) ~ r4age_y_int + female + 
                           hispanic + black + other + ed_cat + r4mstat_cat + 
                           ever_mem + ever_arthritis + ever_stroke + 
                           ever_heart + ever_lung + ever_cancer + ever_hibp + 
                           ever_diabetes + r4BMI + drinking4_cat_impute + 
                           smoker + 
                           ifelse(round((exp(logr9cesd) - 1), 0) > 4, 1, 0)))
          }
          
        }
      }
    }
    
    #---- pooling models ----
    test_pool <- pool(test1)
    summary(test_pool, conf.int = TRUE, conf.level = 0.95, exponentiate = TRUE)
    
    
  }
