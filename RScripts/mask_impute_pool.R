mask_impute_pool <- 
  function(data_wide, exposures, mechanism, method, mask_percent, num_impute, 
           save = "no"){
    #---- create shell for output ----
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
    
    #---- imputation ----
    #---- **predictor matrix ----
    time_updated_vars <- c("mstat_impute", "drinking_impute", "memrye_impute", 
                           "stroke_impute", "hearte_impute", "lunge_impute", 
                           "cancre_impute", "hibpe_impute", "diabe_impute", 
                           "cesd", "BMI", "shlt")
    blocks <- c(apply(expand.grid("r", seq(4, 9), time_updated_vars), 1, 
                      paste, collapse = ""),
                "r4cesd_elevated", "r9cesd_elevated", 
                "total_elevated_cesd", "avg_cesd", 
                "avg_cesd_elevated")
    predict <- matrix(1, length(blocks), ncol = ncol(data_wide)) %>% 
      set_rownames(blocks) %>% 
      set_colnames(colnames(data_wide))
    
    #Don't use these as predictors
    predict[, c("HHIDPN", paste0("r", seq(4, 9), "mstat_cat"), "r3cesd",
                paste0("r", seq(4, 9), "drinking_cat"), 
                paste0("r", seq(3, 9), "conde_impute"), "white", "age_death_y", 
                "observed", "CESD_missing")] <- 0
    
    #---- ****time-updated var models ----
    for(var in time_updated_vars){
      #can't predict itself
      predict[paste0("r", seq(4, 9), var), paste0("r", seq(4, 9), var)] <- 
        (diag(x = 1, nrow = 6, ncol = 6) == 0)*1
      
      #use time-updated predictors
      predictors <- c(time_updated_vars[which(time_updated_vars != var)], 
                      "age_y_int")
      for(predictor in predictors){
        predict[paste0("r", seq(4, 9), var), 
                paste0("r", seq(4, 9), predictor)] <- 
          diag(x = 1, nrow = 6, ncol = 6)
      }
      
      #use relevant wave-specific exposures
      predict[paste0("r", seq(4, 9), var), "r4cesd_elevated"] <- c(1, rep(0, 5))
      predict[paste0("r", seq(4, 9), var), "r9cesd_elevated"] <- c(rep(0, 5), 1)
    }
    
    #---- ****exposure models ----
    #Can't predict 
    predict[c("r4cesd_elevated", "r9cesd_elevated", "total_elevated_cesd", 
              "avg_cesd", "avg_cesd_elevated"),  
            c("r4cesd_elevated", "r9cesd_elevated", "total_elevated_cesd", 
              "avg_cesd", "avg_cesd_elevated")] <- 
      diag(x = 0, nrow = 5, ncol = 5)
    
    #E1a-- can only use wave 4 data
    predictors <- c(time_updated_vars, "age_y_int")
    predict["r4cesd_elevated", 
            apply(expand.grid("r", seq(5, 9), predictors), 1, 
                  paste, collapse = "")] <- 0
    
    #E1b-- can only use wave 9 data
    predict["r9cesd_elevated", 
            apply(expand.grid("r", seq(4, 8), predictors), 1, 
                  paste, collapse = "")] <- 0
    
    #---- **run imputation ----
    if(method == "JMVN"){
      #Joint multivariate normal
      data_imputed <- mice(data = data_wide, m = num_impute, method = "norm", 
                           predictorMatrix = predict, where = is.na(data_wide), 
                           blocks = as.list(rownames(predict)), 
                           seed = 20210126)
    } else if(method == "FCS"){
      #Fully conditional specification
      data_imputed <- mice(data = data_wide, m = num_impute, 
                           method = "polr", 
                           predictorMatrix = predict, where = is.na(data_wide), 
                           blocks = as.list(rownames(predict)), 
                           seed = 20210126)
    }
    
    #---- **save results ----
    if(save == "yes"){
      saveRDS(data_imputed, 
              file = here::here("MI datasets", 
                                paste0(tolower(method), "_", tolower(mechanism), 
                                       as.numeric(sub("%","", mask_percent)))))
    }
    
    for(exposure in exposures){
      #---- fitted models ----
      if(exposure == "CES-D Wave 4"){
        fitted_models <- 
          with(data_imputed, 
               coxph(Surv(survtime, observed) ~ r4age_y_int + female + 
                       hispanic + black + other + ed_cat + r4mstat_cat + 
                       ever_mem + ever_arthritis + ever_stroke + 
                       ever_heart + ever_lung + ever_cancer + ever_hibp + 
                       ever_diabetes + r4BMI + drinking4_cat_impute + 
                       smoker + r4cesd_elevated))
      } else if(exposure == "CES-D Wave 9"){
        fitted_models <- 
          with(data_imputed, 
               coxph(Surv(survtime, observed) ~ r9age_y_int + female + 
                       hispanic + black + other + ed_cat + r9mstat_cat + 
                       ever_mem + ever_arthritis + ever_stroke + 
                       ever_heart + ever_lung + ever_cancer + ever_hibp + 
                       ever_diabetes + r9BMI + drinking9_cat_impute + 
                       smoker + r9cesd_elevated))
      } else if(exposure == "Elevated CES-D Count"){
        fitted_models <- 
          with(data_imputed, 
               coxph(Surv(survtime, observed) ~ r4age_y_int + female + 
                       hispanic + black + other + ed_cat + r4mstat_cat + 
                       ever_mem + ever_arthritis + ever_stroke + 
                       ever_heart + ever_lung + ever_cancer + ever_hibp + 
                       ever_diabetes + r4BMI + drinking4_cat_impute + 
                       smoker + total_elevated_cesd))
      } else{
        fitted_models <- 
          with(data_imputed, 
               coxph(Surv(survtime, observed) ~ r4age_y_int + female + 
                       hispanic + black + other + ed_cat + r4mstat_cat + 
                       ever_mem + ever_arthritis + ever_stroke + 
                       ever_heart + ever_lung + ever_cancer + ever_hibp + 
                       ever_diabetes + r4BMI + drinking4_cat_impute + 
                       smoker + avg_cesd_elevated))
      }
      
      #---- pooling models ----
      pooled <- pool(fitted_models)
      
      #---- storing results ----
      pooled_effect_ests[which(pooled_effect_ests$Exposure == exposure), 
                         c("beta", "LCI", "UCI")] <- 
        summary(pooled, conf.int = TRUE, conf.level = 0.95, 
                exponentiate = TRUE)[nrow(summary(pooled)), 
                                     c("estimate", "2.5 %", "97.5 %")]
    }
    
    #---- return values ----
    return(pooled_effect_ests)
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

# #---- testing ----
# #Single run
# test <- mask_impute_pool(CESD_data_wide, mechanism = "MCAR", method = "FCS",
#                          mask_percent = "10%", num_impute = 5, save = "no")
# #Multiple runs
# test_2 <- replicate(2, mask_impute_pool(CESD_data_wide,
#                                                mechanism = "MCAR",
#                                                method = "JMVN",
#                                                mask_percent = "10%",
#                                                num_impute = 5, save = "no"),
#                     simplify = FALSE)
# 
# #Formatting data
# formatted <- do.call(rbind, test_2)
# 
# #Summarizing results
# results <- formatted %>% group_by(Exposure) %>%
#   summarize_at(.vars = c("beta", "LCI", "UCI"), .funs = mean)

