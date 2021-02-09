mask_impute <- function(data_wide, mechanism, mask_props, num_impute){
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
      assign(paste0("mcar", prop*100, "_mask"), 
             sample.int(n = nrow(data_long), 
                        size = floor(prop*nrow(data_long))))
    }
    
    for(mask in 100*mask_props){
      #mask these values in long data
      data_mask <- data_long
      data_mask[get(paste0("mcar", mask, "_mask")), "cesd"] <- NA
      data_mask %<>% na.omit() %>% 
        pivot_wider(names_from = "wave", values_from = "cesd")
      
      wide_data <- data_wide
      wide_data[which(wide_data$HHIDPN %in% data_mask$HHIDPN), 
                paste0("r", seq(4, 9), "cesd")] <- data_mask[, -1] 
      wide_data[which(!wide_data$HHIDPN %in% data_mask$HHIDPN), 
                paste0("r", seq(4, 9), "cesd")] <- NA
      
      assign(paste0("mcar", mask), wide_data)
    }
  }
  
  # #---- check missings ----
  # #make sure no one is missing every cesd measure
  # missings <- as.data.frame(matrix(nrow = length(mask_props), ncol = 8)) %>% 
  #   set_colnames(c("Mask Prop", seq(0, 6)))
  # missings[, "Mask Prop"] <- 100*mask_props
  # 
  # for(mask in 100*mask_props){
  #   data <- get(paste0("mcar", mask)) %>% 
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
    data <- get(paste0(tolower(mechanism), 100*prop))
    for(wave in 4:9){
      data[, paste0("logr", wave, "cesd")] <- 
        log(1 + data[, paste0("r", wave, "cesd")])
    }
    assign(paste0(tolower(mechanism), 100*prop), data)
  }
  
  #---- imputation ----
  num_impute <- num_impute
  
  #---- **JMVN ----
  #Joint multivariate normal
  #---- ***predictor matrix ----
  predict <- matrix(1, nrow = 6, ncol = ncol(mcar10)) %>% 
    set_rownames(paste0("logr", seq(4, 9), "cesd")) %>% 
    set_colnames(colnames(mcar10))
  #Don't use these as predictors
  predict[, c("HHIDPN", "conde", "age_death_y", "r4cesd_elevated", 
              paste0("r", seq(4, 9), "cesd"), "r9cesd_elevated", 
              "total_elevated_cesd", "avg_cesd", "avg_cesd_elevated", 
              "observed")] <- 0
  
  # #Use values at the current wave to predict-- not sure about this yet
  predict[, paste0("r", seq(4, 9), "BMI")] <- diag(x = 1, nrow = 6, ncol = 6)
  predict[, paste0("r", seq(4, 9), "age_y_int")] <- diag(x = 1, nrow = 6, 
                                                         ncol = 6)
  predict[, paste0("r", seq(4, 9), "shlt")] <- diag(x = 1, nrow = 6, ncol = 6)
  
  #Exclude values that predict themselves
  predict[, paste0("logr", seq(4, 9), "cesd")] <- 
    (diag(x = 1, nrow = 6, ncol = 6) == 0)*1
  
  #---- ***run imputation ----
  for(prop in mask_props){
    data <- get(paste0("mcar", 100*prop))
    assign(paste0("jmvn_mcar", 100*prop), 
           mice(data = data, m = num_impute, method = "norm", 
                predictorMatrix = predict, where = is.na(data), 
                blocks = as.list(paste0("logr", seq(4, 9), "cesd")), 
                seed = 20210126))
    
    #---- ***save results ----
    saveRDS(get(paste0("jmvn_mcar", 100*prop)), 
            file = here("MI datasets", paste0("jmvn_mcar", 100*prop)))
  }
  
  #---- **FCS ----
  #---- ***predictor matrix ----
  predict <- matrix(1, nrow = 6, ncol = ncol(mcar10)) %>% 
    set_rownames(paste0("r", seq(4, 9), "cesd")) %>% 
    set_colnames(colnames(mcar10))
  #Don't use these as predictors
  predict[, c("HHIDPN", "conde", "age_death_y", "r4cesd_elevated", 
              paste0("logr", seq(4, 9), "cesd"), "r9cesd_elevated", 
              "total_elevated_cesd", "avg_cesd", "avg_cesd_elevated", 
              "observed")] <- 0
  
  #Use values at the current wave to predict-- not sure about this yet
  predict[, paste0("r", seq(4, 9), "BMI")] <- diag(x = 1, nrow = 6, ncol = 6)
  predict[, paste0("r", seq(4, 9), "age_y_int")] <- diag(x = 1, nrow = 6, 
                                                         ncol = 6)
  predict[, paste0("r", seq(4, 9), "shlt")] <- diag(x = 1, nrow = 6, ncol = 6)
  
  #Exclude values that predict themselves
  predict[, paste0("r", seq(4, 9), "cesd")] <- 
    (diag(x = 1, nrow = 6, ncol = 6) == 0)*1
  
  #---- ***run imputation ----
  for(prop in mask_props){
    data <- get(paste0("mcar", 100*prop))
    data %<>% mutate_at(paste0("r", seq(4, 9), "cesd"), as.factor)
    assign(paste0("fcs_mcar", 100*prop), 
           mice(data = data, m = num_impute, method = "polr", 
                predictorMatrix = predict, where = is.na(data), 
                blocks = as.list(paste0("r", seq(4, 9), "cesd")), 
                seed = 20210126))
    
    #---- ***save results ----
    saveRDS(get(paste0("fcs_mcar", 100*prop)), 
            file = here("MI datasets", paste0("fcs_mcar", 100*prop)))
  }
  
  #---- **JMVN long ----
  #Longitudinal joint multivariate normal model
  
  #---- **FCS long ----
  #Longitudinal fully conditional specification
  
  #---- effect estimates ----
  #Based on imputations
  model_list <- 
    lapply(model_list <- vector(mode = "list", length(methods)),
           function(x) x <- lapply(x <- vector(mode = "list", length(mask_props)), 
                                   function(x) x <- 
                                     lapply(x <- vector(mode = "list", 
                                                        length = 4), 
                                            function(x) x <- vector(mode = "list", 
                                                                    length = num_impute))))
  
  #naming layers of list
  names(model_list) <- methods
  for(i in 1:length(model_list)){
    names(model_list[[i]]) <- 100*mask_props
    for(j in 1:length(model_list[[i]])){
      names(model_list[[i]][[j]]) <- exposures
    }
  }
  
  for(m in methods[1:2]){
    for(i in 1:length(mask_props)){
      for(j in 1:num_impute){
        imputations <- get(paste0(tolower(m), "_mcar", mask_props[i]*100))
        imputed_data = complete(imputations, action = j)
        if(m == "JMVN"){
          #transform back to original values
          imputed_data[, paste0("r", seq(4, 9), "cesd")] <- 
            round(exp(imputed_data[, paste0("logr", seq(4, 9), "cesd")]) - 1)
        }
        
        if(m == "FCS"){
          imputed_data %<>% mutate_at(paste0("r", seq(4, 9), "cesd"), as.numeric)
        }
        
        #---- **E1a; E1b; E3 ----
        imputed_data %<>% 
          mutate("r4cesd_elevated" = ifelse(r4cesd > 4, 1, 0), 
                 "r9cesd_elevated" = ifelse(r9cesd > 4, 1, 0), 
                 "avg_cesd" = imputed_data %>% 
                   dplyr::select(paste0("r", seq(4, 9, by = 1), "cesd")) %>% 
                   rowMeans(), 
                 "avg_cesd_elevated" = ifelse(avg_cesd > 4, 1, 0))
        
        #---- **E2 ----
        elevated_cesd <- imputed_data %>% 
          dplyr::select(paste0("r", seq(4, 9, by = 1), "cesd"))
        
        elevated_cesd <- (elevated_cesd > 4)*1
        
        imputed_data %<>% mutate("total_elevated_cesd" = rowSums(elevated_cesd))
        
        model_list[[m]][[i]][["CES-D Wave 4"]][[j]] <- 
          coxph(Surv(survtime, observed) ~ r4age_y_int + female + hispanic + 
                  black + other + ed_cat + r4mstat_cat + ever_mem + 
                  ever_arthritis + ever_stroke + ever_heart + ever_lung + 
                  ever_cancer + ever_hibp + ever_diabetes + r4BMI + 
                  drinking4_cat_impute + smoker + r4cesd_elevated, 
                data = imputed_data)
        
        model_list[[m]][[i]][["CES-D Wave 9"]][[j]] <- 
          coxph(Surv(survtime, observed) ~ r9age_y_int + female + hispanic + 
                  black + other + ed_cat + r9mstat_cat + ever_mem + 
                  ever_arthritis + ever_stroke + ever_heart + ever_lung + 
                  ever_cancer + ever_hibp + ever_diabetes + r9BMI + 
                  drinking9_cat_impute + smoker + r9cesd_elevated, 
                data = imputed_data)
        
        model_list[[m]][[i]][["Elevated CES-D Count"]][[j]] <- 
          coxph(Surv(survtime, observed) ~ r4age_y_int + female + hispanic + 
                  black + other + ed_cat + r4mstat_cat + ever_mem + 
                  ever_arthritis + ever_stroke + ever_heart + ever_lung + 
                  ever_cancer + ever_hibp + ever_diabetes + r4BMI + 
                  drinking4_cat_impute + smoker + total_elevated_cesd, 
                data = imputed_data)
        
        model_list[[m]][[i]][["Elevated Average CES-D"]][[j]] <- 
          coxph(Surv(survtime, observed) ~ r4age_y_int + female + hispanic + 
                  black + other + ed_cat + r4mstat_cat + ever_mem + 
                  ever_arthritis + ever_stroke + ever_heart + ever_lung + 
                  ever_cancer + ever_hibp + ever_diabetes + r4BMI + 
                  drinking4_cat_impute + smoker + avg_cesd_elevated, 
                data = imputed_data)
      }
    }
  }
  
  pooled_model_list <- 
    lapply(pooled_model_list <- vector(mode = "list", length(methods)),
           function(x) x <- lapply(x <- vector(mode = "list", length(mask_props)), 
                                   function(x) x <- 
                                     vector(mode = "list", length = 4)))
  
  #naming layers of list
  names(pooled_model_list) <- methods
  for(i in 1:length(pooled_model_list)){
    names(pooled_model_list[[i]]) <- 100*mask_props
    for(j in 1:length(pooled_model_list[[i]])){
      names(pooled_model_list[[i]][[j]]) <- exposures
    }
  }
  
  for(m in methods[1:2]){
    for(i in as.character(mask_props*100)){
      for(j in exposures){
        pooled_model_list[[m]][[i]][[j]] <- 
          summary(pool(model_list[[m]][[i]][[j]]))
      }
    } 
  }
  
  for(m in methods[1:2]){
    for(prop in as.character(100*mask_props)){
      for(exposure in exposures){
        last_term <- nrow(pooled_model_list[[m]][[prop]][[exposure]])
        beta <- 
          exp(pooled_model_list[[m]][[prop]][[exposure]]$estimate[last_term])
        UCI <- 
          exp(pooled_model_list[[m]][[prop]][[exposure]]$estimate[last_term] + 
                pooled_model_list[[m]][[prop]][[exposure]]$std.error[last_term])
        LCI <- 
          exp(pooled_model_list[[m]][[prop]][[exposure]]$estimate[last_term] - 
                pooled_model_list[[m]][[prop]][[exposure]]$std.error[last_term])
        
        table_effect_ests[which(table_effect_ests$Exposure == exposure & 
                                  table_effect_ests$Method == m & 
                                  table_effect_ests$Missingness == 
                                  paste0(prop, "%")), 
                          c("beta", "LCI", "UCI")] <- c(beta, LCI, UCI)
      }
    }
  }
  
}