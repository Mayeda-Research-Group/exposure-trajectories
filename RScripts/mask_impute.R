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
  
  # #---- check missings ----
  #This code isn't going to work anymore b/c I changed the names of the datasets
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
    matrix(1, nrow = 6, ncol = ncol(get(paste0("mask", 100*mask_props[1])))) %>% 
    set_rownames(paste0("logr", seq(4, 9), "cesd")) %>% 
    set_colnames(colnames(get(paste0("mask", 100*mask_props[1]))))
  #Don't use these as predictors
  predict[, c("HHIDPN", "conde", "age_death_y", "r4cesd_elevated", 
              paste0("r", seq(4, 9), "cesd"), "r9cesd_elevated", 
              "total_elevated_cesd", "avg_cesd", "avg_cesd_elevated", 
              "observed")] <- 0
  
  #Use values at the current wave to predict-- not sure about this yet
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
    assign(paste0("impute", 100*prop), 
           mice(data = data, m = num_impute, method = "norm", 
                predictorMatrix = predict, where = is.na(data), 
                blocks = as.list(paste0("logr", seq(4, 9), "cesd")), 
                seed = 20210126))
    
    #---- ***save results ----
    saveRDS(get(paste0("impute", 100*prop)), 
            file = here::here("MI datasets", 
                              paste0("jmvn_", tolower(mechanism), 100*prop)))
  }
  
  #---- **FCS ----
  #---- ***predictor matrix ----
  predict <- 
    matrix(1, nrow = 6, ncol = ncol(get(paste0("mask", 100*mask_props[1])))) %>% 
    set_rownames(paste0("r", seq(4, 9), "cesd")) %>% 
    set_colnames(colnames(get(paste0("mask", 100*mask_props[1]))))
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
    data <- get(paste0("mask", 100*prop))
    data %<>% mutate_at(paste0("r", seq(4, 9), "cesd"), as.factor)
    assign(paste0("impute", 100*prop), 
           mice(data = data, m = num_impute, method = "polr", 
                predictorMatrix = predict, where = is.na(data), 
                blocks = as.list(paste0("r", seq(4, 9), "cesd")), 
                seed = 20210126))
    
    #---- ***save results ----
    saveRDS(get(paste0("impute", 100*prop)), 
            file = here::here("MI datasets", paste0("fcs_", tolower(mechanism), 
                                                    100*prop)))
  }
  
  #---- **JMVN long ----
  #Longitudinal joint multivariate normal model
  
  #---- **FCS long ----
  #Longitudinal fully conditional specification
  
}
