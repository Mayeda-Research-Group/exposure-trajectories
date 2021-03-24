mask_impute_pool <- 
  function(data_wide, exposures, mechanism, method, mask_percent, num_impute, 
           save = "no"){
    #---- create shell for output ----
    pooled_effect_ests <- 
      data.frame("Exposure" = exposures, "beta" = NA, "SD" = NA, "LCI" = NA, 
                 "UCI" = NA, "Method" = method, "Missingness" = mask_percent, 
                 "Type" = mechanism)
    
    #---- create incomplete data ----
    if(mechanism == "MCAR"){
      #---- **MCAR ----
      #it's easier to do this with my own code than the ampute function in MICE, 
      # which requires specifying all possible missing patterns you'd like it to 
      # consider
      mask_prop <- as.numeric(sub("%","", mask_percent))/100
      total_indices <- nrow(data_wide)*6 #6 waves of data per person
      mask_index <- sample(seq(1, total_indices), 
                           size = floor(mask_prop*total_indices), 
                           replace = FALSE)
    }
    
    #masking wave-specific values
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
    #make sure no one is missing every cesd measure
    data_wide %<>% 
      mutate("CESD_missing" = 
               rowSums(data_wide %>% 
                         dplyr::select(paste0("r", seq(4, 9), "cesd")) %>% 
                         mutate_all(function(x) is.na(x)))) %>% 
      filter(CESD_missing < 6)
    
    time_updated_vars <- c("married_partnered", "not_married_partnered", 
                           "widowed", "drinking_cat", "memrye_impute", 
                           "stroke_impute", "hearte_impute", "lunge_impute", 
                           "cancre_impute", "hibpe_impute", "diabe_impute", 
                           "cesd", "BMI", "shlt")
    
    time_invariant_vars <- c("ed_cat", "white", "black", "hispanic", "other", 
                             "female", "survtime", "death2018", "smoker")
    
    #---- convert to long?? ----
    if(method %in% c("2l.norm", "2l.fcs")){
      data_long <- data_wide %>% 
        dplyr::select("HHIDPN", 
                      apply(expand.grid("r", seq(4, 9), 
                                        c(time_updated_vars, "age_y_int")), 1, 
                            paste, collapse = "")) %>% 
        pivot_longer(-c("HHIDPN"), names_to = c("wave", ".value"), 
                     names_pattern = "(^[a-zA-Z][0-9])(.*)")
      
      data_long %<>% 
        left_join(., data_wide %>% 
                    dplyr::select("HHIDPN", all_of(time_invariant_vars))) %>% 
        mutate_at("HHIDPN", as.integer)
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
    
    #---- imputation ----
    #---- **predictor matrix ----
    if(method %in% c("2l.norm", "2l.fcs")){
      blocks <- time_updated_vars
      predict <- matrix(1, length(blocks), ncol = ncol(data_long)) %>% 
        set_rownames(blocks) %>% 
        set_colnames(colnames(data_long))
      
      #Don't use these as predictors
      predict[, c("wave")] <- 0
      
      #Indicated cluster variable
      predict[, "HHIDPN"] <- -2
      
      #Can't predict themselves
      predict[time_updated_vars, time_updated_vars] <- 
        (diag(x = 1, nrow = length(time_updated_vars), 
              ncol = length(time_updated_vars)) == 0)*1
      
    } else{
      blocks <- c(apply(expand.grid("r", seq(4, 9), time_updated_vars), 1, 
                        paste, collapse = ""))
      predict <- matrix(1, length(blocks), ncol = ncol(data_wide)) %>% 
        set_rownames(blocks) %>% 
        set_colnames(colnames(data_wide))
      
      #Don't use these as predictors
      predict[, c("HHIDPN", paste0("r", seq(3, 9), "conde_impute"), 
                  "age_death_y", "observed", "CESD_missing", "r3cesd", 
                  "r4cesd_elevated", "r9cesd_elevated",  "avg_cesd", 
                  "avg_cesd_elevated", "total_elevated_cesd")] <- 0
      
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
      }
    }
    
    #---- **run imputation ----
    max_it <- tibble("Method" = c("FCS", "JMVN", "PMM", "2l.norm", "2l.fcs"), 
                     "10%" = c(5, 5, 5, 5, 5),
                     "20%" = c(5, 5, 5, 5, 5),
                     "30%" = c(5, 5, 5, 5, 5)) %>% 
      column_to_rownames("Method")
    
    #---- ****JMVN ----
    if(method == "JMVN"){
      #Joint multivariate normal
      data_imputed <- mice(data = data_wide, m = num_impute, 
                           maxit = max_it[method, mask_percent], 
                           method = "norm", predictorMatrix = predict, 
                           where = is.na(data_wide), 
                           blocks = as.list(rownames(predict)), seed = 20210126)
      
      # #look at convergence
      #   #10% missing needs maxit = 20
      #   #20% missing needs maxit = 30
      #   #30% missing needs maxit = 40
      # plot(data_imputed)
      
    } else if(method == "FCS"){
      #---- ****FCS ----
      #Fully conditional specification
      data_wide %<>% 
        mutate_at(vars(c(paste0("r", seq(4, 9), "mstat_impute"), "ed_cat",
                         paste0("r", seq(4, 9), "memrye_impute"), 
                         paste0("r", seq(4, 9), "stroke_impute"),
                         paste0("r", seq(4, 9), "hearte_impute"),
                         paste0("r", seq(4, 9), "lunge_impute"), 
                         paste0("r", seq(4, 9), "cancre_impute"), 
                         paste0("r", seq(4, 9), "hibpe_impute"), 
                         paste0("r", seq(4, 9), "diabe_impute"), "smoker", 
                         "hispanic", "black", "other", "female", "death2018")), 
                  as.factor)
      
      #start <- Sys.time()
      data_imputed <- mice(data = data_wide, m = num_impute, 
                           maxit = max_it[method, mask_percent],
                           nnet.MaxNWts = 5000,
                           defaultMethod = 
                             c("norm", "logreg", "polyreg", "polr"),
                           predictorMatrix = predict, where = is.na(data_wide), 
                           blocks = as.list(rownames(predict)),
                           seed = 20210126)
      #end <- Sys.time() - start
      
      # #look at convergence
      #   #10% missing needs maxit = 20
      #   #20% missing needs maxit = 40
      #   #30% missing needs maxit = 40
      #plot(data_imputed)
      
    } else if(method == "PMM"){
      #---- ****PMM ----
      #Predictive Mean Matching
      data_wide %<>% 
        mutate_at(vars(c(paste0("r", seq(4, 9), "mstat_impute"), "ed_cat",
                         paste0("r", seq(4, 9), "memrye_impute"), 
                         paste0("r", seq(4, 9), "stroke_impute"),
                         paste0("r", seq(4, 9), "hearte_impute"),
                         paste0("r", seq(4, 9), "lunge_impute"), 
                         paste0("r", seq(4, 9), "cancre_impute"), 
                         paste0("r", seq(4, 9), "hibpe_impute"), 
                         paste0("r", seq(4, 9), "diabe_impute"), "smoker", 
                         "hispanic", "black", "other", "female", "death2018")), 
                  as.factor)
      
      #start <- Sys.time()
      data_imputed <- mice(data = data_wide, m = num_impute, 
                           maxit = max_it[method, mask_percent], 
                           method = "pmm", predictorMatrix = predict, 
                           where = is.na(data_wide), 
                           blocks = as.list(rownames(predict)), 
                           seed = 20210126)
      #stop <- Sys.time() - start
      
      # #look at convergence
      # #10% missing needs maxit = 40
      # #20% missing needs maxit = 40
      # #30% missing needs maxit = 40
      #  plot(data_imputed)
      
    } else if(method == "2l.norm"){
      #---- ****2l.norm ----
      #2-level heteroskedatic between group variances
      start <- Sys.time()
      data_imputed <- mice(data = data_long, m = num_impute, 
                           #maxit = max_it[method, mask_percent],
                           maxit = 20,
                           method = "2l.norm", predictorMatrix = predict, 
                           where = is.na(data_long), 
                           blocks = as.list(rownames(predict)), 
                           seed = 20210126)
      stop <- Sys.time() - start
      
      #look at convergence
      #10% missing needs maxit = 20
      #20% missing needs maxit = 
      #30% missing needs maxit = 
      plot(data_imputed)
      
    } else if(method == "2l.fcs"){
      #---- ****2l.fcs ----
      #Longitudinal fully conditional specification
      #Fully conditional specification
      data_long %<>% 
        mutate_at(vars(c("mstat_impute", "memrye_impute", "stroke_impute", 
                         "hearte_impute", "lunge_impute", "cancre_impute", 
                         "hibpe_impute", "diabe_impute", "ed_cat", "black", 
                         "hispanic", "other", "female", "death2018", "smoker")), 
                  as.factor)
      
      start <- Sys.time()
      data_imputed <- mice(data = data_long, m = num_impute, 
                           maxit = 20, 
                           defaultMethod = 
                             c("2l.norm", "2l.bin", "2l.norm", "2l.norm"), 
                           predictorMatrix = predict, 
                           where = is.na(data_long), 
                           blocks = as.list(rownames(predict)), 
                           seed = 20210126)
      stop <- Sys.time() - start
      
      # #look at convergence
      # #10% missing needs maxit = 
      # #20% missing needs maxit = 
      # #30% missing needs maxit = 
      plot(data_imputed) 
      
    } 
    
    #---- **save results ----
    if(save == "yes"){
      saveRDS(data_imputed, 
              file = here::here("MI datasets", 
                                paste0(tolower(method), "_", tolower(mechanism), 
                                       as.numeric(sub("%","", mask_percent)))))
    }
    
    #---- fitted models ----
    model_list <- vector(mode = "list", length = length(exposures)) %>% 
      set_names(exposures)
    
    for(i in 1:num_impute){
      complete_data <- complete(data_imputed, action = i)
      
      #---- **post process: dummy vars ----
      for(wave in seq(4, 9)){
        vars <- c("married_partnered", "not_married_partnered", "widowed")
        cols <- apply(expand.grid("r", wave, vars), 1, paste, collapse = "")
        subset <- complete_data[, cols]
        rowmax <- apply(subset, 1, function(x) max(x))
        subset <- subset/rowmax
        subset[subset < 1] <- 0
        subset[subset > 1] <- 1
        
        complete_data[, colnames(subset)] <- subset
      }
      
      #---- **post process: binary vars ----
      waves <- seq(4, 9)
      vars <- c("memrye_impute", "stroke_impute", "hearte_impute", 
                "lunge_impute", "cancre_impute", "hibpe_impute", "diabe_impute")
      cols <- apply(expand.grid("r", waves, vars), 1, paste, collapse = "")
      
      #fix impossible probs
      subset <- complete_data[, cols]
      subset[subset < 0] <- 0
      subset[subset > 1] <- 1
      
      for(col in cols){
        complete_data[, col] <- 
          rbinom(n = nrow(complete_data), size = 1, prob = subset[, col])
      }
      
      #---- **post-process: exposures ----
      complete_data %<>% 
        mutate("r4cesd_elevated" = ifelse(r4cesd > 4, 1, 0), 
               "r9cesd_elevated" = ifelse(r9cesd > 4, 1, 0), 
               "total_elevated_cesd" = 
                 complete_data %>% 
                 dplyr::select(paste0("r", seq(4, 9), "cesd")) %>% 
                 mutate_all(function(x) ifelse(x > 4, 1, 0)) %>% rowSums(), 
               "avg_cesd" = complete_data %>% 
                 dplyr::select(paste0("r", seq(4, 9), "cesd")) %>% rowMeans(), 
               "avg_cesd_elevated" = ifelse(avg_cesd > 4, 1, 0))
      
      # #Sanity check
      # View(complete_data %>% dplyr::select(contains("cesd")))
      
      for(exposure in exposures){
        if(exposure == "CES-D Wave 4"){
          model_list[[exposure]][[i]] <- 
            with(complete_data, 
                 coxph(Surv(survtime, observed) ~ r4not_married_partnered + 
                         r4widowed + ed_cat + r4drinking_cat + r4memrye_impute + 
                         r4stroke_impute + r4hearte_impute + r4lunge_impute + 
                         r4cancre_impute + r4hibpe_impute + r4diabe_impute + 
                         smoker + r4BMI + hispanic + black + other + female + 
                         r4age_y_int + r4shlt + r4cesd_elevated))
        } else if(exposure == "CES-D Wave 9"){
          model_list[[exposure]][[i]] <- 
            with(complete_data, 
                 coxph(Surv(survtime, observed) ~ r9not_married_partnered + 
                         r9widowed + ed_cat + r9drinking_cat + r9memrye_impute + 
                         r9stroke_impute + r9hearte_impute + r9lunge_impute + 
                         r9cancre_impute + r9hibpe_impute + r9diabe_impute + 
                         smoker + r9BMI + hispanic + black + other + female + 
                         r9age_y_int + r9shlt + r9cesd_elevated))
        } else if(exposure == "Elevated CES-D Count"){
          model_list[[exposure]][[i]] <- 
            with(complete_data, 
                 coxph(Surv(survtime, observed) ~ r4not_married_partnered + 
                         r4widowed + ed_cat + r4drinking_cat + r4memrye_impute + 
                         r4stroke_impute + r4hearte_impute + r4lunge_impute + 
                         r4cancre_impute + r4hibpe_impute + r4diabe_impute + 
                         smoker + r4BMI + hispanic + black + other + female + 
                         r4age_y_int + r4shlt + total_elevated_cesd))
        } else{
          model_list[[exposure]][[i]] <- 
            with(complete_data, 
                 coxph(Surv(survtime, observed) ~ r4not_married_partnered + 
                         r4widowed + ed_cat + r4drinking_cat + r4memrye_impute + 
                         r4stroke_impute + r4hearte_impute + r4lunge_impute + 
                         r4cancre_impute + r4hibpe_impute + r4diabe_impute + 
                         smoker + r4BMI + hispanic + black + other + female + 
                         r4age_y_int + r4shlt + avg_cesd_elevated))
        }
      }
    }
    
    #---- pooling models ----
    for(exposure in exposures){
      pooled_model <- summary(pool(model_list[[exposure]]))
      pooled_effect_ests[which(pooled_effect_ests$Exposure == exposure), 
                         c("beta", "SD")] <- 
        pooled_model[nrow(pooled_model), c("estimate", "std.error")]
    }
    
    pooled_effect_ests[, "LCI"] <- 
      pooled_effect_ests$beta - 1.96*pooled_effect_ests$SD
    
    pooled_effect_ests[, "UCI"] <- 
      pooled_effect_ests$beta + 1.96*pooled_effect_ests$SD
    
    #---- return values ----
    return(pooled_effect_ests)
  }

# #---- testing ----
# #Single run
# test <- mask_impute_pool(CESD_data_wide, exposures = exposures,
#                          mechanism = "MCAR", method = "JMVN",
#                          mask_percent = "10%", num_impute = 2, save = "no")
# #Multiple runs
# test_2 <- replicate(2, mask_impute_pool(CESD_data_wide, exposures = exposures,
#                                                mechanism = "MCAR",
#                                                method = "PMM",
#                                                mask_percent = "10%",
#                                                num_impute = 2, save = "no"),
#                     simplify = FALSE)
# 
# #Formatting data
# formatted <- do.call(rbind, test_2)
# 
# #Summarizing results
# results <- formatted %>% group_by(Exposure) %>%
#   summarize_at(.vars = c("beta", "LCI", "UCI"), .funs = mean)
# 
# results2 <- formatted %>% group_by(Exposure) %>%
#   summarize_at(.vars = "beta", ~ quantile(.x, 0.025))


