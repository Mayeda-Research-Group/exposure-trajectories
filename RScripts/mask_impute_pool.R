mask_impute_pool <- 
  function(data_wide, exposures, mechanism, method, mask_percent, truth, 
           save = "no"){
    
    #---- create shell for output ----
    pooled_effect_ests <- 
      data.frame("Exposure" = exposures, "beta" = NA, "SD" = NA, "LCI" = NA, 
                 "UCI" = NA, "Method" = method, "Missingness" = mask_percent, 
                 "Type" = mechanism, "capture_truth" = NA)
    
    #---- create incomplete data ----
    data_wide <- mask(data_wide, mechanism, mask_percent)
    
    time_updated_vars <- c("married_partnered", "not_married_partnered", 
                           "widowed", "drinking_cat", "memrye_impute", 
                           "stroke_impute", "hearte_impute", "lunge_impute", 
                           "cancre_impute", "hibpe_impute", "diabe_impute", 
                           "cesd", "BMI")
    
    time_invariant_vars <- c("ed_cat", "white", "black", "hispanic", "other", 
                             "female", "survtime", "death2018", "smoker", 
                             "observed")
    
    #---- convert to long?? ----
    if(method == "LMM"){
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
    if(method == "LMM"){
      blocks <- time_updated_vars
      predict <- matrix(1, length(blocks), ncol = ncol(data_long)) %>% 
        set_rownames(blocks) %>% 
        set_colnames(colnames(data_long))
      
      #Don't use these as predictors
      predict[, c("wave", "observed", "white")] <- 0
      
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
                  paste0("r", seq(4, 9), "shlt"),"age_death_y", "white", 
                  "observed", "r3cesd", "r4cesd_elevated", "r9cesd_elevated", 
                  "avg_cesd", "avg_cesd_elevated", "total_elevated_cesd")] <- 0
      
      #---- ****time-updated var models ----
      for(var in time_updated_vars){
        #can't predict itself
        predict[paste0("r", seq(4, 9), var), paste0("r", seq(4, 9), var)] <- 
          (diag(x = 1, nrow = 6, ncol = 6) == 0)*1
        
        #can't predict in the same wave (all missing)
        predictors <- time_updated_vars[which(time_updated_vars != var)]
        for(predictor in predictors){
          predict[paste0("r", seq(4, 9), var), 
                  paste0("r", seq(4, 9), predictor)] <- 
            (diag(x = 1, nrow = 6, ncol = 6) == 0)*1
        }
        
        #use time-updated age
        predict[paste0("r", seq(4, 9), var), 
                paste0("r", seq(4, 9), "age_y_int")] <- 
          diag(x = 1, nrow = 6, ncol = 6)
      }
    }
    
    #---- **run imputation ----
    max_it <- tibble("Method" = c("FCS", "JMVN", "PMM", "LMM"), 
                     "10%" = c(20, 20, 20, 5),
                     "20%" = c(25, 20, 20, 5),
                     "30%" = c(25, 25, 25, 5)) %>% 
      column_to_rownames("Method")
    
    #---- ****JMVN ----
    if(method == "JMVN"){
      #Joint multivariate normal
      #start <- Sys.time()
      data_imputed <- mice(data = data_wide, 
                           #m = 1, maxit = 25,
                           m = as.numeric(sub("%","", mask_percent)), 
                           maxit = max_it[method, mask_percent], 
                           method = "norm", predictorMatrix = predict, 
                           where = is.na(data_wide), 
                           blocks = as.list(rownames(predict)), seed = 20210126)
      #stop <- Sys.time() - start
      
      # #look at convergence
      #   #10% missing needs maxit = 20
      #   #20% missing needs maxit = 20
      #   #30% missing needs maxit = 25
      # plot(data_imputed)
      
    } else if(method == "FCS"){
      #---- ****FCS ----
      #Fully conditional specification
      impute_method <- make.method(data_wide)
      impute_method[c(paste0("r", seq(4, 9), "married_partnered"),
                 paste0("r", seq(4, 9), "not_married_partnered"),
                 paste0("r", seq(4, 9), "widowed"),
                 paste0("r", seq(4, 9), "memrye_impute"),
                 paste0("r", seq(4, 9), "stroke_impute"),
                 paste0("r", seq(4, 9), "hearte_impute"),
                 paste0("r", seq(4, 9), "lunge_impute"),
                 paste0("r", seq(4, 9), "cancre_impute"),
                 paste0("r", seq(4, 9), "hibpe_impute"),
                 paste0("r", seq(4, 9), "diabe_impute"))] <- "logreg"
      
      impute_method[c(paste0("r", seq(4, 9), "drinking_cat"),
               paste0("r", seq(4, 9), "BMI"), 
               paste0("r", seq(4, 9), "cesd"))] <- "norm"
      
      impute_method[c(paste0("r", seq(4, 9), "shlt"), "r3cesd",
               "age_death_y", "r4cesd_elevated", "r9cesd_elevated", 
               "total_elevated_cesd", "avg_cesd", "avg_cesd_elevated")] <- ""
      
     impute_method <- impute_method[-which(impute_method == "")]
      
      # data_wide %<>%
      #   mutate_at(vars(c(paste0("r", seq(4, 9), "married_partnered"),
      #                    paste0("r", seq(4, 9), "not_married_partnered"),
      #                    paste0("r", seq(4, 9), "widowed"),
      #                    paste0("r", seq(4, 9), "memrye_impute"),
      #                    paste0("r", seq(4, 9), "stroke_impute"),
      #                    paste0("r", seq(4, 9), "hearte_impute"),
      #                    paste0("r", seq(4, 9), "lunge_impute"),
      #                    paste0("r", seq(4, 9), "cancre_impute"),
      #                    paste0("r", seq(4, 9), "hibpe_impute"),
      #                    paste0("r", seq(4, 9), "diabe_impute"), "smoker",
      #                    "hispanic", "black", "other", "female", "death2018")),
      #             as.factor)
      
      #start <- Sys.time()
      data_imputed <- mice(data = data_wide, 
                           m = as.numeric(sub("%","", mask_percent)), 
                           #m = 1, maxit = 25,
                           maxit = max_it[method, mask_percent],
                           method = impute_method,
                           predictorMatrix = predict, where = is.na(data_wide), 
                           blocks = as.list(rownames(predict)),
                           seed = 20210126)
      #end <- Sys.time() - start
      
      # #look at convergence
      #   #10% missing needs maxit = 20
      #   #20% missing needs maxit = 25
      #   #30% missing needs maxit = 25
      # plot(data_imputed)
      
      } else if(method == "PMM"){
        #---- ****PMM ----
        #Predictive Mean Matching
        impute_method <- make.method(data_wide)
        impute_method[c(paste0("r", seq(4, 9), "married_partnered"),
                        paste0("r", seq(4, 9), "not_married_partnered"),
                        paste0("r", seq(4, 9), "widowed"),
                        paste0("r", seq(4, 9), "memrye_impute"),
                        paste0("r", seq(4, 9), "stroke_impute"),
                        paste0("r", seq(4, 9), "hearte_impute"),
                        paste0("r", seq(4, 9), "lunge_impute"),
                        paste0("r", seq(4, 9), "cancre_impute"),
                        paste0("r", seq(4, 9), "hibpe_impute"),
                        paste0("r", seq(4, 9), "diabe_impute"))] <- "logreg"
        
        impute_method[c(paste0("r", seq(4, 9), "drinking_cat"),
                        paste0("r", seq(4, 9), "BMI"), 
                        paste0("r", seq(4, 9), "cesd"))] <- "norm"
        
        impute_method[c(paste0("r", seq(4, 9), "shlt"), "r3cesd",
                        "age_death_y", "r4cesd_elevated", "r9cesd_elevated", 
                        "total_elevated_cesd", "avg_cesd", "avg_cesd_elevated")] <- ""
        
        impute_method <- impute_method[-which(impute_method == "")]
        
        # data_wide %<>% 
        #   mutate_at(vars(c(paste0("r", seq(4, 9), "married_partnered"),
        #                    paste0("r", seq(4, 9), "not_married_partnered"),
        #                    paste0("r", seq(4, 9), "widowed"),
        #                    paste0("r", seq(4, 9), "memrye_impute"), 
        #                    paste0("r", seq(4, 9), "stroke_impute"),
        #                    paste0("r", seq(4, 9), "hearte_impute"),
        #                    paste0("r", seq(4, 9), "lunge_impute"), 
        #                    paste0("r", seq(4, 9), "cancre_impute"), 
        #                    paste0("r", seq(4, 9), "hibpe_impute"), 
        #                    paste0("r", seq(4, 9), "diabe_impute"), "smoker", 
        #                    "hispanic", "black", "other", "female", "death2018")), 
        #             as.factor)
        
        #start <- Sys.time()
        data_imputed <- mice(data = data_wide, 
                             #m = as.numeric(sub("%","", mask_percent)), 
                             #maxit = max_it[method, mask_percent], 
                             m = 2, maxit = 5,
                             method = "pmm", donors = 5, 
                             predictorMatrix = predict, 
                             where = is.na(data_wide), 
                             blocks = as.list(rownames(predict)), 
                             seed = 20210126)
        #stop <- Sys.time() - start
        
        # #look at convergence
        # #10% missing needs maxit = 20
        # #20% missing needs maxit = 20
        # #30% missing needs maxit = 25
        # plot(data_imputed)
        
      } else if(method == "LMM"){
        #---- ****LMM ----
        #impute with Bayesian longitudinal model (allowing for heteroskedasticity)
        start <- Sys.time()
        data_imputed <- mice(data = data_long, 
                             #m = as.numeric(sub("%","", mask_percent)), 
                             #maxit = max_it[method, mask_percent],
                             m = 2, maxit = 2,
                             method = "2l.lmer", predictorMatrix = predict, 
                             where = is.na(data_long), 
                             blocks = as.list(rownames(predict)), 
                             seed = 20210126)
        stop <- Sys.time() - start
        
        # } else if(method == "2l.fcs"){
        #   #---- ****2l.fcs ----
        #   #Longitudinal fully conditional specification
        #   #Fully conditional specification
        #   data_long %<>% 
        #     mutate_at(vars(c("married_partnered", "not_married_partnered", 
        #                      "widowed", "memrye_impute", "stroke_impute", 
        #                      "hearte_impute", "lunge_impute", "cancre_impute", 
        #                      "hibpe_impute", "diabe_impute", "black", "hispanic", 
        #                      "other", "female", "death2018", "smoker")), 
        #               as.factor)
        #   
        #   start <- Sys.time()
        #   data_imputed <- mice(data = data_long, m = num_impute, 
        #                        maxit = 10, 
        #                        defaultMethod = 
        #                          c("2l.norm", "2l.bin", "2l.norm", "2l.norm"), 
        #                        predictorMatrix = predict, 
        #                        where = is.na(data_long), 
        #                        blocks = as.list(rownames(predict)), 
        #                        seed = 20210126)
        #   stop <- Sys.time() - start
        #   
        #   # #look at convergence
        #   # #10% missing needs maxit = 
        #   # #20% missing needs maxit = 
        #   # #30% missing needs maxit = 
        #   plot(data_imputed) 
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
    
    for(i in 1:(as.numeric(sub("%","", mask_percent)))){
      complete_data <- complete(data_imputed, action = i)
      
      if(method == "JMVN"){
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
                  "lunge_impute", "cancre_impute", "hibpe_impute", 
                  "diabe_impute")
        cols <- apply(expand.grid("r", waves, vars), 1, paste, collapse = "")
        
        #fix impossible probs
        subset <- complete_data[, cols]
        subset[subset < 0] <- 0
        subset[subset > 1] <- 1
        
        for(col in cols){
          complete_data[, col] <- 
            rbinom(n = nrow(complete_data), size = 1, prob = subset[, col])
        }
      }
      
      if(method %in% c("FCS", "PMM")){
        #---- **post process: dummy vars ----
        for(wave in seq(4, 9)){
          vars <- c("married_partnered", "not_married_partnered", "widowed")
          cols <- apply(expand.grid("r", wave, vars), 1, paste, collapse = "")
          subset <- complete_data[, cols] %>% as.matrix() %>% 
            apply(2, as.numeric)
          probs <- 1/rowSums(subset)
          
          for(row in which(probs != 1)){
            cats <- as.numeric(which(subset[row, ] == 1))
            this_one <- sample(cats, size = 1, 
                               prob = rep(probs[row], floor(1/probs[row])))
            subset[row, ] <- 0
            subset[row, this_one] <- 1
          }
          complete_data[, colnames(subset)] <- subset
        }
      }
      
      if(method == "LMM"){
        #---- **long --> wide ----
        complete_data %<>% 
          pivot_wider(id_cols = c("HHIDPN", all_of(time_invariant_vars)), 
                      names_from = wave, 
                      values_from = c("age_y_int", all_of(time_updated_vars)), 
                      names_glue = "{wave}{.value}")
        
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
                  "lunge_impute", "cancre_impute", "hibpe_impute", 
                  "diabe_impute")
        cols <- apply(expand.grid("r", waves, vars), 1, paste, collapse = "")
        
        #fix impossible probs
        subset <- complete_data[, cols]
        subset[subset < 0] <- 0
        subset[subset > 1] <- 1
        
        for(col in cols){
          complete_data[, col] <- 
            rbinom(n = nrow(complete_data), size = 1, 
                   prob = as.matrix(subset[, col]))
        }
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
                         r4age_y_int + r4cesd_elevated))
        } else if(exposure == "CES-D Wave 9"){
          model_list[[exposure]][[i]] <- 
            with(complete_data, 
                 coxph(Surv(survtime, observed) ~ r9not_married_partnered + 
                         r9widowed + ed_cat + r9drinking_cat + r9memrye_impute + 
                         r9stroke_impute + r9hearte_impute + r9lunge_impute + 
                         r9cancre_impute + r9hibpe_impute + r9diabe_impute + 
                         smoker + r9BMI + hispanic + black + other + female + 
                         r9age_y_int + r9cesd_elevated))
        } else if(exposure == "Elevated CES-D Count"){
          model_list[[exposure]][[i]] <- 
            with(complete_data, 
                 coxph(Surv(survtime, observed) ~ r4not_married_partnered + 
                         r4widowed + ed_cat + r4drinking_cat + r4memrye_impute + 
                         r4stroke_impute + r4hearte_impute + r4lunge_impute + 
                         r4cancre_impute + r4hibpe_impute + r4diabe_impute + 
                         smoker + r4BMI + hispanic + black + other + female + 
                         r4age_y_int + total_elevated_cesd))
        } else{
          model_list[[exposure]][[i]] <- 
            with(complete_data, 
                 coxph(Surv(survtime, observed) ~ r4not_married_partnered + 
                         r4widowed + ed_cat + r4drinking_cat + r4memrye_impute + 
                         r4stroke_impute + r4hearte_impute + r4lunge_impute + 
                         r4cancre_impute + r4hibpe_impute + r4diabe_impute + 
                         smoker + r4BMI + hispanic + black + other + female + 
                         r4age_y_int + avg_cesd_elevated))
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
    
    #---- truth coverage ----
    pooled_effect_ests$capture_truth <- 
      (truth$beta > pooled_effect_ests$LCI)*
      (truth$beta < pooled_effect_ests$UCI)
    
    #---- return values ----
    return(pooled_effect_ests)
    }

# #---- testing ----
# #Single run
# test <- mask_impute_pool(CESD_data_wide, exposures = exposures,
#                          mechanism = "MCAR", method = "PMM", truth = truth,
#                          mask_percent = "10%", save = "yes")
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
# 
# #---- code optimization ----
# library("profvis")
# profvis::profvis(
#   mask_impute_pool(CESD_data_wide, exposures = exposures,
#                    mechanism = "MNAR", method = "PMM",
#                    mask_percent = "10%", truth = truth, save = "no"))
# 
# 
