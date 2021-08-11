mask_impute_pool <- 
  function(data_wide, exposures, mechanism, method, mask_percent, beta_0_table,
           beta_mat, truth, save = "no"){
    
    #---- create shell for output ----
    pooled_effect_ests <- 
      data.frame("Exposure" = exposures, "beta" = NA, "SD" = NA, "LCI" = NA, 
                 "UCI" = NA, "Method" = method, "Missingness" = mask_percent, 
                 "Type" = mechanism, "capture_truth" = NA)
    
    #---- create incomplete data ----
    data_wide <- 
      mask(data_wide, mechanism, mask_percent, beta_0_table, beta_mat)
    
    time_updated_vars <- c("married_partnered", "not_married_partnered", 
                           "widowed", "drinking_cat", "memrye_impute", 
                           "stroke_impute", "hearte_impute", "lunge_impute", 
                           "cancre_impute", "hibpe_impute", "diabe_impute", 
                           "cesd", "BMI")
    
    time_invariant_vars <- c("ed_cat", "black", "hispanic", "other", 
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
      predict[, c("wave", "observed")] <- 0
      
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
      
      #Don't use these as predictors
      predict[, c("HHIDPN", paste0("r", seq(4, 9), "married_partnered"), 
                  paste0("r", seq(3, 9), "conde_impute"), "white", "r3cesd", 
                  paste0("r", seq(3, 9), "shlt"), "age_death_y", 
                  "r4cesd_elevated", "r9cesd_elevated", "total_elevated_cesd", 
                  "avg_cesd", "avg_cesd_elevated", "observed", 
                  paste0("r", seq(4, 9), "cesd_death2018"), 
                  paste0("r", seq(3, 8), "cesd_conde_impute"))] <- 0
      
      # #Sanity check
      # colSums(predict)
    }
    
    #---- **run imputation ----
    max_it <- tibble("Method" = c("FCS", "JMVN", "PMM", "LMM"), 
                     "10%" = c(10, 10, 10, 5),
                     "20%" = c(15, 10, 10, 5),
                     "30%" = c(15, 15, 15, 5)) %>% 
      column_to_rownames("Method")
    
    #---- ****JMVN ----
    if(method == "JMVN"){
      #Joint multivariate normal
      #start <- Sys.time()
      data_imputed <- fast_impute(predictor_matrix = predict, data_wide, method, 
                                  mechanism, mask_percent, 
                                  m = 2, maxit = 5,
                                  # m = as.numeric(sub("%","", mask_percent)), 
                                  # maxit = max_it[method, mask_percent], 
                                  save = save)
      
      #stop <- Sys.time() - start
      
      #this is from mice package
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
      
      impute_method[c(paste0("r", seq(3, 9), "shlt"), "r3cesd", 
                      "r3cesd_conde_impute",
                      #paste0("r", seq(4, 9), "married_partnered"),
                      "age_death_y", "r4cesd_elevated", "r9cesd_elevated", 
                      "total_elevated_cesd", "avg_cesd", 
                      "avg_cesd_elevated")] <- ""
      
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
                           #m = 2, maxit = 5,
                           m = as.numeric(sub("%","", mask_percent)),
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
      #start <- Sys.time()
      data_imputed <- mice(data = data_wide, 
                           m = as.numeric(sub("%","", mask_percent)),
                           maxit = max_it[method, mask_percent],
                           #m = 2, maxit = 5,
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
      #start <- Sys.time()
      data_imputed <- mice(data = data_long, 
                           # m = 2, maxit = 5,
                           m = as.numeric(sub("%","", mask_percent)),
                           maxit = max_it[method, mask_percent],
                           method = "2l.lmer", predictorMatrix = predict, 
                           where = is.na(data_long), 
                           blocks = as.list(rownames(predict)), 
                           seed = 20210126)
      
      #stop <- Sys.time() - start
      
      # #look at convergence
      # #10% missing needs maxit = 
      # #20% missing needs maxit = 
      # #30% missing needs maxit = 
      # plot(data_imputed)
      
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
      #for(i in 1:2){
      complete_data <- complete(data_imputed, action = i)
      
      if(method == "LMM"){
        #---- **LMM: long --> wide ----
        complete_data %<>% 
          pivot_wider(id_cols = c("HHIDPN", all_of(time_invariant_vars)), 
                      names_from = wave, 
                      values_from = c("age_y_int", all_of(time_updated_vars)), 
                      names_glue = "{wave}{.value}")
      }
      
      if(method %in% c("JMVN", "LMM")){
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
        subset[is.na(subset)] <- 0.5
        
        for(col in cols){
          complete_data[, col] <- 
            rbinom(n = nrow(complete_data), size = 1, 
                   prob = unlist(subset[, col]))
        }
      }
      
      #---- **post process: dummy vars ----
      for(wave in seq(4, 9)){
        vars <- c("married_partnered", "not_married_partnered", "widowed")
        cols <- apply(expand.grid("r", wave, vars), 1, paste, collapse = "")
        subset <- complete_data[, cols]
        
        #fix impossible probs
        subset[subset < 0] <- 0
        subset[subset > 1] <- 1
        
        for(column in 1:3){
          subset[, column] <- 
            rbinom(n = nrow(subset), size = 1, prob = unlist(subset[, column]))
        }
        
        subset[, "sum"] <- rowSums(subset, na.rm = TRUE)
        
        for(row in 1:nrow(subset)){
          if(subset[row, "sum"] == 1){
            next
          } else if(subset[row, "sum"] %in% c(0, 3)){
            this_cat <- sample(c(1, 2, 3), size = 1)
            subset[row, ] <- 0
            subset[row, this_cat] <- 1
          } else{
            which_cols <- which(subset[row, ] == 1)
            this_cat <- sample(which_cols, size = 1)
            subset[row, ] <- 0
            subset[row, this_cat] <- 1
          }
        }
        complete_data[, colnames(subset[, 1:3])] <- subset[, 1:3]
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
