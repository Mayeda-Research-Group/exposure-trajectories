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
    mask_wave_specific <- c("mstat_impute", "drinking_impute", "memrye_impute", 
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
                      paste, collapse = ""))
    predict <- matrix(1, length(blocks), ncol = ncol(data_wide)) %>% 
      set_rownames(blocks) %>% 
      set_colnames(colnames(data_wide))
    
    #Don't use these as predictors
    predict[, c("HHIDPN", paste0("r", seq(4, 9), "mstat_cat"), 
                paste0("r", seq(4, 9), "drinking_cat"),
                paste0("r", seq(3, 9), "conde_impute"), "white", "age_death_y", 
                "observed", "CESD_missing", "r3cesd", "r4cesd_elevated", 
                "r9cesd_elevated",  "avg_cesd", "avg_cesd_elevated", 
                "total_elevated_cesd")] <- 0
    
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
    
    #---- **run imputation ----
    #---- ****JMVN ----
    if(method == "JMVN"){
      #Joint multivariate normal
      data_imputed <- mice(data = data_wide, m = num_impute, maxit = 20, 
                           method = "norm", predictorMatrix = predict, 
                           where = is.na(data_wide), 
                           blocks = as.list(rownames(predict)), seed = 20210126)
      
      #look at convergence
        #10% missing needs maxit = 20
        #20% missing needs maxit =  
        #30% missing needs maxit = 
      plot(data_imputed)
      
    } else if(method == "FCS"){
      #---- ****FCS ----
      #Fully conditional specification
      data_wide %<>% 
        mutate_all(as.factor) %>% 
        mutate_at(vars(c(paste0("r", seq(4, 9), "BMI"), "avg_cesd")), 
                  as.numeric)
      data_imputed <- mice(data = data_wide, m = num_impute, 
                           nnet.MaxNWts = 5000,
                           defaultMethod = 
                             c("norm", "logreg", "polyreg", "polr"),
                           predictorMatrix = predict, where = is.na(data_wide), 
                           blocks = as.list(rownames(predict)), 
                           seed = 20210126)
    } else if(method == "FCS Long"){
      #---- ****FCS Long ----
      #level-1 outcomes
      Y <- data_wide %>% 
        dplyr::select(c("HHIDPN", 
                        rownames(predict)[!rownames(predict) %in% 
                                            c("r4cesd_elevated", 
                                              "r9cesd_elevated", 
                                              "total_elevated_cesd", "avg_cesd", 
                                              "avg_cesd_elevated")]))
      Y %<>% cbind(matrix(rep(Y$HHIDPN, 6), ncol = 6, byrow = FALSE) %>% 
                     set_colnames(paste0("r", seq(4, 9), "HHIDPN")), .) %>%
        dplyr::select(-c("HHIDPN")) %>%
        mutate_all(as.factor) %>%
        pivot_longer(everything(), names_to = ".value", 
                     names_prefix = "r\\d") %>%
        mutate_at(c("cesd", "BMI", "shlt"), as.numeric) 
      
      #level-2 outcomes
      Y2 <- data_wide %>% 
        dplyr::select("HHIDPN", c("r4cesd_elevated", 
                                  "r9cesd_elevated", 
                                  "total_elevated_cesd", "avg_cesd", 
                                  "avg_cesd_elevated")) %>% 
        left_join(Y, by = "HHIDPN") %>% 
        mutate_at(c("r4cesd_elevated", "r9cesd_elevated", "avg_cesd_elevated"), 
                  as.factor) %>% 
        dplyr::select(c("r4cesd_elevated", 
                        "r9cesd_elevated", 
                        "total_elevated_cesd", "avg_cesd", 
                        "avg_cesd_elevated"))
      
      #level-1 predictors (non-missing)
      X <- data_wide %>% dplyr::select(paste0("r", seq(4, 9), "age_y_int")) %>% 
        pivot_longer(everything()) %>% dplyr::select("value") %>% cbind(1, .)
      
      #level-2 predictors (non-missing)
      X2 <- data_wide %>% 
        dplyr::select("HHIDPN", "ed_cat", "hispanic", "black", "other", 
                      "female", "death2018", "survtime") %>% 
        left_join(Y, by = "HHIDPN") %>% 
        dplyr::select("ed_cat", "hispanic", "black", "other", "female", 
                      "death2018", "survtime") %>% 
        mutate_all(as.factor) %>% mutate_at("survtime", as.numeric) %>% 
        cbind(1, .)
      
      #random effects (intercept only)
      Z <- matrix(1, ncol = 1, nrow = nrow(X))
      
      #get rid of HHIDPN from Y
      Y %<>% dplyr::select(-"HHIDPN")
      
      # data_imputed <- 
      #   jomo2(Y = Y, Y2 = Y2, X = X, X2 = X2, Z = Z, clus = data_wide$HHIDPN, 
      #         nburn = 10, nbetween = 10, nimp = num_impute)
      
      start <- Sys.time()
      data_imputed <- 
        jomo2com(Y.con = Y[, c("cesd", "BMI", "shlt")], 
                 Y.cat = Y[, c("mstat_impute", "drinking_impute", 
                               "memrye_impute", "stroke_impute", 
                               "hearte_impute", "lunge_impute", "cancre_impute", 
                               "hibpe_impute", "diabe_impute")],
                 Y.numcat = c(rep(3, 2), rep(2, 7)),
                 Y2.con = Y2[, c("total_elevated_cesd", "avg_cesd")], 
                 Y2.cat = Y2[, c("r4cesd_elevated", "r9cesd_elevated", 
                                 "avg_cesd_elevated")], 
                 Y2.numcat = c(rep(2, 3)), X = X, X2 = X2, Z = Z, 
                 clus = data_wide$HHIDPN, 
                 nburn = 10, nbetween = 10, nimp = num_impute)
      stop <- Sys.time() - start
      
    } else if(method == "PMM"){
      #---- ****PMM ----
      #Predictive Mean Matching
      #start <- Sys.time()
      data_imputed <- mice(data = data_wide, m = num_impute, maxit = 20, 
                           method = "pmm", predictorMatrix = predict, 
                           where = is.na(data_wide), 
                           blocks = as.list(rownames(predict)), 
                           seed = 20210126)
      #stop <- Sys.time() - start
      
      #look at convergence
      #10% missing needs maxit = 20
      #20% missing needs maxit =  
      #30% missing needs maxit = 
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
      
      #---- **post-processing ----
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
                 coxph(Surv(survtime, observed) ~ r4mstat_impute + ed_cat + 
                         r4drinking_impute + r4memrye_impute + r4stroke_impute + 
                         r4hearte_impute + r4lunge_impute + r4cancre_impute + 
                         r4hibpe_impute + r4diabe_impute + smoker + r4BMI + 
                         hispanic + black + other + female + r4age_y_int + 
                         r4shlt + r4cesd_elevated))
        } else if(exposure == "CES-D Wave 9"){
          model_list[[exposure]][[i]] <- 
            with(complete_data, 
                 coxph(Surv(survtime, observed) ~ r9mstat_impute + ed_cat + 
                         r9drinking_impute + r9memrye_impute + r9stroke_impute + 
                         r9hearte_impute + r9lunge_impute + r9cancre_impute + 
                         r9hibpe_impute + r9diabe_impute + smoker + r9BMI + 
                         hispanic + black + other + female + r9age_y_int + 
                         r9shlt + r9cesd_elevated))
        } else if(exposure == "Elevated CES-D Count"){
          model_list[[exposure]][[i]] <- 
            with(complete_data, 
                 coxph(Surv(survtime, observed) ~ r4mstat_impute + ed_cat + 
                         r4drinking_impute + r4memrye_impute + r4stroke_impute + 
                         r4hearte_impute + r4lunge_impute + r4cancre_impute + 
                         r4hibpe_impute + r4diabe_impute + smoker + r4BMI + 
                         hispanic + black + other + female + r4age_y_int + 
                         r4shlt + total_elevated_cesd))
        } else{
          model_list[[exposure]][[i]] <- 
            with(complete_data, 
                 coxph(Surv(survtime, observed) ~ r4mstat_impute + ed_cat + 
                         r4drinking_impute + r4memrye_impute + r4stroke_impute + 
                         r4hearte_impute + r4lunge_impute + r4cancre_impute + 
                         r4hibpe_impute + r4diabe_impute + smoker + r4BMI + 
                         hispanic + black + other + female + r4age_y_int + 
                         r4shlt + avg_cesd_elevated))
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


