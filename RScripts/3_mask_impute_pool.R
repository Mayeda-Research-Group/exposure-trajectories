#This is the main analysis function, formatted for runs on the Hoffman cluster
#
#Input: CESD_data_wide.csv
#Output: filename depends on arguments
# format for filename [mechanism]_[method]_[percent]_([sens]).csv

# #for testing
# directory <- "/Users/CrystalShaw/Dropbox/Projects/"

mask_impute_pool <- 
  function(mechanism, method, mask_percent, directory, seed, save = "no", 
           sens = "no"){
    #---- load packages ----
    if (!require("pacman")){
      install.packages("pacman", repos='http://cran.us.r-project.org')
    }
    
    #make sure that devtools::install_github(repo = "amices/mice") was setup
    #on the cluster prior to running this script for LMM
    
    p_load("tidyverse", "magrittr", "broom", "ResourceSelection", "survival", 
           "openxlsx", "lubridate", "future.apply", "lme4", "devtools", 
           "miceFast", "mice")
    
    #No scientific notation
    options(scipen = 999)
    
    #---- source scripts ----
    source(paste0(directory, "exposure_trajectories/RScripts/mask.R"))
    source(paste0(directory, "exposure_trajectories/RScripts/fast_impute.R"))
    
    #---- set seed ----
    set.seed(seed)
    
    #---- read in data ----
    data_wide <- 
      read_csv(paste0(directory, 
                      "exposure_trajectories/data/CESD_data_wide.csv")) %>% 
      as.data.frame()
    
    #---- **exposure based on analysis ----
    if(sens == "yes"){
      exposure_cols <- c("r4cesd_elevated", "r9cesd_elevated", 
                         "avg_cesd_elevated", "prop_elevated_cesd")
      
      data_wide %<>% dplyr::select(-all_of(exposure_cols)) %>% 
        rename_at(vars(paste0(exposure_cols, "_sens")), ~ exposure_cols)
    } 
    
    beta_0_table <- 
      read_csv(paste0(directory, "exposure_trajectories/data/beta_0_table.csv"))
    beta_mat <- 
      read_csv(paste0(directory, "exposure_trajectories/data/beta_mat.csv")) %>% 
      as.data.frame() %>% set_rownames("beta")
    if(sens == "yes"){
      truth <- 
        read_csv(paste0(directory, 
                        "exposure_trajectories/data/truth_sens.csv")) %>% 
        dplyr::mutate("Type" = mechanism)
    } else{
      truth <- 
        read_csv(paste0(directory, "exposure_trajectories/data/truth.csv")) %>% 
        dplyr::mutate("Type" = mechanism)
    }
    
    #---- exposures ----
    exposures <- c("CES-D Wave 4", "CES-D Wave 9", "Elevated Average CES-D", 
                   "Elevated CES-D Prop")
    
    #---- create shell for output ----
    pooled_effect_ests <- 
      data.frame("Exposure" = exposures, "beta" = NA, "SE" = NA, "LCI" = NA, 
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
    
    # #Sanity check-- table of num missings per masking proportion
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
      
    } else if(!method %in% c("LMM", "CC", "Exposed", "Unexposed")){
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
      predict[, c("HHIDPN", "intercept", 
                  paste0("r", seq(3, 9), "conde_impute"), "white", "r3cesd", 
                  paste0("r", seq(3, 9), "shlt"), "age_death_y", 
                  "r4cesd_elevated", "r9cesd_elevated", "total_elevated_cesd",
                  "prop_elevated_cesd", "avg_cesd", "avg_cesd_elevated", 
                  "observed", paste0("r", seq(4, 9), "cesd_death2018"), 
                  paste0("r", seq(3, 8), "cesd_conde_impute"))] <- 0
      
      # #Sanity check
      # colSums(predict)
    }
    
    #---- **start time ----
    start <- Sys.time()
    
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
                                  #m = 2, maxit = 5,
                                  m = as.numeric(sub("%","", mask_percent)),
                                  maxit = max_it[method, mask_percent],
                                  save = save, directory = directory)
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
      impute_method <- mice::make.method(data_wide)
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
                      "age_death_y", "r4cesd_elevated", "r9cesd_elevated", 
                      "total_elevated_cesd", "avg_cesd", 
                      "avg_cesd_elevated")] <- ""
      
      impute_method <- impute_method[-which(impute_method == "")]
      #start <- Sys.time()
      
      data_imputed <- 
        mice::mice(data = data_wide, 
                   #m = 2, maxit = 2,
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
      data_imputed <- fast_impute(predictor_matrix = predict, data_wide, method,
                                  mechanism, mask_percent,
                                  #m = 2, maxit = 5,
                                  m = as.numeric(sub("%","", mask_percent)),
                                  maxit = max_it[method, mask_percent],
                                  save = save, directory = directory)
      #stop <- Sys.time() - start
      
      #from the mice package
      # #look at convergence
      # #10% missing needs maxit = 20
      # #20% missing needs maxit = 20
      # #30% missing needs maxit = 25
      # plot(data_imputed)
      
    } else if(method == "LMM"){
      #---- ****LMM ----
      #start <- Sys.time()
      data_imputed <- mice::mice(data = data_long, 
                                 #m = 2, maxit = 2,
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
    } 
    
    #---- **save results ----
    if(save == "yes" & sens == "no"){
      saveRDS(data_imputed, 
              file = paste0(directory, "exposure_trajectories/results/", 
                            "MI_datasets/", tolower(method), "_", 
                            tolower(mechanism), 
                            as.numeric(sub("%","", mask_percent))))
    }
    
    #---- fitted models ----
    #---- **CC ----
    if(method %in% c("CC", "Exposed", "Unexposed")){
      #---- ****post-process: exposures ----
      complete_data <- data_wide 
      
      if(sens == "yes"){
        complete_data %<>% 
          mutate("r4cesd_elevated" = ifelse(r4cesd >= 1, 1, 0), 
                 "r9cesd_elevated" = ifelse(r9cesd >= 1, 1, 0))
      } else{
        complete_data %<>% 
          mutate("r4cesd_elevated" = ifelse(r4cesd >= 4, 1, 0), 
                 "r9cesd_elevated" = ifelse(r9cesd >= 4, 1, 0))
      }
      
      #indicate where to take averages
      indicator <- complete_data %>% 
        dplyr::select(paste0("r", seq(4, 9), "cesd")) %>% is.na() %>% 
        rowSums() > 2
      complete_data %<>% mutate("avg_indicator" = 1 - (1*indicator))
      
      complete_data[, "avg_cesd"] <- 
        rowMeans(complete_data %>% 
                   dplyr::select(paste0("r", seq(4, 9), "cesd")), na.rm = TRUE) 
      complete_data[which(complete_data$avg_indicator == 0), "avg_cesd"] <- NA
      
      if(sens == "yes"){
        complete_data %<>% 
          mutate("avg_cesd_elevated" = ifelse(avg_cesd >= 1, 1, 0))
      } else{
        complete_data %<>% 
          mutate("avg_cesd_elevated" = ifelse(avg_cesd >= 4, 1, 0))
      }
      
      if(sens == "yes"){
        complete_data[, "prop_elevated_cesd"] <- 
          rowMeans(complete_data %>% 
                     dplyr::select(paste0("r", seq(4, 9), "cesd")) %>% 
                     mutate_all(function(x) ifelse(x >= 1, 1, 0)), na.rm = TRUE) 
      } else{
        complete_data[, "prop_elevated_cesd"] <- 
          rowMeans(complete_data %>% 
                     dplyr::select(paste0("r", seq(4, 9), "cesd")) %>% 
                     mutate_all(function(x) ifelse(x >= 4, 1, 0)), na.rm = TRUE) 
      }
      
      complete_data[which(complete_data$avg_indicator == 0), 
                    "prop_elevated_cesd"] <- NA
      
      if(method == "Exposed"){
        complete_data[which(is.na(complete_data$r4cesd_elevated)), 
                      "r4cesd_elevated"] <- 1
        complete_data[which(is.na(complete_data$r9cesd_elevated)), 
                      "r9cesd_elevated"] <- 1
        complete_data[which(is.na(complete_data$avg_cesd_elevated)), 
                      "avg_cesd_elevated"] <- 1
        complete_data[which(is.na(complete_data$prop_elevated_cesd)), 
                      "prop_elevated_cesd"] <- 1
      }
      
      if(method == "Unexposed"){
        complete_data[which(is.na(complete_data$r4cesd_elevated)), 
                      "r4cesd_elevated"] <- 0
        complete_data[which(is.na(complete_data$r9cesd_elevated)), 
                      "r9cesd_elevated"] <- 0
        complete_data[which(is.na(complete_data$avg_cesd_elevated)), 
                      "avg_cesd_elevated"] <- 0
        complete_data[which(is.na(complete_data$prop_elevated_cesd)), 
                      "prop_elevated_cesd"] <- 0
      }
      
      #---- ****models ----
      model_list <- vector(mode = "list", length = length(exposures)) %>% 
        set_names(exposures)
      
      for(exposure in exposures){
        if(exposure == "CES-D Wave 4"){
          model_list[[exposure]] <- 
            with(complete_data, 
                 coxph(Surv(survtime, observed) ~ r4not_married_partnered + 
                         r4widowed + ed_cat + r4drinking_cat + r4memrye_impute + 
                         r4stroke_impute + r4hearte_impute + r4lunge_impute + 
                         r4cancre_impute + r4hibpe_impute + r4diabe_impute + 
                         smoker + r4BMI + hispanic + black + other + female + 
                         r4age_y_int + r4cesd_elevated))
        } else if(exposure == "CES-D Wave 9"){
          model_list[[exposure]] <- 
            with(complete_data, 
                 coxph(Surv(survtime, observed) ~ r9not_married_partnered + 
                         r9widowed + ed_cat + r9drinking_cat + r9memrye_impute + 
                         r9stroke_impute + r9hearte_impute + r9lunge_impute + 
                         r9cancre_impute + r9hibpe_impute + r9diabe_impute + 
                         smoker + r9BMI + hispanic + black + other + female + 
                         r9age_y_int + r9cesd_elevated))
        } else if(exposure == "Elevated CES-D Prop"){
          model_list[[exposure]] <- 
            with(complete_data, 
                 coxph(Surv(survtime, observed) ~ r4not_married_partnered + 
                         r4widowed + ed_cat + r4drinking_cat + r4memrye_impute + 
                         r4stroke_impute + r4hearte_impute + r4lunge_impute + 
                         r4cancre_impute + r4hibpe_impute + r4diabe_impute + 
                         smoker + r4BMI + hispanic + black + other + female + 
                         r4age_y_int + prop_elevated_cesd))
        } else{
          model_list[[exposure]] <- 
            with(complete_data, 
                 coxph(Surv(survtime, observed) ~ r4not_married_partnered + 
                         r4widowed + ed_cat + r4drinking_cat + r4memrye_impute + 
                         r4stroke_impute + r4hearte_impute + r4lunge_impute + 
                         r4cancre_impute + r4hibpe_impute + r4diabe_impute + 
                         smoker + r4BMI + hispanic + black + other + female + 
                         r4age_y_int + avg_cesd_elevated))
        }
      }
      
    } else{
      model_list <- vector(mode = "list", length = length(exposures)) %>% 
        set_names(exposures)
      
      for(i in 1:(as.numeric(sub("%","", mask_percent)))){
        #for(i in 1:2){
        #---- **complete data ----
        if(method %in% c("JMVN", "PMM")){
          complete_data <- data_imputed[[i]]
        } else{
          complete_data <- complete(data_imputed, action = i)
        }
        
        if(method == "LMM"){
          #---- ****LMM: long --> wide ----
          complete_data %<>% 
            pivot_wider(id_cols = c("HHIDPN", all_of(time_invariant_vars)), 
                        names_from = wave, 
                        values_from = c("age_y_int", all_of(time_updated_vars)), 
                        names_glue = "{wave}{.value}")
        }
        
        if(method %in% c("JMVN", "LMM")){
          #---- ****post process: binary vars ----
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
                     prob = unlist(subset[, col]))
          }
        }
        
        #---- ****post process: dummy vars ----
        for(wave in seq(4, 9)){
          vars <- c("married_partnered", "not_married_partnered", "widowed")
          cols <- apply(expand.grid("r", wave, vars), 1, paste, collapse = "")
          subset <- complete_data[, cols]
          
          if(method %in% c("JMVN", "LMM")){
            #fix impossible probs
            subset[subset < 0] <- 0
            subset[subset > 1] <- 1
            
            for(column in 1:3){
              subset[, column] <- 
                rbinom(n = nrow(subset), size = 1, 
                       prob = unlist(subset[, column]))
            }
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
        
        #---- ****post-process: exposures ----
        if(sens == "yes"){
          complete_data %<>% 
            mutate("r4cesd_elevated" = ifelse(r4cesd >= 1, 1, 0), 
                   "r9cesd_elevated" = ifelse(r9cesd >= 1, 1, 0), 
                   "prop_elevated_cesd" = 
                     complete_data %>% 
                     dplyr::select(paste0("r", seq(4, 9), "cesd")) %>% 
                     mutate_all(function(x) ifelse(x >= 1, 1, 0)) %>% 
                     rowMeans(), 
                   "avg_cesd" = complete_data %>% 
                     dplyr::select(paste0("r", seq(4, 9), "cesd")) %>% 
                     rowMeans(), 
                   "avg_cesd_elevated" = ifelse(avg_cesd >= 1, 1, 0))
        } else{
          complete_data %<>% 
            mutate("r4cesd_elevated" = ifelse(r4cesd >= 4, 1, 0), 
                   "r9cesd_elevated" = ifelse(r9cesd >= 4, 1, 0), 
                   "prop_elevated_cesd" = 
                     complete_data %>% 
                     dplyr::select(paste0("r", seq(4, 9), "cesd")) %>% 
                     mutate_all(function(x) ifelse(x >= 4, 1, 0)) %>% 
                     rowMeans(), 
                   "avg_cesd" = complete_data %>% 
                     dplyr::select(paste0("r", seq(4, 9), "cesd")) %>% 
                     rowMeans(), 
                   "avg_cesd_elevated" = ifelse(avg_cesd >= 4, 1, 0))
        }
        
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
          } else if(exposure == "Elevated CES-D Prop"){
            model_list[[exposure]][[i]] <- 
              with(complete_data, 
                   coxph(Surv(survtime, observed) ~ r4not_married_partnered + 
                           r4widowed + ed_cat + r4drinking_cat + r4memrye_impute + 
                           r4stroke_impute + r4hearte_impute + r4lunge_impute + 
                           r4cancre_impute + r4hibpe_impute + r4diabe_impute + 
                           smoker + r4BMI + hispanic + black + other + female + 
                           r4age_y_int + prop_elevated_cesd))
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
    }
    
    #---- pooling models ----
    for(exposure in exposures){
      if(method %in% c("CC", "Exposed", "Unexposed")){
        pooled_model <- broom::tidy(model_list[[exposure]])
      } else{
        pooled_model <- summary(mice::pool(model_list[[exposure]]))
      }
      
      pooled_effect_ests[which(pooled_effect_ests$Exposure == exposure), 
                         c("beta", "SE")] <- 
        pooled_model[nrow(pooled_model), c("estimate", "std.error")]
    }
    
    pooled_effect_ests[, "LCI"] <- 
      pooled_effect_ests$beta - 1.96*pooled_effect_ests$SE
    
    pooled_effect_ests[, "UCI"] <- 
      pooled_effect_ests$beta + 1.96*pooled_effect_ests$SE
    
    #---- truth coverage ----
    pooled_effect_ests$capture_truth <- 
      (truth$beta > pooled_effect_ests$LCI)*
      (truth$beta < pooled_effect_ests$UCI)
    
    #---- total time in minutes ----
    pooled_effect_ests %<>% 
      mutate("time" = as.numeric(difftime(Sys.time(), start, units = "mins")))
    
    #---- seed used ----
    pooled_effect_ests %<>% 
      mutate("seed" = seed)
    
    #---- return values ----
    if(sens == "yes"){
      write_csv(pooled_effect_ests, 
                paste0(directory, "exposure_trajectories/results/", mechanism, 
                       "_", method, "_", mask_percent, "_sens.csv"), 
                append = TRUE)
    } else{
      write_csv(pooled_effect_ests, 
                paste0(directory, "exposure_trajectories/results/", mechanism, 
                       "_", method, "_", mask_percent, ".csv"), append = TRUE)
    }
  }

## Read in the arguments listed in the:
## R CMD BATCH --no-save --no-restore '--args mechanism="MNAR" method="JMVN" mask_percent="10%"' test_args2.R
## expression:
args=(commandArgs(TRUE))

## args is now a list of character vectors

## Check to see if arguments are passed and set default values if not,
## then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  mechanism = 'MNAR'
  method = 'JMVN'
  mask_percent = '10%'
  seed = 123
  save = 'no'
  sens = 'no'
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
## Now print values just to make sure:
print(mechanism)
print(method)
print(mask_percent)

mask_impute_pool(mechanism = mechanism, method = method, 
                 mask_percent = mask_percent, directory = "/u/home/c/cshaw343/", 
                 seed = seed, save = save, sens = sens)
