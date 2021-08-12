fast_impute <- 
  function(predictor_matrix, data_wide, method, mechanism, mask_percent, m, 
           maxit, save = "no"){
    
    #---- where matrix ----
    where <- is.na(data_wide)*1  
    impute_vars <- rownames(predictor_matrix)
    
    #---- imputation 0: mean imputation ----
    avgs <- colMeans(data_wide[, -1], na.rm = TRUE)
    
    for(var in impute_vars){
      data_wide[where[, var] == 1 , var] <- avgs[var]
    }
    
    #---- **clean: drinking cat ----
    data_wide[, impute_vars[grep("drinking", impute_vars)]] <- 
      round(data_wide[, impute_vars[grep("drinking", impute_vars)]])
    
    #---- **draw: binary vars ----
    data_wide[, 
              impute_vars[-grep(pattern = "cesd|BMI|drinking", impute_vars)]] <- 
      apply(data_wide[, impute_vars[-grep(pattern = "cesd|BMI|drinking", 
                                          impute_vars)]] , 2, 
            function(x) rbinom(n = length(x), size = 1, prob = x))
    
    #---- **clean: marriage cats ----
    for(wave in seq(4, 9)){
      vars <- c("married_partnered", "not_married_partnered", "widowed")
      cols <- apply(expand.grid("r", wave, vars), 1, paste, collapse = "")
      subset <- data_wide[, cols]
      
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
      data_wide[, colnames(subset[, 1:3])] <- subset[, 1:3]
    }
    
    #---- pre-allocate chain storage ---- 
    if(save == "yes"){
      #colnames are imputation number:iteration number:stat number
      #stat number: 1 = mean; 2 = sd
      trace_data <- matrix(nrow = length(impute_vars), ncol = 2*m*maxit) %>% 
        set_colnames(apply(expand.grid(seq(1:m), seq(1:maxit), seq(1, 2)), 
                           1, paste, collapse = ":")) %>% 
        set_rownames(impute_vars)
    }
    
    #---- pre-allocate list of imputed datasets ----
    impute_list <- list()
    
    #---- imputation loop ----
    for(run in 1:m){
      imputed_data <- data_wide
      for(iter in 1:maxit){
        for(var in impute_vars){
          imputed_data[where[, var] == 1, var] <- NA
          #---- **PMM ----
          if(method == "PMM"){
            #fastMice will only do PMM is the outcome variable is factored
            imputed_data[, var] <- as.factor(imputed_data[, var])
            
            imputed_data[, var] <- 
              as.numeric(as.character(
                fill_NA_N(imputed_data, model = "pmm", posit_y = var, k = 10,
                          posit_x = c("intercept", names(
                            predictor_matrix[
                              var, which(predictor_matrix[var, ] == 1)])))))
          } else{
            imputed_data[, var] <- 
              as.numeric(fill_NA(imputed_data, model = "lm_bayes", 
                                 posit_y = var, 
                                 posit_x = c("intercept", names(
                                   predictor_matrix[
                                     var, 
                                     which(predictor_matrix[var, ] == 1)]))))
          }
          
          if(exists("trace_data")){
            #in the third slot: 1 = mean, 2 = sd
            trace_data[var, paste0(c(run, iter, 1), collapse = ":")] <- 
              mean(imputed_data[where[, var] == 1, var])
            trace_data[var, paste0(c(run, iter, 2), collapse = ":")] <- 
              sd(imputed_data[where[, var] == 1, var])
          }
        }
      }
      impute_list[[run]] <- imputed_data
    }
    
    #---- save results ----
    if(save == "yes"){
      #where matrix
      write_csv(as.data.frame(where), 
                file = here::here("MI datasets", 
                                  paste0("where_", tolower(method), "_", 
                                         tolower(mechanism), 
                                         as.numeric(sub("%","", 
                                                        mask_percent)), 
                                         ".csv")))
      
      #trace_data plots data
      write_csv(as.data.frame(trace_data), 
                file = here::here("MI datasets", 
                                  paste0("trace_data_", tolower(method), "_", 
                                         tolower(mechanism), 
                                         as.numeric(sub("%","", 
                                                        mask_percent)),
                                         ".csv")))
    }
    
    #---- return ----
    return(impute_list)
  }

# #---- function testing ----
# test <- fast_impute(predictor_matrix = predict, data_wide, method = "JMVN",
#                     mechanism = "MNAR", mask_percent = "10%", m = 2, maxit = 5,
#                     save = "yes")
















