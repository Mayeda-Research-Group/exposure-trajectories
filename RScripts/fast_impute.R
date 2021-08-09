fast_impute <- function(predictor_matrix, data_wide, method, m, maxit){
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
  data_wide[, impute_vars[-grep(pattern = "cesd|BMI|drinking", impute_vars)]] <- 
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
  
  #---- imputation ----
  for(run in 1:m){
    for(iter in 1:maxit){
      for(var in impute_vars){
        data_wide[where[, var] == 1, var] <- NA
        test <- 
          fill_NA_N(data_wide, model = "pmm", posit_y = var, 
                    posit_x = names(
                      predictor_matrix[var, 
                                       which(predictor_matrix[var, ] == 1)]))
        
        test_mice <- mice(data = data_wide, 
                          m = 1, maxit = 5, method = "pmm", donors = 5, 
                          predictorMatrix = predict, 
                          where = is.na(data_wide), 
                          blocks = as.list(rownames(predict)), 
                          seed = 20210126)
      }
    }
  }
  
  
  
  
  
  
  
  #make list to store imputed matrices
  for(run in 1:m){
    for(var in rownames(predict)){
      wherever where = 1
    }
    list[m] <- imputed matrix
  }
  
  return(list(imputed matrices))
  
}