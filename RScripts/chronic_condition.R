chronic_condition <- function(condition, 
                              condition_vars, condition_rx_vars, dataset){
  #---- get subset ----
  if(sum(is.na(condition_rx_vars)) != 0){
    subset <- dataset %>% 
      dplyr::select(all_of(condition_vars))
  } else{
    subset <- dataset %>% 
      dplyr::select(c(all_of(condition_vars), all_of(condition_rx_vars)))
    
  }
  
  #---- cleaning ----
  #Recoding self-report condition 
  for(i in 1:nrow(subset)){
    for(j in 1:length(condition_vars)){
      if(is.na(subset[i, j]) | subset[i, j] %in% c(0, 1, 2)){
        next
      } else{
        #Dispute and don't know
        if(subset[i, j] == 5){
          subset[i, 1:j] <- NA
          #Dispute and does not have condition
        } else if(subset[i, j] == 4){
          subset[i, 1:j] <- 0
          #Dispute and has condition
        } else if(subset[i, j] == 3){
          subset[i, 1:j] <- 1
        } 
      }
    }
  }
  
  #Recoding rx
  if(sum(is.na(condition_rx_vars)) == 0){
    condition_rx <- subset %>% dplyr::select(all_of(condition_rx_vars))
    #5 = "No"; 8 = "Don't know"; 9 = "Refused"
    condition_rx[condition_rx == 5] <- 0
    condition_rx[condition_rx == 8 | condition_rx == 9] <- NA
    
    #Update subset
    subset[, colnames(condition_rx)] <- condition_rx
  }
  
  #---- code condition ----
  subset %<>% 
    mutate("any" = rowSums(subset, na.rm = TRUE)) %>%
    mutate("ever_condition" = ifelse(any > 0, 1, 0))
  
  #update hrs_samp
  dataset[, head(colnames(subset), -2)] <- subset[, head(colnames(subset), -2)]
  dataset[, paste0("ever_", condition)] <- subset[, "ever_condition"]
  
  return(dataset)
}

# #Testing
# test <- chronic_condition("diabetes", paste0("r", seq(1, 13), "diab"),
#                           c(paste0("diabetes_rx_insulin", seq(4, 9)),
#                             paste0("diabetes_rx_swallowed", seq(4, 9))),
#                           hrs_samp)
# 
# test <- chronic_condition("hibp", paste0("r", seq(1, 13), "hibp"),
#                           paste0("bp_rx", seq(4, 9)),
#                           hrs_samp)
# 
# test_no_rx <- chronic_condition("mem", paste0("r", seq(5, 9), "memry"),
#                                 NA,
#                                 hrs_samp)
# 
# #Sanity check
# for(var in c(paste0("r", seq(1, 13), "hibp"), paste0("bp_rx", seq(4, 9)))){
#   print(var)
#   print(table(test[, var], useNA = "ifany"))
# }
# 
# table(test$ever_hibp, useNA = "ifany")
# 
# for(var in paste0("r", seq(5, 9), "memry")){
#   print(var)
#   print(table(test_no_rx[, var], useNA = "ifany"))
# }
# 
# table(test_no_rx$ever_mem, useNA = "ifany")

chronic_condition_per_wave <- 
  function(condition, condition_vars, condition_rx_vars, dataset){
    
    if(sum(is.na(condition_rx_vars)) != 0){
      #---- get subset ----
      subset <- dataset %>% 
        dplyr::select(all_of(condition_vars))
      #---- Single imputation using previous complete observation ----
      for (i in length(condition_vars)){
        hrs_samp %<>% 
          mutate("r9condition_impute" = 
                   ifelse(is.na(r4mstat), 
                          hrs_samp %>% 
                            dplyr::select(paste0("r", seq(4, 8), "mstat")) %>% 
                            apply(., 1, function(x) x[max(which(!is.na(x)))]), 
                          r4mstat)) %>% 
      }
    } else{
      #---- get subset ----
      subset <- dataset %>% 
        dplyr::select(c(all_of(condition_vars), all_of(condition_rx_vars)))
      #---- cleaning rx----
      #Recoding rx
      if(sum(is.na(condition_rx_vars)) == 0){
        condition_rx <- subset %>% dplyr::select(all_of(condition_rx_vars))
        #5 = "No"; 8 = "Don't know"; 9 = "Refused"
        condition_rx[condition_rx == 5] <- 0
        condition_rx[condition_rx == 8 | condition_rx == 9] <- NA
        
        # Update subset
        if (length(condition_rx_vars) > length(years)){
          # If there are one or more therapy for the same disease
          rx_new <-
            matrix(nrow = nrow(condition_rx), ncol = length(years))
          colnames(rx_new) <- c(paste0("r", years, condition, "_med"))
          for(i in 1:nrow(condition_rx)){
            for(j in 1:length(years)){
              rx_new[i, j] <- 
                ifelse(is.na(condition_rx[i, j]) & 
                         is.na(condition_rx[i, j + length(years)]), NA, 
                       ifelse((condition_rx[i, j] > 0 | 
                                 condition_rx[i, j + length(years)] > 0), 1, 0))
            }
          }
          subset %<>% dplyr::select(-condition_rx_vars)
          subset[, colnames(rx_new)] <- rx_new
        }
        else{
          subset[, colnames(condition_rx)] <- condition_rx
        }
      }
      
      #---- code condition ----
      condi_new <- matrix(nrow = nrow(dataset), ncol = length(condition_vars))
      colnames(condi_new) <- c(paste0(condition_vars, "_derived"))
      
      for(i in 1:nrow(subset)){
        for(j in 1:length(condition_vars)){
          condi_new[i, j]  <- 
            ifelse(is.na(subset[i,j]) & 
                     is.na(subset[i, j + length(condition_vars)]), NA, 
                   ifelse((subset[i,j] > 0 | 
                             subset[i,j + length(condition_vars)] > 0), 
                          1, 0))
        }
      }
      subset %<>% dplyr::select(-condition_vars)
      
    }
    #update hrs_samp
    dataset[, colnames(subset)] <- subset
    
    return(dataset)
  }

# Test
    
test <- chronic_condition_per_wave("diabetes", paste0("r", seq(4, 9), "diabe"),
                          c(paste0("diabetes_rx_insulin", seq(4, 9)),
                            paste0("diabetes_rx_swallowed", seq(4, 9))),
                          hrs_samp)


test <- chronic_condition_per_wave("hibpe", paste0("r", seq(4, 9), "hibpe"),
                          paste0("bp_rx", seq(4, 9)),
                          hrs_samp)

test_no_rx <- chronic_condition_per_wave("mem", paste0("r", seq(5, 9), "memry"),
                                NA,
                                hrs_samp)

#Sanity check
for(var in c(paste0("r", seq(1, 13), "hibp"), paste0("bp_rx", seq(4, 9)))){
  print(var)
  print(table(test[, var], useNA = "ifany"))
}

table(test$ever_hibp, useNA = "ifany")

for(var in paste0("r", seq(5, 9), "memry")){
  print(var)
  print(table(test_no_rx[, var], useNA = "ifany"))
}

table(test_no_rx$ever_mem, useNA = "ifany")

