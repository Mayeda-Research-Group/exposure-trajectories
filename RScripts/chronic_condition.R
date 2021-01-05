chronic_condition <- function(condition, 
                              condition_vars, condition_rx_vars, dataset){
  #---- get subset ----
  if(is.na(condition_rx_vars)){
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
  if(!is.na(condition_rx_vars)){
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




