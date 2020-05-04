fu_time <- function(vector){
  waves <- vector[1:5]
  ages <- vector[6:10]
  idx <- !is.na(waves) 
  
  if(any(idx)){
    first_slot <- min(which(idx))
    last_slot <- max(which(idx))
    
    return(ages[last_slot] - ages[first_slot])
  } else{
    return(NA)
  }
}

# #Test df
# test <- analytic_df %>% 
#   dplyr::select(contains(c("CYSC_ADJ", "AGE")))
# vector <- test[1, ]
