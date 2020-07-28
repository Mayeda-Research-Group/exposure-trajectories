impute_ages <- function(vector){
  #For complete vectors, do nothing
  if(sum(is.na(vector)) == 0){
    return(vector)
    } else{
      
    last_age <- max(which(!is.na(vector)))
    missing_age <- which(is.na(vector))
    
    for(slot in missing_age){
      if(slot < last_age){
        vector[slot] <- vector[last_age] - 24*(last_age - slot)
      } else{
        vector[slot] <- vector[last_age] + 24*(slot - last_age)
      }
    }
    return(vector)
  }
}

# #test vectors
# vector <- c(NA, 75, 77, 79, 81)
# vector <- c(75, 77, 79, 81, NA)
# vector <- c(75, 77, 79, 81, 83)
