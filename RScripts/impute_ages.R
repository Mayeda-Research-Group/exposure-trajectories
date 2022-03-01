#This function imputes a vector of ages in months by adding or subtracting 
# in 24-month increments (HRS has visits every two years)

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
# # converted ages in years to months
# vector <- c(NA, 75*12, 77*12, 79*12, 81*12)
# impute_ages(vector)
# vector <- c(75*12, 77*12, 79*12, 81*12, NA)
# impute_ages(vector)
# vector <- c(75*12, NA, NA, 81*12, 83*12)
# impute_ages(vector)
# vector <- c(75*12, 77*12, 79*12, 81*12, 83*12)
# impute_ages(vector)
