impute_ages <- function(vector){
  #No one in our analytic sample has a vector of ages completely missing 
  #because everyone has age at Wave O
  last_age <- length(vector)
  missing_age <- which(vector == 999)
  
  for(slot in missing_age){
    vector[slot] <- vector[last_age] - 2*(last_age - slot)
  }
  
  return(vector)
}

# #test vectors
# vector <- c(999, 75, 77, 79, 81)
# vector <- c(75, 999, 999, 81, 83)
