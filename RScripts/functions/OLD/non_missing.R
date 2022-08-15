#This function finds the first or last non-missing element of a vector
#Use to find the first or last CYSC measure for a participant

non_missing <- function(vector, first = TRUE) {
  idx <- !is.na(vector) 
  
  if(first == TRUE){
    return(ifelse(any(idx), vector[min(which(idx))], NA))
  } else{
    return(ifelse(any(idx), vector[max(which(idx))], NA))
  }
}

# #Test vectors
# vector <- c(NA, 1, 2, NA, 3, NA)
# vector <- rep(NA, 5)
