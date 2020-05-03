fu_time <- function(vector){
  idx <- !is.na(vector) 
  
  if(any(idx)){
    first_slot <- min(which(idx))
    last_slot <- max(which(idx))
    
    return(2*(last_slot - first_slot))
  } else{
    return(NA)
  }
}

# #Test vectors
# vector <- c(NA, 1, NA, 2, NA, NA)
# vector <- rep(NA, 5)
