#Want to detect whether someone was observed in the dataset at age 70 or older
#This will ensure that the person survived up to at least age 70

detect_70 <- function(vector){
  
  at_least_70 <- which(vector >= 70)
  
  keep = ifelse(sum(at_least_70 >= 1), 1, 0)
  
  return(keep)
}

# #test vectors
# vector <- c(999, 75, 77, 79, 81)
# vector <- c(60, 61, 62, 63, 64)
# vector <- c(68, 69, 70, 71, 72)
