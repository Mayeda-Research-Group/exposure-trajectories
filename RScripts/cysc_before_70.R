cysc_before_70 <- function(vector){
  waves <- vector[1:5]
  ages <- vector[6:10]
  idx <- !is.na(waves) 
  
  if(sum(idx) >= 1){
    count <- sum(ages[idx] < 70)
    return((count >= 1)*1)
  } else{
    return(0)
  }
}

#Test df
test <- hrs_samp %>%
  dplyr::select(contains(c("CYSC_ADJ", "age_y")))
vector <- test[2, ]
