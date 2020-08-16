#This function creates a derived variable based on variables that can be 
#measured or self-reported within the RAND HRS dataset. The function picks up 
#all measured variables and fills in missing values with self-reported variables
#when available. This function also creates a set of columns [measure]_measured
#which indicates that the value is measured (vs. self-reported)

measured_self_report <- function(data, measured_cols, self_cols, 
                                 derived_variable){
  #Get measured values
  derived <- cbind(data[, measured_cols], 
                   1 - is.na(data[, measured_cols]))
  #Set NAs to 0 so we can pick up remaining NAs from self-report
  derived[is.na(derived)] <- 0 
  colnames(derived) <- c(paste0(letter_waves, derived_variable), 
                         paste0(letter_waves, derived_variable, "_measured"))
  
  #Pick up self-reported measure for empty (= 0) cells
  pick_up <- (derived[, paste0(letter_waves, derived_variable)] == 0)*1*
    data[, self_cols]
  
  #Fill in self-reported data
  derived[, paste0(letter_waves, derived_variable)] <- 
    derived[, paste0(letter_waves, derived_variable)] + pick_up
  
  #---- return values ----
  return(derived)
}

# #---- unit testing ----
# test <- measured_self_report(hrs_samp, paste0("r", number_waves, "pmbmi"), 
#                              paste0("r", number_waves, "bmi"), "BMI")
# 
# colnames(test) == colnames(BMI)
# sum(is.na(test)) == sum(is.na(BMI))
# nrow(test)*ncol(test) == sum(is.na(test)) + 
#   sum((test == BMI) == TRUE, na.rm = TRUE)



