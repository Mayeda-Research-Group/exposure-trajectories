#This function creates a derived variable based on variables that can be either
# measured or self-reported within the RAND HRS dataset. The function picks up 
# all measured variables and fills in missing values with self-reported variables
# when available. This function also creates a set of columns [measure]_measured
# which indicates that the value is measured (vs. self-reported)

measured_self_report <- function(data, measured_cols, self_cols, 
                                 derived_variable, measured_waves_start, 
                                 all_waves_end){
  #Get measured values
  derived <- cbind(data[, measured_cols], 
                   1 - is.na(data[, measured_cols]))
  #Set NAs to 0 so we can pick up remaining NAs from self-report
  #NAs in self-reported data will be coded at -1 at first b/c I had a problem
  #with NA-ing out existing values otherwise
  derived[is.na(derived)] <- 0 
  colnames(derived) <- 
    c(paste0(seq(measured_waves_start, all_waves_end, by = 1), 
             derived_variable), 
      paste0(seq(measured_waves_start, all_waves_end, by = 1), derived_variable, 
             "_measured"))
  
  self_report <- data[, self_cols]
  self_report[is.na(self_report)] <- -1
  
  #Pick up self-reported measure for empty (= 0) cells
  pick_up <- (derived[, paste0(seq(measured_waves_start, all_waves_end, by = 1), 
                               derived_variable)] == 0)*1*
    self_report[, tail(self_cols, 
                       n = (all_waves_end - measured_waves_start + 1))]
  
  #Fill in self-reported data
  derived[, paste0(seq(measured_waves_start, all_waves_end, by = 1), 
                   derived_variable)] <- 
    derived[, paste0(seq(measured_waves_start, all_waves_end, by = 1), 
                     derived_variable)] + pick_up
  
  #Turn -1 to NA
  derived[derived < 0] <- NA
  
  #Attach self-report only waves
  derived[, paste0(seq(1, (measured_waves_start - 1), by = 1), 
                   derived_variable)] <- 
    data[, head(self_cols, n = (measured_waves_start - 1))]
  derived[, paste0(seq(1, (measured_waves_start - 1), by = 1), 
                   derived_variable, "_measured")] <- 0
  
  #---- return values ----
  return(derived[, c(paste0(seq(1, all_waves_end, by = 1), derived_variable), 
                     paste0(seq(1, all_waves_end, by = 1), derived_variable, 
                            "_measured"))])
}

# #---- unit testing ----
# test <- measured_self_report(hrs_samp, paste0("r", number_waves, "pmbmi"),
#                              paste0("r", number_waves, "bmi"), "BMI")
# 
# colnames(test) == colnames(BMI)
# sum(is.na(test)) == sum(is.na(BMI))
# nrow(test)*ncol(test) == sum(is.na(test)) +
#   sum((test == BMI) == TRUE, na.rm = TRUE)
# 
# data = hrs_samp
# measured_cols = paste0("r", seq(8, 13, by = 1), "pmwght")
# self_cols = paste0("r", number_waves, "weight")
# derived_variable = "weight"
# measured_waves_start = 8
# all_waves_end = 13




