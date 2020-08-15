#This function creates a derived variable based on variables that can be 
#measured or self-reported within the RAND HRS dataset. The function picks up 
#all measured variables and fills in missing values with self-reported variables
#when available. This function also creates a set of columns [measure]_measured
#which indicates that the value is measured (vs. self-reported)

#Get measured BMI
BMI <- cbind(hrs_samp[, paste0("r", number_waves, "pmbmi")], 
             1 - is.na(hrs_samp[, paste0("r", number_waves, "pmbmi")]))
#Set NAs to 0 so we can pick up remaining NAs from self-report
BMI[is.na(BMI)] <- 0 
colnames(BMI) <- c(paste0(letter_waves, "BMI"), 
                   paste0(letter_waves, "BMI_measured"))

#Pick up self-reported BMI for empty (= 0) cells
pick_up <- (BMI[, paste0(letter_waves, "BMI")] == 0)*1*
  hrs_samp[, paste0("r", number_waves, "bmi")]

#Fill in self-reported BMI
BMI[, paste0(letter_waves, "BMI")] <- 
  BMI[, paste0(letter_waves, "BMI")] + pick_up

hrs_samp %<>% cbind(BMI)
