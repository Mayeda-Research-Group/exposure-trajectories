# This function finds the most recent non-missing wave value for a variable
#   for ever/never chronic condition per wave

impute_chronic_condition <- 
  function(condition, condition_vars, waves, dataset) {
    # getting the subset
    subset <- dataset %>%
      dplyr::select(all_of(condition_vars))
    
    for(j in length(condition_vars):1){
      wave <- waves[j] 
      if (j != 1) {
        subset[, paste0("r", waves[j], condition, "_impute")] <- 
          ifelse(is.na(subset[, j]), 
                 subset %>% 
                   dplyr::select(paste0("r", waves[1:j-1], condition)) %>%
                   apply(., 1, function(x) x[max(which(!is.na(x)))]),
                 subset[, j])
      } else{
        subset[, paste0("r", waves[j], condition, "_impute")] <- 
          subset[, j]
      }
    }
    subset %<>% select(contains("impute"))
    dataset[, colnames(subset)] <- subset
    return(dataset)
  }


# Test
# test2 <- impute_chronic_condition("diabe", paste0("r", seq(1, 9), "diabe"),
#                                    seq(1, 9), hrs_samp)
# table(test2$r1diabe_impute, test2$r1diabe, useNA = "ifany")
# table(test2$r2diabe_impute, test2$r2diabe, useNA = "ifany")
# table(test2$r3diabe_impute, test2$r3diabe, useNA = "ifany")
# table(test2$r4diabe_impute, test2$r4diabe, useNA = "ifany")
# table(test2$r5diabe_impute, test2$r5diabe, useNA = "ifany")
# table(test2$r6diabe_impute, test2$r6diabe, useNA = "ifany")
# table(test2$r7diabe_impute, test2$r7diabe, useNA = "ifany")
# table(test2$r8diabe_impute, test2$r8diabe, useNA = "ifany")
# table(test2$r9diabe_impute, test2$r9diabe, useNA = "ifany")
# View(test2 %>% dplyr::select(contains(paste0("r", seq(1, 9), "diabe")),
#                              -contains("rx"))%>%
#        filter(!rowSums(is.na(.)) == 0))
# 
# test3 <- impute_chronic_condition("hibpe", paste0("r", seq(4, 9), "hibpe"),
#                                    seq(4, 9), hrs_samp)
# table(test3$r5hibpe_impute, test3$r5hibpe, useNA = "ifany")
# table(test3$r6hibpe_impute, test3$r6hibpe, useNA = "ifany")
# table(test3$r7hibpe_impute, test3$r7hibpe, useNA = "ifany")
# table(test3$r8hibpe_impute, test3$r8hibpe, useNA = "ifany")
# table(test3$r9hibpe_impute, test3$r9hibpe, useNA = "ifany")
# View(test3 %>% dplyr::select(contains(paste0("r", seq(4, 9), "hibpe")),
#                              -contains("rx"))%>%
#        filter(!rowSums(is.na(.)) == 0))

impute_status <- function(status, status_vars, waves, impute_waves, dataset) {
  # getting the subset
  subset <- dataset %>%
    dplyr::select(all_of(status_vars))
  
  for(j in 1:length(impute_waves)){
    imp_wave <- impute_waves[j] 
    subset[, paste0("r", imp_wave, status, "_impute")] <- 
        ifelse(is.na(subset[, paste0("r", imp_wave, status)]), 
               subset %>% 
                 dplyr::select(paste0("r", waves[1:imp_wave - 1], status)) %>%
                 apply(., 1, function(x) x[max(which(!is.na(x)))]),
               subset[, paste0("r", imp_wave, status)])
      subset[, paste0("r", imp_wave, status, "_impute")] <-
        ifelse(is.na(subset[, paste0("r", imp_wave, status, "_impute")]),
               subset %>% 
                 dplyr::select(
                   paste0("r", (imp_wave + 1):length(waves), status)) %>%
                 apply(., 1, function(x) x[min(which(!is.na(x)))]),
               subset[, paste0("r", imp_wave, status, "_impute")])
    }
  subset %<>% select(contains("impute"))
  dataset[, colnames(subset)] <- subset
  return(dataset)
}

# # Test
# test <- impute_status("mstat", paste0("r", seq(1, 13), "mstat"),
#                                    seq(1, 13), seq(4, 9),  hrs_samp)
# sum(test$r4mstat != test$r4mstat_impute, na.rm = TRUE)
# View(test %>%
#        dplyr::select("HHIDPN", paste0("r", number_waves, "mstat"),
#                      "r4mstat_impute"))
# View(test %>%
#        dplyr::select("HHIDPN", paste0("r", number_waves, "mstat"),
#                      "r4mstat_impute") %>% filter(is.na(r4mstat)))
# View(test %>%
#        dplyr::select("HHIDPN", paste0("r", number_waves, "mstat"),
#                      "r5mstat_impute") %>% filter(is.na(r5mstat)))

# #Test
# test2 <- impute_status("drinkd", paste0("r", seq(3, 13), "drinkd"),
#                        seq(3, 13), seq(4, 9),  hrs_samp)
# sum(test2$r4drinkd != test2$r4drinkd_impute, na.rm = TRUE)
# View(test2 %>%
#        dplyr::select("HHIDPN", paste0("r", seq(3,13), "drinkd"),
#                      "r4drinkd_impute"))
# View(test2 %>%
#        dplyr::select("HHIDPN", paste0("r", seq(3,13), "drinkd"),
#                      "r4drinkd_impute") %>% filter(is.na(r4drinkd)))
# View(test2 %>%
#        dplyr::select("HHIDPN", paste0("r", seq(3,13), "drinkd"),
#                      "r5drinkd_impute") %>% filter(is.na(r5drinkd)))
