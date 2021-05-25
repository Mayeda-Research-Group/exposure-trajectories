#---- optimized betas table ----
mechanisms <- c("MAR", "MNAR")
percents <- c(10, 20, 30)

beta_0_table <- expand_grid(mechanisms, percents) %>% 
  mutate("beta0" = 0)

for(mechanism in mechanisms){
  for(percent in percents){
    optimized <- 
      read_rds(file = paste0(path_to_dropbox, "/exposure_trajectories/data/", 
                             "optimized_masking_intercepts/optim_", mechanism, 
                             percent, ".RDS"))
    
    beta_0_table[which(beta_0_table$mechanisms == mechanism & 
                         beta_0_table$percents == percent), "beta0"] <- 
      optimized$minimum
  }
}

#---- beta matrix ----
beta_mat <- #effect sizes
  matrix(c(log(1.1), log(1.05), log(1.05), log(1.25), log(1.1), log(1.25)), 
         nrow = 1) %>% 
  #MAR
  set_colnames(c("cesdpre", "condepre", "cesdpre_condepre",
                 #MNAR
                 "death2018", "cesdcurrent", "death2018_cesdcurrent")) %>% 
  set_rownames(c("beta"))

#---- Expit function ----
expit <- function(x) {
  output <- (exp(x)/(1+exp(x)))
  return(output)
}

logit <- function(x){
  output <- log(x/(1-x))
  return(output)
}

#---- Mask function ----
mask <- function(data_wide, mechanism, mask_percent, beta_0_table, beta_mat){
  
  mask_prop <- as.numeric(sub("%","", mask_percent))/100
  
  #---- create incomplete data ----
  if(mechanism == "MCAR"){
    #---- MCAR ----
    #it's easier to do this with my own code than the ampute function in MICE, 
    # which requires specifying all possible missing patterns you'd like it to 
    # consider
    total_indices <- nrow(data_wide)*6 #6 waves of data per person
    mask_index <- sample(seq(1, total_indices), 
                         size = floor(mask_prop*total_indices), 
                         replace = FALSE)
  } else{
    subset <- data_wide
    
    if(mechanism == "MNAR"){
      #---- MNAR ----
      for(wave in seq(4, 9)){
        subset %<>% 
          mutate(!!paste0("r", wave, "pcesd") := 
                   expit(beta_mat["beta", "death2018"]*death2018 + 
                           beta_mat["beta", "cesdcurrent"]*
                           !!sym(paste0("r", wave, "cesd")) + 
                           beta_mat["beta", "death2018_cesdcurrent"]*
                           !!sym(paste0("r", wave, "cesd_death2018")) + 
                           as.numeric(beta_0_table[which(
                             beta_0_table$mechanisms == mechanism & 
                               beta_0_table$percents == mask_prop*100), 
                             "beta0"])))
      }
    } else if(mechanism == "MAR"){
      #---- MAR ----
      for(wave in seq(4, 9)){
        subset %<>% 
          mutate(!!paste0("r", wave, "pcesd") := 
                   expit(beta_mat["beta", "cesdpre"]*
                           !!sym(paste0("r", wave - 1, "cesd")) + 
                           beta_mat["beta", "condepre"]*
                           !!sym(paste0("r", wave - 1, "conde_impute")) + 
                           beta_mat["beta", "cesdpre_condepre"]*
                           !!sym(paste0("r", wave - 1, "cesd_conde_impute")) + 
                           as.numeric(beta_0_table[which(
                             beta_0_table$mechanisms == mechanism & 
                               beta_0_table$percents == mask_prop*100), 
                             "beta0"])))
      }
    }
    
    subset %<>% dplyr::select(contains("pcesd", ignore.case = FALSE))
    
    for (j in 1:ncol(subset)){
      subset[, paste0("r", j + 3, "cesd_missing")] <- 
        rbinom(nrow(subset), size = 1, prob = subset[[j]])
    }
    
    subset_long <- subset %>% select(contains("cesd_missing")) %>%
      pivot_longer(everything(), names_to = "orig_varname", 
                   values_to = "cesd_missing")
    
    mask_index <- which(subset_long$cesd_missing == 1)
  }
  
  #---- masking wave-specific values ----
  mask_wave_specific <- c("married_partnered", "not_married_partnered", 
                          "widowed", "drinking_cat", "memrye_impute", 
                          "stroke_impute", "hearte_impute", "lunge_impute", 
                          "cancre_impute", "hibpe_impute", "diabe_impute", 
                          "cesd", "BMI", "shlt")
  
  for(var in mask_wave_specific){
    #mask values
    data_long <- data_wide %>% 
      dplyr::select("HHIDPN", paste0("r", seq(4, 9), var)) %>% 
      pivot_longer(-"HHIDPN")
    
    data_long[mask_index, "value"] <- NA
    
    #need to get rid of rows with NA values to pivot back to wide
    data_long %<>% na.omit() %>%
      pivot_wider(names_from = "name", values_from = "value")
    
    #mask in wide dataset with masked values
    data_wide[which(data_wide$HHIDPN %in% data_long$HHIDPN),
              paste0("r", seq(4, 9), var)] <- data_long[, -1]
    data_wide[which(!data_wide$HHIDPN %in% data_long$HHIDPN),
              paste0("r", seq(4, 9), var)] <- NA
  }
  
  #masking derived values
  data_wide %<>% 
    mutate("r4cesd_elevated" = ifelse(r4cesd > 4, 1, 0), 
           "r9cesd_elevated" = ifelse(r9cesd > 4, 1, 0), 
           "total_elevated_cesd" = 
             rowSums(data_wide %>% 
                       dplyr::select(paste0("r", seq(4, 9), "cesd")) %>% 
                       mutate_all(function(x) x > 4)), 
           "avg_cesd" = 
             rowMeans(data_wide %>% 
                        dplyr::select(paste0("r", seq(4, 9), "cesd"))), 
           "avg_cesd_elevated" = ifelse(avg_cesd > 4, 1, 0))
  
  #---- check missings ----
  # make sure no one is missing every cesd measure
  data_wide %<>%
    mutate("CESD_missing" =
             rowSums(data_wide %>%
                       dplyr::select(paste0("r", seq(4, 9), "cesd")) %>%
                       mutate_all(function(x) is.na(x)))) %>%
    filter(CESD_missing < 6) %>% dplyr::select(-CESD_missing)
  
  # Return the dataset
  # return(mask_index)
  return(data_wide)
}

# # Test
# path_to_box <- "C:/Users/yingyan_wu"
# path_to_dropbox <- "C:/Users/yingyan_wu/Dropbox"
# 
# #---- read in analytical sample ----
# CESD_data_wide <-
#   read_csv(paste0(path_to_dropbox,
#                   "/exposure_trajectories/data/",
#                   "CESD_data_wide.csv"))
# # geting the coefficients in MAR and MNAR models
# coef_df <- data.frame(
#   intercept = exp(beta_0),
#   age_j = exp(beta_age),
#   CESD_j_1 = exp(beta_cesdpre),
#   conde_j_1 = exp(beta_condepre),
#   cesd_j = exp(beta_cesdcurrent)
# ) %>% print()

# MCAR
# test <- mask(CESD_data_wide, "MCAR", "20%")
# {index <- mask(CESD_data_wide, "MCAR", "20%")
#   length(index)/(9445*6)
# }
# test %>%
#   select(contains("cesd"),
#          -contains(c("elevated", "avg", "r3"))) %>%
#   pivot_longer(everything(),
#                names_to = "origvar",
#                values_to = "CESD") %>%
#     dplyr::count(CESD) %>%
#     mutate(prop = prop.table(n))
# 
# view(test[is.na(test$r9cesd), ])
# 
# test_func <- function(data_wide, mechanism, mask_percent){
#   mask_prop <- as.numeric(sub("%","", mask_percent))/100
#   n <- nrow(mask(data_wide, mechanism, mask_percent))
#   return(n)
# }
# 
# size = 1000
# system.time(
#   test_sample_size <- replicate(size, test_func(CESD_data_wide, "MCAR", "30%")))
# summary(test_sample_size)
# 
# MAR
# test1 <- mask(CESD_data_wide, "MAR", "20%")
# test1 %>%
#   select(contains("cesd"),
#          -contains(c("elevated", "avg", "r3"))) %>%
#   pivot_longer(everything(),
#                names_to = "origvar",
#                values_to = "CESD") %>%
#   dplyr::count(CESD) %>%
#   mutate(prop = prop.table(n))
# 
# test1 %>%
#   select(contains("r5")) %>%
#   filter(is.na(r5cesd))
# 
# size = 1000
# system.time(
#   test1_sample_size <- replicate(size, test_func(CESD_data_wide, "MAR", "30%")))
# summary(test1_sample_size)

# # MNAR
# test2 <- mask(CESD_data_wide, "MNAR", "30%")
# test2 %>%
#   select(contains("cesd"),
#          -contains(c("elevated", "avg", "r3"))) %>%
#   pivot_longer(everything(),
#                names_to = "origvar",
#                values_to = "CESD") %>%
#   dplyr::count(CESD) %>%
#   mutate(prop = prop.table(n))
# test2 %>%
#   select(contains("r5")) %>%
#   filter(is.na(r5cesd))
# 
# size = 1000
# system.time(
#   test2_sample_size <- replicate(size, test_func(CESD_data_wide, "MNAR", "30%")))
# summary(test2_sample_size)
