#---- effect estimates ----
#Based on imputations
model_list <- 
  lapply(model_list <- vector(mode = "list", length(methods)),
         function(x) x <- lapply(x <- vector(mode = "list", length(mask_props)), 
                                 function(x) x <- 
                                   lapply(x <- vector(mode = "list", 
                                                      length = 4), 
                                          function(x) x <- vector(mode = "list", 
                                                                  length = num_impute))))

#naming layers of list
names(model_list) <- methods
for(i in 1:length(model_list)){
  names(model_list[[i]]) <- 100*mask_props
  for(j in 1:length(model_list[[i]])){
    names(model_list[[i]][[j]]) <- exposures
  }
}

for(m in methods[1:2]){
  for(i in 1:length(mask_props)){
    for(j in 1:num_impute){
      imputations <- get(paste0(tolower(m), "_mcar", mask_props[i]*100))
      imputed_data = complete(imputations, action = j)
      if(m == "JMVN"){
        #transform back to original values
        imputed_data[, paste0("r", seq(4, 9), "cesd")] <- 
          round(exp(imputed_data[, paste0("logr", seq(4, 9), "cesd")]) - 1)
      }
      
      if(m == "FCS"){
        imputed_data %<>% mutate_at(paste0("r", seq(4, 9), "cesd"), as.numeric)
      }
      
      #---- **E1a; E1b; E3 ----
      imputed_data %<>% 
        mutate("r4cesd_elevated" = ifelse(r4cesd > 4, 1, 0), 
               "r9cesd_elevated" = ifelse(r9cesd > 4, 1, 0), 
               "avg_cesd" = imputed_data %>% 
                 dplyr::select(paste0("r", seq(4, 9, by = 1), "cesd")) %>% 
                 rowMeans(), 
               "avg_cesd_elevated" = ifelse(avg_cesd > 4, 1, 0))
      
      #---- **E2 ----
      elevated_cesd <- imputed_data %>% 
        dplyr::select(paste0("r", seq(4, 9, by = 1), "cesd"))
      
      elevated_cesd <- (elevated_cesd > 4)*1
      
      imputed_data %<>% mutate("total_elevated_cesd" = rowSums(elevated_cesd))
      
      model_list[[m]][[i]][["CES-D Wave 4"]][[j]] <- 
        coxph(Surv(survtime, observed) ~ r4age_y_int + female + hispanic + 
                black + other + ed_cat + r4mstat_cat + ever_mem + 
                ever_arthritis + ever_stroke + ever_heart + ever_lung + 
                ever_cancer + ever_hibp + ever_diabetes + r4BMI + 
                drinking4_cat_impute + smoker + r4cesd_elevated, 
              data = imputed_data)
      
      model_list[[m]][[i]][["CES-D Wave 9"]][[j]] <- 
        coxph(Surv(survtime, observed) ~ r9age_y_int + female + hispanic + 
                black + other + ed_cat + r9mstat_cat + ever_mem + 
                ever_arthritis + ever_stroke + ever_heart + ever_lung + 
                ever_cancer + ever_hibp + ever_diabetes + r9BMI + 
                drinking9_cat_impute + smoker + r9cesd_elevated, 
              data = imputed_data)
      
      model_list[[m]][[i]][["Elevated CES-D Count"]][[j]] <- 
        coxph(Surv(survtime, observed) ~ r4age_y_int + female + hispanic + 
                black + other + ed_cat + r4mstat_cat + ever_mem + 
                ever_arthritis + ever_stroke + ever_heart + ever_lung + 
                ever_cancer + ever_hibp + ever_diabetes + r4BMI + 
                drinking4_cat_impute + smoker + total_elevated_cesd, 
              data = imputed_data)
      
      model_list[[m]][[i]][["Elevated Average CES-D"]][[j]] <- 
        coxph(Surv(survtime, observed) ~ r4age_y_int + female + hispanic + 
                black + other + ed_cat + r4mstat_cat + ever_mem + 
                ever_arthritis + ever_stroke + ever_heart + ever_lung + 
                ever_cancer + ever_hibp + ever_diabetes + r4BMI + 
                drinking4_cat_impute + smoker + avg_cesd_elevated, 
              data = imputed_data)
    }
  }
}

pooled_model_list <- 
  lapply(pooled_model_list <- vector(mode = "list", length(methods)),
         function(x) x <- lapply(x <- vector(mode = "list", length(mask_props)), 
                                 function(x) x <- 
                                   vector(mode = "list", length = 4)))

#naming layers of list
names(pooled_model_list) <- methods
for(i in 1:length(pooled_model_list)){
  names(pooled_model_list[[i]]) <- 100*mask_props
  for(j in 1:length(pooled_model_list[[i]])){
    names(pooled_model_list[[i]][[j]]) <- exposures
  }
}

for(m in methods[1:2]){
  for(i in as.character(mask_props*100)){
    for(j in exposures){
      pooled_model_list[[m]][[i]][[j]] <- 
        summary(pool(model_list[[m]][[i]][[j]]))
    }
  } 
}

for(m in methods[1:2]){
  for(prop in as.character(100*mask_props)){
    for(exposure in exposures){
      last_term <- nrow(pooled_model_list[[m]][[prop]][[exposure]])
      beta <- 
        exp(pooled_model_list[[m]][[prop]][[exposure]]$estimate[last_term])
      UCI <- 
        exp(pooled_model_list[[m]][[prop]][[exposure]]$estimate[last_term] + 
              pooled_model_list[[m]][[prop]][[exposure]]$std.error[last_term])
      LCI <- 
        exp(pooled_model_list[[m]][[prop]][[exposure]]$estimate[last_term] - 
              pooled_model_list[[m]][[prop]][[exposure]]$std.error[last_term])
      
      table_effect_ests[which(table_effect_ests$Exposure == exposure & 
                                table_effect_ests$Method == m & 
                                table_effect_ests$Missingness == 
                                paste0(prop, "%")), 
                        c("beta", "LCI", "UCI")] <- c(beta, LCI, UCI)
    }
  }
}
