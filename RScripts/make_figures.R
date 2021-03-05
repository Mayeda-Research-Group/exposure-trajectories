#---- package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "ghibli", "openxlsx", "magrittr")

#No scientific notation
options(scipen = 999)

#---- values ----
methods <- c("JMVN")
mask_props <- c(.10, .25, .40)

#---- note ----
# Since the difference between win and OS, put substituted directory here
# Yingyan's directory: C:/Users/yingyan_wu
#                      C:/Users/yingyan_wu/Dropbox
# Crystal's directory: /Users/CrystalShaw
#                     ~/Dropbox/Projects

#Changing directories here will change them throughout the script
path_to_box <- "/Users/CrystalShaw"
path_to_dropbox <- "~/Dropbox/Projects"

#---- read in data ----
CESD_data_wide <- 
  read_csv(paste0(path_to_dropbox, 
                  "/exposure_trajectories/data/", 
                  "CESD_data_wide.csv"), 
           col_types = cols(.default = col_double(), HHIDPN = col_character(), 
                            death2018 = col_integer(), 
                            ed_cat = col_factor(), 
                            r4mstat_cat = col_factor(), 
                            r9mstat_cat = col_factor(),
                            drinking4_cat_impute = col_factor(),
                            drinking9_cat_impute = col_factor(),
                            female = col_factor(), hispanic = col_factor(), 
                            black = col_factor(), other = col_factor(), 
                            smoker = col_integer()))

for(method in tolower(methods)){
  for(prop in 100*mask_props){
    assign(paste0(method, "_mcar", prop), 
           readRDS(here("MI datasets", paste0(method, "_mcar", prop))))
  }
}

table_effect_ests <- read.xlsx(paste0(path_to_dropbox, 
                                      "/exposure_trajectories/manuscript/", 
                                      "tables/main_text_tables.xlsx"))

#---- diagnostics: trace plots ----
#trace plots-- can plot these in ggplot if we want by accessing chainMean and 
# chainVar in imputation object. Right now not all the variables show in the 
# saved image
png(paste0("/Users/CrystalShaw/Dropbox/Projects/exposure_trajectories/",
           "manuscript/figures/mcar10_jmvn_traceplot.png"), 
    width = 7, height = 4.5, units = "in", res = 300)
plot(jmvn10)
dev.off()

#---- visualizations ----
#---- **effect estimates ----
table_effect_ests %<>% mutate_at(c("Missingness"), as.factor) %>% 
  mutate("Missingness Type" = "MCAR")
table_effect_ests$Method <- 
  factor(table_effect_ests$Method, 
         levels = c("Truth", "JMVN", "FCS", "JMVN Long", "FCS Long"))

ggplot(table_effect_ests, 
       aes(x = beta, y = Missingness, color = Method, shape = Method)) +
  geom_point(size = 3.5, position = position_dodge(0.60)) + 
  scale_shape_manual(values = c(rep("square", (nrow(table_effect_ests))))) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = .2, 
                position = position_dodge(0.60)) + theme_minimal() + 
  theme(legend.position = "bottom", legend.direction = "horizontal") + 
  scale_color_ghibli_d("LaputaMedium", direction = -1) + 
  scale_y_discrete(limits = rev(levels(table_effect_ests$Missingness))) + 
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") + 
  facet_grid(rows = vars(`Missingness Type`), cols = vars(Exposure)) 

ggsave(paste0(path_to_dropbox, "/exposure_trajectories/",
              "manuscript/figures/effect_ests.jpeg"), device = "jpeg", 
       dpi = 300, width = 9, height = 5, units = "in")
  
#---- **individual imputations ----
#---- ***fcs ----
fcs_mean_imputation_10 <- vector(mode = "list", length = 6)
for(i in 1:length(fcs_mean_imputation_10)){
  wave = i + 3
  fcs_mean_imputation_10[[i]] = 
    rowMeans((as.data.frame(fcs_mcar10$imp[[c(paste0("r", wave, "cesd"))]])) %>% 
               mutate_all(as.numeric))
}

#Do we want to stratify by wave?-- aggregate for now
num_missing = 0
for(i in 1:length(fcs_mean_imputation_10)){
  num_missing = num_missing + length(fcs_mean_imputation_10[[i]])
}

plot_data <- as.data.frame(matrix(nrow = num_missing, ncol = 2)) %>% 
  set_colnames(c("Observed", "Imputed"))
plot_data[, "Imputed"] <- unlist(fcs_mean_imputation_10)

for(i in 4:9){
  index <- as.numeric(rownames(fcs_mcar10[["imp"]][[paste0("r", i, "cesd")]]))
  start <- min(which(is.na(plot_data$Observed)))
  
  plot_data[start:(start + length(index) - 1), "Observed"] <- 
    CESD_data_wide[index, paste0("r", i, "cesd")]
}

#plot
ggplot(data = plot_data, aes(x = Observed, y = Imputed)) + 
  geom_point(color = "#B4DAE5FF") + 
  geom_smooth(method = lm, se = FALSE, color = "#F0D77BFF") +
  geom_abline(slope = 1, intercept = 0, color = "#5C5992FF", lty = "dashed", 
              size = 1) +
  ggtitle(paste0("Missingness Pattern: MCAR 10% \n", 
                 "Imputation Strategy: Fully Conditional Specification")) + 
  theme_minimal()  

ggsave(paste0("/Users/CrystalShaw/Dropbox/Projects/exposure_trajectories/",
              "manuscript/figures/mcar10_fcs_obs_pred.jpeg"), 
       device = "jpeg", width = 7, height = 4.5, units = "in", dpi = 300)

#---- **derived exposures ----





mean_imputation_mcar <- 
  lapply(mean_imputation_mcar <- vector(mode = "list", length(mask_props)),
         function(x) x <- lapply(x <- vector(mode = "list", 4),
                                 function(x) x <- vector(mode = "list", 
                                                         length = 6)))
#naming layers of list
names(mean_imputation_mcar) <- mask_props*100
for(i in 1:length(mean_imputation_mcar)){
  names(mean_imputation_mcar[[i]]) <- c("jmvn", "fcs", "jmvn_long", "fcs_long")
}

for(i in 1:length(mean_imputation_mcar)){
  data <- get(paste0("jmvn_mcar", mask_props[i]*100))
  for(j in 1:6){
    wave = j + 3
    mean_imputation_mcar[[i]][["jmvn"]][[j]] <- 
      rowMeans(data$imp[[c(paste0("logr", wave, "cesd"))]])
  }
}

plot_data <- as.data.frame(matrix(nrow = num_missing, ncol = 2)) %>% 
  set_colnames(c("Observed", "Imputed"))
plot_data[, "Imputed"] <- unlist(mean_imputation)

observed_data <- 
  cbind(CESD_data_wide %>% 
          dplyr::select(paste0("r", seq(4, 9), "cesd")) %>% 
          pivot_longer(everything(), names_to = "rwave", 
                       values_to = "complete"), 
        mcar10 %>% 
          dplyr::select(paste0("logr", seq(4, 9), "cesd")) %>% 
          pivot_longer(everything(), names_to = "logrwave", 
                       values_to = "masked")) %>% 
  arrange(rwave) %>% filter(is.na(masked))

plot_data[, "Observed"] <- observed_data$complete
plot_data %<>% mutate("log_observed" = log(1 + Observed), 
                      "trans_imputed" = round(exp(Imputed) - 1))

# #Sanity check
# head(CESD_data_wide$r4cesd[as.numeric(names(mean_imputation[[1]]))])
# head(plot_data)
# tail(CESD_data_wide$r9cesd[as.numeric(names(mean_imputation[[6]]))])
# tail(plot_data)

#plot
ggplot(data = plot_data, aes(x = log_observed, y = Imputed)) + 
  geom_point(color = "#B4DAE5FF") + 
  geom_smooth(method = lm, se = FALSE, color = "#F0D77BFF") +
  geom_abline(slope = 1, intercept = 0, color = "#5C5992FF", lty = "dashed", 
              size = 1) + labs(x = "Log(1 + Observed)") + 
  ggtitle(paste0("Missingness Pattern: MCAR 10% \n", 
                 "Imputation Strategy: Joint Multivariate Normal")) + 
  theme_minimal()  

ggsave(paste0("/Users/CrystalShaw/Dropbox/Projects/exposure_trajectories/",
              "manuscript/figures/mcar10_jmvn_obs_pred_log.jpeg"), 
       device = "jpeg", width = 7, height = 4.5, units = "in", dpi = 300)

ggplot(data = plot_data, aes(x = Observed, y = trans_imputed)) + 
  geom_point(color = "#B4DAE5FF") + 
  geom_smooth(method = lm, se = FALSE, color = "#F0D77BFF") +
  geom_abline(slope = 1, intercept = 0, color = "#5C5992FF", lty = "dashed", 
              size = 1) + labs(y = "Exp(Imputed) - 1") + 
  ggtitle(paste0("Missingness Pattern: MCAR 10% \n", 
                 "Imputation Strategy: Joint Multivariate Normal")) + 
  theme_minimal()  

ggsave(paste0("/Users/CrystalShaw/Dropbox/Projects/exposure_trajectories/",
              "manuscript/figures/mcar10_jmvn_obs_pred.jpeg"), 
       device = "jpeg", width = 7, height = 4.5, units = "in", dpi = 300)







#---- ***visualize imputations ----



# #Sanity check
# head(CESD_data_wide$r4cesd[as.numeric(names(mean_imputation[[1]]))])
# head(plot_data)
# tail(CESD_data_wide$r9cesd[as.numeric(names(mean_imputation[[6]]))])
# tail(plot_data)


