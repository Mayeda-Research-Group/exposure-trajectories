#Start with a simple model CysC ~ age + sex/gender + race/ethnicity + 
#                                   sex/gender*age + race/ethnicity*age

#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "lme4", "tidyverse", "magrittr")

#---- Read in data ----
analytic_df <- read_csv(here::here("Data", "analytic_df.csv"))

#---- wide --> long ----
long_df <- analytic_df %>% 
  dplyr::select("HHIDPN", "female", "hispanic", "black", "other", 
                contains("AGE", ignore.case = FALSE), contains("CYSC_ADJ")) %>% 
  pivot_longer(cols = c(contains("AGE"), contains("CYSC_ADJ")), 
               names_to = c("wave", ".value"), 
               names_pattern = "(.)(.)") %>% 
  set_colnames(c("HHIDPN", "female", "hispanic", "black", "other", "wave", 
                 "Age", "CYSC"))

#---- look at the outcome variable ----
hist(long_df$CYSC)
boxplot(long_df$CYSC)

#CYSC is really skewed-- take the log of the outcome
long_df[, "log_CYSC"] <- log(long_df$CYSC)

#Look at transformed outcome-- super symmetric!
hist(long_df$log_CYSC)
boxplot(long_df$log_CYSC)

#---- Model ----
#When we use the untransformed outcome, variance of residuals has fanning 
#behavior; violation of assumptions

#Model with random intercept for HHIDPN

lmm1 <- lmer(log_CYSC ~ female + hispanic + black + other + Age + 
               female*Age + hispanic*Age + black*Age + other*Age + 
               (1|HHIDPN), 
             data = long_df)

#Check model assumptions
plot(lmm1)
hist(resid(lmm1))
qqnorm(resid(lmm1))
qqline(resid(lmm1))

summary(long_df$log_CYSC)

#Number of people
#I think it's this: 15,653
#Try running model with random slope in people with 3
length(which(analytic_df %>% dplyr::select(contains("CYSC_ADJ")) %>% rowSums(., na.rm = TRUE) != 0))

