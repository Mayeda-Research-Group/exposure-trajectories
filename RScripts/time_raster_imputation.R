#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "splines", "mice")

set.seed(20200520)

#---- Read in and format data ----
analytic_df <- read_csv(here::here("Data", "analytic_df.csv"))

#wide --> long
long_df <- analytic_df %>% 
  dplyr::select("HHIDPN", "female", "hispanic", "black", "other", 
                contains("AGE", ignore.case = FALSE), contains("CYSC_ADJ")) %>% 
  pivot_longer(cols = c(contains("AGE"), contains("CYSC_ADJ")), 
               names_to = c("wave", ".value"), 
               names_pattern = "(.)(.)") %>% 
  set_colnames(c("HHIDPN", "female", "hispanic", "black", "other", "wave", 
                 "Age", "CYSC")) %>% mutate_at("HHIDPN", as.numeric)

long_df[, "log_CYSC"] <- log(long_df$CYSC)

#Only consider those who are "age eligible" (50 or older)
long_df %<>% filter(Age >= 50)

#Get rid of all of the missing observations 
#We have 15,404 people with at least one measure
complete_data <- long_df %>% na.omit()
#length(unique(complete_data$HHIDPN))

#---- Take 10% of the data ----
test <- sample_frac(complete_data, size = 0.10)

#---- Identify first observation for individuals ----
id_first <- test %>% distinct(HHIDPN, .keep_all = TRUE) %>%
  mutate("first" = TRUE)

test %<>% left_join(., id_first, by = NULL) %>% arrange(HHIDPN) %>% 
  #Create a column for the imputation-- will induce missingness here
  mutate("pred_log_CYSC" = log_CYSC)

#---- Induce 25% missingness ----
#Specify that this only happens when the observation isn't "first"
not_first <- which(is.na(test$first))
missing <- sample(not_first, size = 0.25*length(not_first))

test[missing, "pred_log_CYSC"] <- NA

#---- Creating break points ----
youngest <- min(test$Age)
oldest <- max(test$Age)

#Specify break ages-- the assumption is linearity trend between break points
brk <- seq(45, oldest, by = 1)
k <- length(brk)

#Calculate B-spline
#Add a little noise before 45 because we eventually want a prediction there
X <- bs(test$Age, knots = brk, B = c(brk[1] - 0.0001, brk[k]),
        degree = 1)
X <- X[, -(k + 1)]
dimnames(X)[[2]] <- paste("x", 1:ncol(X), sep = "")
test <- cbind(test, X)

#Time warping
warp.setup <- data.frame(Age = brk,
                         Age2 = seq(45, oldest, length.out = k))
warp.model <- lm(Age2 ~ bs(Age, knots = brk[c(-1, -k)], degree = 1) - 1,
                 data = warp.setup, x = TRUE, y = TRUE)
warped.knots <- warp.model$y
maxage <- max(warped.knots)
Age2  <- predict(warp.model, newdata = test)
test <- cbind(test, Age2 = Age2)

#---- Imputation model ----
Y <- "pred_log_CYSC"
meth <- make.method(test)
meth[1:length(meth)] <- ""
meth[Y] <- "2l.pan" #Requesting 2-level imputation (for repeated measure)

pred <- make.predictorMatrix(test)
pred[1:nrow(pred), 1:ncol(pred)] <- 0

#Set HHIDPN as the class variable
pred[Y, "HHIDPN"] <- (-2)

#Set the fixed effects
pred[Y, c("female", "hispanic", "black", "other", "Age")] <- 1

#Set B-spline bases as random effect
pred[Y, paste("x", 1:60, sep = "")] <- 2

test_impute <- mice(test, meth = meth, pred = pred, m = 2, 
                    maxit = 10, seed = 52711, print = FALSE)
