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
                 "Age", "CYSC"))

long_df[, "log_CYSC"] <- log(long_df$CYSC)

#Only consider those who are "age eligible" (50 or older)
long_df %<>% filter(Age >= 50)

#Get rid of all of the missing observations 
#We have 15,404 people with at least one measure
complete_data <- long_df %>% na.omit()
#length(unique(complete_data$HHIDPN))

#---- Take 10% of the data ----
test <- sample_frac(complete_data, size = 0.10)

#---- Imputation model ----
youngest <- min(MCAR_25_train$Age)
oldest <- max(MCAR_25_train$Age)

#Define break points
brk <- seq(45, oldest, by = 5)
k <- length(brk)

#Time warping
warp.setup <- data.frame(Age = brk,
                         Age2 = seq(45, oldest, length.out = k))
warp.model <- lm(Age2 ~ bs(Age, knots = brk[c(-1, -k)], degree = 1) - 1,
                 data = warp.setup, x = TRUE, y = TRUE)
warped.knots <- warp.model$y
maxage <- max(warped.knots)
Age2  <- predict(warp.model, newdata = MCAR_25_train)
MCAR_25_train <- cbind(MCAR_25_train, Age2 = Age2)

id <- unique(MCAR_25_train$HHIDPN)
MCAR_25_train2 <- appendbreak(MCAR_25_train, brk, id = id, 
                              warp.model = warp.model, typ = "sup")
table(data2$typ)

