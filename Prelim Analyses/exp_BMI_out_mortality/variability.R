#---- Package + dataset loading  ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "readr", "tidyverse", "magrittr", "plyr", "haven", "labelled", 
       "lubridate")

#Changing directories here will change them throughout the script
# Yingyan's directory: C:/Users/yingyan_wu
#                      C:/Users/yingyan_wu/Dropbox
# Crystal's directory: /Users/CrystalShaw
#                     ~/Dropbox/Projects
path_to_box <- "C:/Users/yingyan_wu"
path_to_dropbox <- "C:/Users/yingyan_wu/Dropbox"
hrs_samp <- read_csv(paste0(path_to_dropbox,
                            "/exposure_trajectories/data/",
                            "hrs_samp_6BMI_waves4-9.csv"))
#----BMI trajectory ----
bmi_data <- hrs_samp %>% 
  dplyr::select("HHIDPN", contains("age_y_int"), contains("BMI"), 
                                      -contains("BMI_measured"))
bmi_temp <- bmi_data%>%
  tidyr::gather(oldvar, age, `4age_y_int`:`9age_y_int`, factor_key=TRUE )%>%
    mutate("wave" = ifelse(oldvar =="4age_y_int", 4, 
                      ifelse(oldvar=="5age_y_int",5,
                        ifelse(oldvar=="6age_y_int",6,
                          ifelse(oldvar=="7age_y_int",7,
                            ifelse(oldvar=="8age_y_int",8,
                              ifelse(oldvar=="9age_y_int",9,0)))))))%>%
      dplyr::select(-"oldvar",-contains("BMI"))

bmi_temp_2 <- bmi_data%>%
  tidyr::gather(oldvar, BMI, `4BMI`:`9BMI`, factor_key=TRUE)%>%
    mutate("wave" = ifelse(oldvar =="4BMI", 4, 
                      ifelse(oldvar=="5BMI",5,
                        ifelse(oldvar=="6BMI",6,
                          ifelse(oldvar=="7BMI",7,
                            ifelse(oldvar=="8BMI",8,
                              ifelse(oldvar=="9BMI",9,0)))))))%>%
      dplyr::select(-"oldvar",-contains("age_y_int"))

bmi_long <- merge(bmi_temp, bmi_temp_2, by=c("HHIDPN", "wave") )

p <- ggplot(bmi_long, aes(x=age, y=BMI, group=HHIDPN))
p+geom_point()+geom_line()

#----CESD trajectory ----
cesd_data <- hrs_samp %>% 
  dplyr::select("HHIDPN", contains("age_y_int"), contains("cesd"))
cesd_temp <- cesd_data%>%
  tidyr::gather(oldvar, age, `4age_y_int`:`9age_y_int`, factor_key=TRUE )%>%
  mutate("wave" = ifelse(oldvar =="4age_y_int", 4, 
                   ifelse(oldvar=="5age_y_int",5,
                     ifelse(oldvar=="6age_y_int",6,
                       ifelse(oldvar=="7age_y_int",7,
                         ifelse(oldvar=="8age_y_int",8,
                           ifelse(oldvar=="9age_y_int",9,0)))))))%>%
  dplyr::select(-"oldvar",-contains("cesd"))

cesd_temp_2 <- cesd_data%>%
  tidyr::gather(oldvar, cesd, `r4cesd`:`r9cesd`, factor_key=TRUE)%>%
  mutate("wave" = ifelse(oldvar =="r4cesd", 4, 
                  ifelse(oldvar=="r5cesd",5,
                    ifelse(oldvar=="r6cesd",6,
                      ifelse(oldvar=="r7cesd",7,
                        ifelse(oldvar=="r8cesd",8,
                          ifelse(oldvar=="r9cesd",9,0)))))))%>%
  dplyr::select(-"oldvar",-contains("age_y_int"))

cesd_long <- merge(cesd_temp, cesd_temp_2, by=c("HHIDPN", "wave") )

p_cesd <- ggplot(cesd_long, aes(x=age, y=cesd, group=HHIDPN))
p_cesd+geom_point()+geom_line()


#---- test ----
# bmi_test <- data.frame(HHIDPN = c(1,2,3,4,5), `4BMI` = c(15,14,16,17,18), 
#                        `5BMI` = c(20,21,22,25,27),`6BMI` = c(1,2,3,4,5),
#                        `7BMI` = c(6,7,8,9,10), `8BMI` = c(11,12,13,14,110),
#                        `9BMI` = c(31,32,33,34,35),
#                        `4age` = c(55,57,58,50,70), `5age` = c(57,59,60,52,72),
#                        `6age` = c(59,61,62,54,74),`7age` = c(61,63,64,56,76),
#                        `8age` = c(63,65,66,58,78),`9age` = c(65,67,68,60,80))
# bmi_temp <- bmi_test%>%
#   tidyr::gather(oldvar, age, X4age:X9age, factor_key=TRUE )%>%
#     mutate("wave" = ifelse(oldvar =="X4age", 4, 
#                       ifelse(oldvar=="X5age",5,
#                         ifelse(oldvar=="X6age",6,
#                          ifelse(oldvar=="X7age",7,
#                            ifelse(oldvar=="X8age",8,
#                              ifelse(oldvar=="X9age",9,0)))))))%>%
#       dplyr::select(-"oldvar",-contains("BMI"))
# 
# bmi_temp_2 <- bmi_test%>%
#   tidyr::gather(oldvar, BMI, X4BMI:X9BMI, factor_key=TRUE)%>%
#     mutate("wave" = ifelse(oldvar =="X4BMI", 4, 
#                       ifelse(oldvar=="X5BMI",5,
#                        ifelse(oldvar=="X6BMI",6,
#                         ifelse(oldvar=="X7BMI",7,
#                          ifelse(oldvar=="X8BMI",8,
#                           ifelse(oldvar=="X9BMI",9,0)))))))%>%
#     dplyr::select(-"oldvar",-contains("age"))
# 
# bmi_long <- merge(bmi_temp, bmi_temp_2, by=c("HHIDPN", "wave") )
# 
# p <- ggplot(bmi_long, aes(x=age, y=BMI, group=HHIDPN))+geom_point()+geom_line()
# p
