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
#---- Obtain a random subset ----
set.seed(62283)
sub_hrs_samp <- dplyr::sample_n(hrs_samp, 1000)

#---- Shuffle the order of the colors----
require(scales)
n <- length(levels(as.factor(sub_hrs_samp$HHIDPN)))# number of colors
cols <- hue_pal(h = c(0, 360) + 15, 
                c = 100, l = 65, 
                h.start = 0, direction = 1)(n)[order(sample(1:n, n))] 

#----BMI trajectory ----
bmi_data <- sub_hrs_samp %>% 
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
bmi_long$HHIDPN <- as.factor(bmi_long$HHIDPN)
p <- ggplot(bmi_long, aes(x=age, y=BMI, group=HHIDPN, color = HHIDPN))
p+geom_point()+geom_line()+theme(legend.position = "none")+
  scale_color_manual(values = cols)

#----CESD trajectory ----
cesd_data <- sub_hrs_samp %>% 
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

p_cesd <- ggplot(cesd_long, aes(x=age, y=cesd, group=HHIDPN, color = as.factor(HHIDPN)))
p_cesd+geom_point()+geom_line()+theme(legend.position = "none")+
  scale_color_manual(values = cols) 


