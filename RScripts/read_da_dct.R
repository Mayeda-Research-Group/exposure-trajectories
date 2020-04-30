#---- Function Overview ----
#This function reads in and formats data using a data file (.da) 
#and a data dictionary file (.dct). 
#This can be used to read in HRS data files such as the tracker file and 
#senstive data files

#Package dependencies: readr

#Inputs: 
# data_path: path to the .da file
# dict_path: path to the .dct file
# HHIDPN: defaults to TRUE; merges the HHID and PN variables

#Ouput: formatted dataframe 

#---- Load package ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("readr")

#---- Function ----
read_da_dct <- function(data_path, dict_path, HHIDPN = TRUE){
  
  # Read the dictionary file
  df_dict <- read.table(dict_path, skip = 2, fill = TRUE,
                        stringsAsFactors = FALSE)
  
  #Set column names for dictionary dataframe
  colnames(df_dict) <- c("col.num", "col.type", "col.name", "col.width",
                         "col.lbl")
  
  #Remove last row which only contains a closing}
  df_dict <- df_dict[-nrow(df_dict), ]
  
  #Extract numeric value from column width field
  df_dict$col.width <- as.integer(sapply(df_dict$col.width, gsub,
                                         pattern = "[^0-9\\.]",
                                         replacement = ""))
  
  #Convert column types to format to be used with read_fwf function
  df_dict$col.type <-
    sapply(df_dict$col.type,
           function(x) ifelse(x %in% c("int","byte","long"), "i",
                              ifelse(x == "float", "n",
                                     ifelse(x == "double", "d", "c"))))
  
  #Read the data file into a dataframe
  data <- readr::read_fwf(file = data_path,
                          readr::fwf_widths(widths = df_dict$col.width,
                                            col_names = df_dict$col.name),
                          col_types = paste(df_dict$col.type, collapse = ""))
  
  # Add column labels to headers
  attributes(data)$variable.labels <- df_dict$col.lbl
  
  #Merge HHID and PN
  if(HHIDPN == TRUE){
    data %<>% as.data.frame() %>% unite("HHIDPN", c("HHID", "PN"), sep = "")
  } else{
    data %<>% as.data.frame()
  }
  
  return(data)
}

