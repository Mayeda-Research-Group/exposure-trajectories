#This function reads in and formats data so that we can reads in the 
#HRS tracker file and HRS biomarker data

read_data <- function(data_path, dict_path){
  
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
  
  data %<>% as.data.frame() %>% unite("HHIDPN", c("HHID", "PN"), sep = "")
  
  return(data)
}

