# This function is to read in results with colnames set.
# Created by Crystal


read_results <- function(paths){
  data.table::fread(paths, fill = TRUE) %>% na.omit() %>%
    set_colnames(c("Exposure", "Beta", "SE", "LCI", "UCI", "Method",
                   "Percent", "Mechanism", "Truth Capture", "Time", "Seed"))
}

# test <- data.table::fread(main_paths[2], fill = TRUE) %>%
#   set_colnames(c("Exposure", "Beta", "SE", "LCI", "UCI", "Method",
#                  "Percent", "Mechanism", "Truth Capture", "Time")) %>% na.omit()