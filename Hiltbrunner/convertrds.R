# Load the RDS file
rds_data <- readRDS("Data/sce-clustering-allsamples-original_final.rds")

# Iterate over each element
for (i in seq_along(rds_data)) {
  # Generate a file name for each element
  csv_file <- paste0("element_", i, ".csv")
  
  # Check the class of the element
  if (inherits(rds_data[[i]], "data.frame")) {
    # If it's a data frame, write it to CSV
    write.csv(rds_data[[i]], csv_file, row.names = FALSE)
  } else if (inherits(rds_data[[i]], "list")) {
    # If it's a list, convert it to data frame and then write it to CSV
    df <- as.data.frame(do.call(rbind, rds_data[[i]]))
    write.csv(df, csv_file, row.names = FALSE)
  } else {
    # If it's neither a data frame nor a list, write it to a CSV with a single row
    write.csv(t(data.frame(rds_data[[i]])), csv_file, row.names = FALSE)
  }
}

