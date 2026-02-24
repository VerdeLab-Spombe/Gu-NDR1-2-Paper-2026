rm(list=ls())
library(tidyverse)
library(splines)
library(boot)

# Function to process a single CSV file
process_csv <- function(file_path) {
  # Read the CSV file
  data <- read.csv(file_path)
  
  # Determine the number of columns
  num_cols <- ncol(data)
  
  # Generate X axis values based on the number of columns
  x_values <- seq(-200, 200, length.out = num_cols)
  
  # Rename columns to match expected naming convention
  colnames(data) <- paste0('Lag_', seq(-20, 20, length.out = num_cols))
  
  # Initialize data frames to store results
  fitted_curves <- data.frame(matrix(ncol = num_cols, nrow = 0))
  colnames(fitted_curves) <- colnames(data)
  ci_95 <- data.frame(matrix(ncol = num_cols, nrow = 0))
  colnames(ci_95) <- colnames(data)
  
  # Aggregate all rows to get a single y_values
  y_values <- colMeans(data, na.rm = TRUE)
  
  # Remove missing or infinite values
  valid_idx <- which(!is.na(y_values) & is.finite(y_values))
  y_values <- y_values[valid_idx]
  x_values_valid <- x_values[valid_idx]
  
  # Debugging prints
  cat("Processing file:", file_path, "\n")
  cat("Number of valid points:", length(valid_idx), "\n")
  cat("x_values_valid length:", length(x_values_valid), "\n")
  cat("y_values length:", length(y_values), "\n")
  
  # Only fit the spline if there are enough valid points
  if (length(y_values) > 3) {
    # Fit a smoothing spline
    spline_fit <- smooth.spline(x_values_valid, y_values)
    
    # Perform non-parametric bootstrapping to estimate confidence intervals
    num_bootstrap <- 1000
    bootstrap_samples <- matrix(NA, nrow = num_bootstrap, ncol = length(x_values))
    
    for (b in 1:num_bootstrap) {
      resample_idx <- sample(length(y_values), length(y_values), replace = TRUE)
      bootstrap_samples[b, valid_idx] <- predict(spline_fit, x_values_valid[resample_idx])$y
    }
    
    # Calculate the 95% confidence interval
    ci_lower <- apply(bootstrap_samples, 2, quantile, probs = 0.025)
    ci_upper <- apply(bootstrap_samples, 2, quantile, probs = 0.975)
    
    # Store the results
    fitted_curve <- predict(spline_fit, x_values)$y
    fitted_curves <- rbind(fitted_curves, fitted_curve)
    ci_95 <- rbind(ci_95, ci_lower, ci_upper)
  } else {
    cat("Skipping file:", file_path, "due to insufficient valid points.\n")
  }
  
  return(list(fitted_curves = fitted_curves, ci_95 = ci_95))
}

# Main function to process all CSV files in a directory
process_all_csv_files <- function(directory_path) {
  # Get a list of all CSV files in the directory
  files <- list.files(directory_path, pattern = "\\.csv$", full.names = TRUE)
  
  # Process each file
  for (file_path in files) {
    processed_data <- process_csv(file_path)
    
    # Save the fitted curves and 95% CI data to separate files
    fitted_curves_file <- sub(".csv$", "_fitted_curves.csv", file_path)
    ci_95_file <- sub(".csv$", "_ci_95.csv", file_path)
    
    write.csv(processed_data$fitted_curves, fitted_curves_file, row.names = FALSE)
    write.csv(processed_data$ci_95, ci_95_file, row.names = FALSE)
  }
}

# Example usage:
# process_all_csv_files('path/to/your/csv/files')


# Example usage:
getwd()
process_all_csv_files("C:/Users/Lab_6130A/Desktop/New version-manuscript-02182026/Codes for paper/5-Data extract-Fig 1 D-J-Fig S2/R code_needed for data analysis/example_2/Rawdata")

# Function to read the second row of each CSV file
read_second_row <- function(file) {
  df <- read.csv(file, header = F)
  return(df[2, ])
}

# List all CSV files in the current working directory with the name pattern "fitted_curves"
csv_files <- list.files(pattern = "fitted_curves.*\\.csv")

dir.create("output")
# Read the second row from each CSV file and combine into a single DataFrame
data_list <- lapply(csv_files, read_second_row)
combined_df <- bind_rows(data_list)

# Ensure all columns are numeric and handle NA values
combined_df <- combined_df %>%
  mutate(across(everything(), ~ as.numeric(as.character(.)))) %>%
  drop_na()

# Define x-values explicitly
x_values <- seq(-200, 200, by = 10)

# Check if the number of x_values matches the number of columns in combined_df
if (length(x_values) != ncol(combined_df)) {
  stop("The number of x-values does not match the number of columns in the combined data.")
}

# Check if combined_df is not empty before proceeding
if (nrow(combined_df) > 0) {
  
  # Calculate Mean, SD, SEM for each column
  summary_stats <- combined_df %>%
    summarise(across(everything(), list(
      Mean = ~ mean(., na.rm = TRUE),
      SD = ~ sd(., na.rm = TRUE),
      SEM = ~ sd(., na.rm = TRUE) / sqrt(n())
    )))
  
  # Function to perform bootstrapping and calculate confidence intervals
  bootstrap_ci <- function(data, indices) {
    d <- data[indices, ]
    return(colMeans(d, na.rm = TRUE))
  }
  
  # Perform non-parametric bootstrapping to estimate CI95
  boot_results <- boot(
    data = as.matrix(combined_df),
    statistic = bootstrap_ci,
    R = 1000
  )
  
  # Calculate 95% confidence intervals for each column
  ci95_list <- lapply(1:ncol(combined_df), function(i) {
    ci <- boot.ci(boot_results, type = "perc", index = i)
    return(c(ci$percent[4], ci$percent[5]))
  })
  
  # Combine CI95 results into a DataFrame
  ci95_df <- do.call(rbind, ci95_list) %>%
    as.data.frame() %>%
    setNames(c("CI95_lower", "CI95_upper"))
  
  # Extract Means, SD, SEM
  means <- summary_stats %>% select(ends_with("Mean"))
  sds <- summary_stats %>% select(ends_with("SD"))
  sems <- summary_stats %>% select(ends_with("SEM"))
  
  # Save Means, SDs, SEMs, CI95 to separate CSV files
  write.csv(means, "output/fitted_curves_means.csv", row.names = FALSE)
  write.csv(sds, "output/ffitted_curves_sds.csv", row.names = FALSE)
  write.csv(sems, "output/ffitted_curves_sems.csv", row.names = FALSE)
  write.csv(ci95_df, "output/ffitted_curves_ci95.csv", row.names = FALSE)
  
  # Save the raw data to a separate CSV file
  write.csv(combined_df, "output/ffitted_curves_raw_data.csv", row.names = FALSE)
  
  # Fit a new curve using smoothing splines
  mean_values <- colMeans(combined_df, na.rm = TRUE)
  
  fit <- smooth.spline(x_values, mean_values)
  
  # Save the fitted curve data to a CSV file
  fitted_curve_df <- data.frame(x = fit$x, y = fit$y)
  write.csv(fitted_curve_df, "output/ffitted_curve.csv", row.names = FALSE)
  
  # Print the final results
  print("Raw data saved to fitted_curves_raw_data.csv")
  print("Means saved to fitted_curves_means.csv")
  print("SDs saved to fitted_curves_sds.csv")
  print("SEMs saved to fitted_curves_sems.csv")
  print("CI95 saved to fitted_curves_ci95.csv")
  print("Fitted curve data saved to fitted_curve.csv")
  
} else {
  print("No valid data found in the provided CSV files.")
}



















