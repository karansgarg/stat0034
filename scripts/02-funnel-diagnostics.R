# Title: 02-funnel-diagnostics.R
# Description: Analysing results of sampling: creating summary statistics (bias)

################################################################################
# Setup
################################################################################
# Load utils functions
source('/Users/karansgarg/Documents/MSc Statistics/STAT0034/Project/git/stat0034/src/utils.R')
setup() # Import relevant scripts

################################################################################
# Import data
################################################################################
# Define the folder path containing sampling data
data_folder <- "./out/samples/"

# Get a list of all CSV files in the folder
csv_files <- list.files(path = data_folder, pattern = "*.csv", full.names = TRUE)

# Initialize an empty list to store dataframes
sampling_data <- list()

# Loop over the files and read them into dataframes
for (file in csv_files) {
  
  # Extract the filename without extension
  file_name <- tools::file_path_sans_ext(basename(file))
  
  # Read the CSV file into a dataframe
  data <- read_csv(file) %>% data.frame()
  
  # Store the dataframe in the list with the filename as the list item name
  sampling_data[[file_name]] <- data
}

# Import NUTS/Stan data separately
#stan_samples <- read_csv(file='./out/_old/funnel_exploration/samples/stan.csv')
#stan_samples <- stan_samples$`stan_y$y`

################################################################################
# Extract y values (first dimension of funnel samples)
################################################################################
# The focus of the diagnostics will be whether the range of y values were appropriately
# sampled. The neck of the funnel corresponds to where y takes negative values, so
# we would expect to see worse performance in sampling negative y values if there
# were issues with the sampler.

y_data <- lapply(sampling_data, function(sampled_data_df){return(sampled_data_df[,1])})

################################################################################
# Compute rolling bias
################################################################################
print(Sys.time())
y_rolling_bias <- lapply(y_data, rolling_bias)
print(Sys.time())
y_rolling_bias_sq <- lapply(y_data, rolling_bias_sq)
print(Sys.time())
# stan_rolling_bias <- rolling_bias(stan_samples)
# print(Sys.time())
# stan_rolling_bias_sq <- rolling_bias_sq(stan_samples)
# print(Sys.time())

# Export rolling bias data
# Directory to save the CSV files
output_dir <- "./out/diagnostics/rolling_bias/"
output_dir_sq <- "./out/diagnostics/rolling_bias_sq/"

# Loop over the list and write each dataframe to a CSV file
for(name in names(y_rolling_bias)){
  
  # Create the full path for the output file and write dataframe to CSV
  output_file <- file.path(output_dir, paste0(name, ".csv"))
  write.csv(y_rolling_bias[[name]], file = output_file, row.names = FALSE)
  
  # Create the full path for the output file and write dataframe to CSV (for sq)
  output_file_sq <- file.path(output_dir_sq, paste0(name, ".csv"))
  write.csv(y_rolling_bias_sq[[name]], file = output_file_sq, row.names = FALSE)
}

# output_file <- file.path(output_dir, paste0("stan", ".csv"))
# write.csv(stan_rolling_bias, file = output_file, row.names = FALSE)
# 
# output_file_sq <- file.path(output_dir_sq, paste0("stan", ".csv"))
# write.csv(stan_rolling_bias_sq, file = output_file_sq, row.names = FALSE)