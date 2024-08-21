# Title: utils.R
# Description: Miscellaneous useful functions

################################################################################
# Imports
################################################################################
#XXX

################################################################################
# Set up function
################################################################################
setup <- function(){
  # Set working directory
  setwd('/Users/karansgarg/Documents/MSc Statistics/STAT0034/Project/git/stat0034')
  
  # Load relevant scripts and import packages needed for general analysis
  source('./src/algorithms.R') #Import functions for HMC algorithms
  source('./src/targets.R') # Import functions for testing different targets
  source('./src/diagnostics.R') # Import functions to perform diagnostics
  library(tidyverse) # General purpose data tools (mainly ggplot2)
}

seed <- function(){
  # Set desired seed
  set.seed(06031953)
}

################################################################################
# Euclidean norm (for assessing tail behaviour)
################################################################################
norm <- function(x){return(sqrt(sum(x^2)))}

################################################################################
# Import data
################################################################################

import_data <- function(folder_path="./out/samples/"){
  
  # Get a list of all CSV files in the folder
  csv_files <- list.files(path = folder_path, pattern = "*.csv", full.names = TRUE)
  
  # Initialize an empty list to store dataframes
  data_list <- list()
  
  # Loop over the files and read them into dataframes
  for (file in csv_files) {
    
    # Extract the filename without extension
    file_name <- tools::file_path_sans_ext(basename(file))
    
    # Read the CSV file into a dataframe
    data <- read_csv(file) %>% data.frame()
    
    # Store the dataframe in the list with the filename as the list item name
    data_list[[file_name]] <- data
  }
  
  return(data_list)
}