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
  # Load relevant scripts and import packages needed for general analysis
  source('./src/algorithms.R') #Import functions for HMC algorithms
  source('./src/targets.R') # Import functions for testing different targets
  source('./src/diagnostics.R') # Import functions to perform diagnostics
  library(tidyverse)
}

seed <- function(){
  # Set desired seed
  set.seed(06031953)
}